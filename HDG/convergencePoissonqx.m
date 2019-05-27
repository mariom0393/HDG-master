% Academic 2D HDG code for solving the Poisson equation with Dirichlet 
% boundary conditions.
%
% First version by Sonia Fernandez-Mendez LaCaN UPC-BarcelonaTech 2016
% Second version by Ruben Sevilla, Swansea University 2017
% Third version by Matteo Giacomini, UPC-BarcelonaTech 2019
%
% www.lacan.upc.edu
%
% Main data variables:
%  X: nodal coordinates
%  T: mesh connectivitity matrix
%  F: faces (here sides) for each element 
%    (Faces are numbered so that interior faces are first)
%  elemInfo: element type
%  infoFaces.intFaces: [elem1 face1 elem2 face2 rotation] for each face
%    (Each face is assumed to have the orientation given by the 
%    first element, it is flipped when seen from the second element)
%  infoFaces.extFaces: [elem1 face1] for each exterior face
%  referenceElement: integration points and shape functions (volume and
%    boundary sides)
%

clearvars
close all
setpath
global mu
mu=1;

% Degree of polynomials to be studied
pDeg = [1 2 3 4];
% Number of mesh for each polynomial degree
hRef = [3 3 3 3];

% bool compute VS load
computeError = 1;
% Error storage
uErr = zeros(sum(hRef,2),1);
uErrStar = zeros(sum(hRef,2),1);
uErrqx = zeros(sum(hRef,2),1);
hSize = zeros(sum(hRef,2),1);
% Index to store the errors
iErr = 1;

col = {'b-o', 'r-s', 'g-d', 'k-^', 'm->'};
colStar = {'b--o', 'r--s', 'g--d', 'k--^', 'm-->'};
leg = {};


%% Compute the error for different configurations
% Loop over polynomial degree
for iDeg=1:length(pDeg)
    
    degree=pDeg(iDeg); 
    nOfMesh = hRef(iDeg);
    
    % Loop over computational meshes
    for iMesh=1:nOfMesh

        fprintf(' == Mesh #%d - Polynomials of degree #%d == \n',iMesh,degree)
        
        % Load computational mesh 
        load(sprintf('mesh%d_P%d.mat',iMesh,degree));
        nOfElements = size(T,1);
        hMax = computeMeshSizeTri2D(T,X,nOfElements);
        
        % Setup physical coefficients
        % Viscosity parameter
        muElem = mu*ones(nOfElements,1);

        % Stabilization parameter
        tau = ones(nOfElements,3); %Uniform tau parameter (all faces)
        %tau(:,2:3) = 0; %Non-null only for 1st face in each element
        
        if computeError==1
  
            % HDG preprocess
            [F, infoFaces] = hdg_preprocess(T,X);
            nOfFaces = max(max(F)); 

            % Computation
            % Loop in elements
            referenceElement=createReferenceElementTri(degree);
            disp('Loop in elements...')
            [K,f, QQ, UU, Qf, Uf] = hdgMatrixPoisson(muElem,X,T,F,referenceElement,infoFaces,tau);

            %Dirichlet BC
            %Dirichlet face nodal coordinates
            nOfFaceNodes = degree+1; 
            nOfInteriorFaces = size(infoFaces.intFaces,1);
            nOfExteriorFaces_D = size(infoFaces.extFaces_D,1);
            nOfExteriorFaces_N = size(infoFaces.extFaces_N,1);
            nOfExteriorFaces_R = size(infoFaces.extFaces_R,1);

            uDirichlet = computeProjectionFaces(@analyticalPoisson,infoFaces.extFaces_D,X,T,referenceElement);
            dofDirichlet= (nOfInteriorFaces+nOfExteriorFaces_N+nOfExteriorFaces_R)*nOfFaceNodes + (1:nOfExteriorFaces_D*nOfFaceNodes);
            dofUnknown = 1:(nOfInteriorFaces*nOfFaceNodes+nOfExteriorFaces_N*nOfFaceNodes+nOfExteriorFaces_R*nOfFaceNodes);

            % System reduction (Dirichlet faces  are set to prescribed value)
            f = f(dofUnknown)-K(dofUnknown,dofDirichlet)*uDirichlet;
            K = K(dofUnknown,dofUnknown);

            % Face solution
            disp('Solving linear system...')
            lambda = K\f; 
            uhat = [lambda(1:(nOfInteriorFaces+nOfExteriorFaces_N+nOfExteriorFaces_R)*nOfFaceNodes); uDirichlet];

            % Elemental solution
            disp('Calculating element by element solution...')
            [u,q]=computeElementsSolution(uhat,UU,QQ,Uf,Qf,F);
            
            figure(1),clf,
            plotMesh(X,T);

            figure(3),clf
            plotDiscontinuosSolution(X,T,q(:,1),referenceElement,20)
            colorbar, title('HDG solution: q_x')
            set(gca,'fontsize', 18);

            % Errors
            % analytical solution
            u_ex = @analyticalPoisson;
            qx_ex = @analyticalPoissonqx;

            %Error for the HDG solution
            Error=computeL2Norm(referenceElement,X,T,u,u_ex);
            fprintf('Error HDG = %e\n',Error);
            
            %Error for the HDG solution qx
            Errorqx=computeL2Norm(referenceElement,X,T,q(:,1),qx_ex);
            fprintf('Error HDG qx= %e\n',Error);
            

            % Store solution
            solFile = sprintf('poissonDir_mesh%d_P%d_tau%d.mat',iMesh,degree,tau(1));
            save(solFile);
            
            uErrqx(iErr) = Errorqx;
            hSize(iErr) = hMax;
            iErr = iErr + 1;
        else
            
            solFile = sprintf('poissonDir_mesh%d_P%d_tau%d.mat',iMesh,degree,tau(1));
            load(solFile);
            computeError = 0; % do not change this variable with the loaded value
            uErrqx(iErr) = Errorqx;
            hSize(iErr) = hMax;
            iErr = iErr + 1;
        end
    end
end


%% Postprocess
for iDeg=1:length(pDeg)
    
    degree=pDeg(iDeg); 
    nOfMesh = hRef(iDeg);
    fprintf(' == Convergence for polynomials of degree #%d == \n',degree);
    oldMesh = hRef(1:iDeg-1);
    meshIni = sum(oldMesh,2) + 1;
    meshFin = meshIni + nOfMesh - 1;

    % Rate of convergence
    denom = log10(hSize(meshIni+1:meshFin,:)) - log10(hSize(meshIni:meshFin-1,:));
    slopesL2 = ( log10(uErrqx(meshIni+1:meshFin)) - log10(uErrqx(meshIni:meshFin-1)) )./denom(:,1);
    
    % Plot
    figure(100), hold on
    plot(log10(hSize(meshIni:meshFin,1)), log10(uErrqx(meshIni:meshFin)), col{iDeg}, 'LineWidth', 2, 'MarkerSize', 8)
    
    leg{iDeg} = sprintf('{\\boldmath{$q_x$}}, $k$=%d',degree);
end

box on
grid
xlabel('log$_{10}(h)$','Interpreter','latex','FontName','cmr12')
ylabel('log$_{10}(||E||_{L^2(\Omega)})$','Interpreter','latex','FontName','cmr12')
set(gca,'FontSize',22,'FontName','cmr12')
h = legend(leg, 'Location', 'Southeast','FontName','cmr12');
set(h,'Interpreter','latex');