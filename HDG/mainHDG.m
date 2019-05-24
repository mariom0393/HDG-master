% Academic 2D HDG code for solving the Poisson equation with Neuman 
% boundary conditions.
%
% Developed by: Marito, Agus, Carlitos y Eugenito.
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
  
%% Load computational mesh 
degree=2; 
load(['mesh4_P',num2str(degree)]); 
 
nOfElements = size(T,1);

figure(1),clf,plotMesh(X,T);

%% HDG preprocess
[F, infoFaces] = hdg_preprocess(T,X);
nOfFaces = max(max(F)); 

%% Viscosity parameter
muElem = mu*ones(nOfElements,1);

%% Stabilization parameter
tau = ones(nOfElements,3); %Uniform tau parameter (all faces)
%tau(:,2:3) = 0; %Non-null only for 1st face in each element

%% Computation
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

figure(2),clf
plotDiscontinuosSolution(X,T,u,referenceElement,20)
colorbar, title('HDG solution: u')

% figure(3),clf
% plotDiscontinuosSolution(X,T,q(:,1),referenceElement,20)
% colorbar, title('HDG solution: q_x')
% 
% figure(4),clf
% plotDiscontinuosSolution(X,T,q(:,2),referenceElement,20)
% colorbar, title('HDG solution: q_y')

% figure(5),clf
% plotjump(X,T,q(:,2),referenceElement,20)
% colorbar, title('HDG solution: q_y')

figure(3),clf
plotDiscontinuosq(X,T,q,referenceElement,20)
colorbar, title('HDG solution: q')

% Local postprocess for superconvergence 
 disp('Performing local postprocess...')
 referenceElement_star = createReferenceElementTriStar(referenceElement);
 u_star = HDGpostprocess(X,T,u,-q,referenceElement_star);

%Plots postprocess solution
figure(22),clf
plotPostprocessedSolution(X,T,u_star,referenceElement_star,20);
colorbar, title('HDG solution: u*')

%% Errors

% analytical solution
u_ex = @analyticalPoisson;

%Error for the HDG solution
Error=computeL2Norm(referenceElement,X,T,u,u_ex);
fprintf('Error HDG = %e\n',Error);

% Error for the postprocessed solution
ErrorPost=computeL2NormPostprocess(referenceElement_star,X,T,u_star,u_ex);
fprintf('Error HDG postprocessed = %e\n',ErrorPost);
disp(' ')