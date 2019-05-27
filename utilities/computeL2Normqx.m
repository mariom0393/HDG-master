function L2Norm = computeL2Normq(referenceElement,X,T,q,q0,varargin)

% The function computeL2Norm allows to compute the L2 Norm using an
% analytical function handle. This function has to depend on x-y
% coordinates.
%
% Input:
%  referenceElement: information of the reference element
%  X: nodal coordinates
%  T: connectivity matrix
%  q: FEM solution (nodal)
%  q0: analytical function handle. This function has to have the x-y
%      coordinates matrix [xi yi] as a first input argument.
%  varargin: ordered list of input arguments needed to execute u0 but x-y
%            coordinates.
% Output:
%  L2Norm: computed L2 norm = sqrt(integral[(u-u0)^2]) over the domain defined by
%          T,X

%Number of elements and number of mesh nodes
[nOfElements,nOfElementNodes] = size(T);

%Loop in 2D elements
L2Norm = 0;
for iElem = 1:nOfElements 
    Te = T(iElem,:); 
    Xe = X(Te,:);
    ind = (iElem-1)*nOfElementNodes+1:iElem*nOfElementNodes;
    qe = q(ind);
    L2Norm = L2Norm + ElementalL2Norm(Xe,referenceElement,qe,q0,varargin{:}); 
end
L2Norm = sqrt(L2Norm);
end

%_______________________________________________________________________
function elemL2Norm = ElementalL2Norm(Xe,referenceElement,qe,q0,varargin)

    %Information of the reference element
    IPw = referenceElement.IPweights; 
    N = referenceElement.N;
    Nxi = referenceElement.Nxi;
    Neta = referenceElement.Neta;

    %Number of Gauss points
    ngauss = length(IPw);

    % x and y coordinates of the element nodes
    xe = Xe(:,1); ye = Xe(:,2);

    %Jacobian (straight faces)
    % v1 = Xe(1,:);  v2 = Xe(2,:);  v3 = Xe(3,:);
    % J = [(v2-v1)/2 ; (v3-v1)/2];
    % detJ = det(J);

    %Compute elemental L2 Norm
    elemL2Norm = 0;
    for g = 1:ngauss
    %Values at current integration point
    N_g = N(g,:);
    Nxi_g = Nxi(g,:);
    Neta_g = Neta(g,:);
    xy_g = N_g*Xe;
    q0_g = q0(xy_g,varargin{:});
    %Jacobian
    J = [Nxi_g*xe	  Nxi_g*ye   
         Neta_g*xe  Neta_g*ye];
    
    detJ=det(J);
    invJ11= J(2,2)/detJ;
    invJ12= -J(1,2)/detJ;
    invJ21= -J(2,1)/detJ;
    invJ22= J(1,1)/detJ;
     
    % xy-derivatives
    Nx = invJ11*Nxi_g + invJ12*Neta_g;
    Ny = invJ21*Nxi_g + invJ22*Neta_g;
    
    %Integration weight
    dvolu=IPw(g)*det(J);
    
    qe_g=N_g*(qe);
    %Contribution of the current integration point to the elemental L2 Norm 
    elemL2Norm = elemL2Norm + (qe_g-q0_g)^2*dvolu;
    end
end