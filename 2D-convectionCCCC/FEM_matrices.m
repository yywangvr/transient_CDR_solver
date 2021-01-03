function [M,K,C] = FEM_matrices(X,T,Conv,referenceElement) 


% reference element information
nen = referenceElement.nen; 
ngaus = referenceElement.ngaus; 
wgp = referenceElement.GaussWeights; 
N = referenceElement.N; 
Nxi = referenceElement.Nxi; 
Neta = referenceElement.Neta; 

% Number of nodes and elements
nPt = size(X,1); 
nElem = size(T,1); 

M = zeros(nPt);
K = zeros(nPt); 
C = zeros(nPt);

% Loop on elements
for ielem=1:nElem
    Te = T(ielem,:); 
    Xe = X(Te,:); 
    Conve = Conv(Te,:);
    [Me,Ke,Ce] = EleMat(Xe,Conve,nen,ngaus,wgp,N,Nxi,Neta);
    % Assembly
    M(Te,Te) = M(Te,Te) + Me; 
    K(Te,Te) = K(Te,Te) + Ke; 
    C(Te,Te) = C(Te,Te) + Ce; 
end



function [Me,Ke,Ce] = EleMat(Xe,Conve,nen,ngaus,wgp,N,Nxi,Neta)
%

Me = zeros(nen);
Ke = zeros(nen);
Ce = zeros(nen);


% Loop on Gauss points 
for ig = 1:ngaus
    N_ig = N(ig,:);
    Nxi_ig = Nxi(ig,:);
    Neta_ig = Neta(ig,:);
    Jacob = [Nxi_ig*(Xe(:,1))	Nxi_ig*(Xe(:,2))
             Neta_ig*(Xe(:,1))	Neta_ig*(Xe(:,2))];
    dvolu = wgp(ig)*det(Jacob);
    res = Jacob\[Nxi_ig;Neta_ig];
    Nx = res(1,:);
    Ny = res(2,:);
    a = N_ig*Conve; ax = a(1); ay = a(2); 
    
    Me = Me + N_ig'*N_ig*dvolu; 
    aGradN = (ax*Nx) + (ay*Ny);
    
    Ke = Ke + aGradN'*aGradN*dvolu;
    Ce = Ce + N_ig'*aGradN*dvolu;
end

