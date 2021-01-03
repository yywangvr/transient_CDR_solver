function [K,f] = FEM_system(X,T,referenceElement,example)
% [K,f] = FEM_system(X,T,referenceElement,example)
% Matrix K and r.h.s vector f obtained after discretizong a 2D convection-diffusion equation
%
% X:            nodal coordinates
% T:            connectivities (elements)
% referenceElement: reference element properties (quadrature, shape functions...)
% example: example properties


velo = example.velo;
a = example.a;
nu = example.nu;
sigma=example.sigma;
method = example.method; 

Te = T(1,:); Xe = X(Te,:); 
hx = max(Xe(:,1)) - min(Xe(:,1)); 
hy = max(Xe(:,2)) - min(Xe(:,2)); 
h = (hx+hy)/2; 


nen = referenceElement.nen;
ngaus = referenceElement.ngaus;
wgp = referenceElement.GaussWeights;
N = referenceElement.N;
Nxi = referenceElement.Nxi;
Neta = referenceElement.Neta;

% Number of elements and number of nodes in the mesh
nElem = size(T,1);
nPt = size(X,1);

K = zeros(nPt,nPt);
f = zeros(nPt,1);

if method == 0
    % Galerkin
    tau = 0;

else
    % SUPG GLS
    Pe = a*h/(2*nu);
    tau_p = h*(1 + 9/Pe^2+((h*sigma)/(2*a))^2)^(-1/2)/(2*a);
    disp(strcat('Recommended stabilization parameter = ',num2str(tau_p)));
    tau = cinput('Stabilization parameter',tau_p);
    if isempty(tau)
        tau = tau_p;
    end

end

% Loop on elements
for ielem = 1:nElem
    % Te: global number of the nodes in the current element
    Te = T(ielem,:);
    % Xe: coordenates of the nodes in the current element
    Xe = X(Te,:);
    % Element matrices
    [Ke,fe] = EleMat(Xe,nen,ngaus,wgp,N,Nxi,Neta,method,tau,velo,nu,sigma);
    % Assemble the element matrices
    K(Te,Te) = K(Te,Te) + Ke;
    f(Te) = f(Te) + fe;
end


function [Ke,fe] = EleMat(Xe,nen,ngaus,wgp,N,Nxi,Neta,method,tau,velo,nu,sigma)
%

ax = velo(1);
ay = velo(2);

Ke = zeros(nen,nen);
fe = zeros(nen,1);


% Loop on Gauss points (computation of integrals on the current element)
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
    
    if method == 0
        % Galerkin
        Ke = Ke + (nu*(Nx'*Nx+Ny'*Ny) + N_ig'*(ax*Nx+ay*Ny)+sigma*(N_ig'*N_ig))*dvolu;
        aux = N_ig*Xe; 
        f_ig = SourceTerm(aux);
        fe = fe + N_ig'*(f_ig*dvolu);
    elseif method == 1
        % GLS method
        Ke = Ke + (nu*(Nx'*Nx+Ny'*Ny) + N_ig'*(ax*Nx+ay*Ny) +sigma*(N_ig'*N_ig)+...
            tau*(ax*Nx+ay*Ny+sigma*N_ig)'*(ax*Nx+ay*Ny+sigma*N_ig))*dvolu;
        aux = N_ig*Xe; 
        f_ig = SourceTerm(aux);
        fe = fe + (N_ig+tau*(ax*Nx+ay*Ny+sigma*N_ig))'*(f_ig*dvolu);
        
    else
        % SUPG
        Ke = Ke + (nu*(Nx'*Nx+Ny'*Ny) + N_ig'*(ax*Nx+ay*Ny) +sigma*(N_ig'*N_ig)+...
            tau*(ax*Nx+ay*Ny)'*(ax*Nx+ay*Ny+sigma*N_ig))*dvolu;
        aux = N_ig*Xe; 
        f_ig = SourceTerm(aux);
        fe = fe + (N_ig+tau*(ax*Nx+ay*Ny))'*(f_ig*dvolu);
    end
end

