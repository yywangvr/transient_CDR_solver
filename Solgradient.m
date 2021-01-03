function  [ugrad]= Solgradient(X,T,u,referenceElement )

% reference element information
nen = referenceElement.nen; 
ngaus = referenceElement.ngaus; 
wgp = referenceElement.GaussWeights; 
N = referenceElement.N; 
Nxi = referenceElement.Nxi; 
Neta = referenceElement.Neta;

% Number of nodes and elements
[nPt,nDim] = size(X); 

nElem = size(T,1); 

ugrad=zeros(nPt,nDim);

M = zeros(nPt);
GradNx = zeros(nPt,1);
GradNy = zeros(nPt,1);



for ielem =1: nElem
    
    Te = T(ielem,:); 
    Xe = X(Te,:);
    ue = u(Te,:);
    
    

    [Me,GradNxe,GradNye]=EleGrad(Xe,ue,nen,ngaus,wgp,N,Nxi,Neta);
    
     M(Te,Te) = M(Te,Te) + Me;
     GradNx(Te) = GradNx(Te)+GradNxe;
     GradNy(Te) = GradNy(Te)+GradNye;
    
    
end

%[L,U]=lu(diag(sum(M,2)));
[L,U]=lu(M);
ugrad(:,1)= U\(L\GradNx);
ugrad(:,2)= U\(L\GradNy);


end



function [Me,GradNx,GradNy]=EleGrad(Xe,ue,nen,ngaus,wgp,N,Nxi,Neta)

Me = zeros(nen);
GradNx= zeros(nen,1);
GradNy= zeros(nen,1);



for ig=1:ngaus
    
    N_ig = N(ig,:);
    Nxi_ig = Nxi(ig,:);
    Neta_ig = Neta(ig,:);
    Jacob = [Nxi_ig*(Xe(:,1))	Nxi_ig*(Xe(:,2))
        Neta_ig*(Xe(:,1))	Neta_ig*(Xe(:,2))];
    dvolu = wgp(ig)*det(Jacob);
   
    
    res = Jacob\[Nxi_ig;Neta_ig];
    Nx = res(1,:);
    Ny = res(2,:);
    
    Me = Me + N_ig'*N_ig*dvolu;
    GradNx= GradNx + N_ig'*(Nx*ue)* dvolu;
    GradNy= GradNy + N_ig'*(Ny*ue)* dvolu;
    
    
end  
    
end



