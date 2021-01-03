function [M,K,C,Sm,Sk,f] = FEM_matrices(X,T,Conv,u_aux,referenceElement,ugrad,para) 


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
Sm=zeros(nPt);
Sk=zeros(nPt);
f = zeros(nPt,1);

if strcmp(para.method,'OSS')
    
    Cs=zeros(nPt);
    Cr= zeros(nPt);
    Cl= zeros(nPt);
    Ks= zeros(nPt);
    L= zeros(nPt);
end





% Loop on elements
for ielem=1:nElem
    Te = T(ielem,:); 
    Xe = X(Te,:); 
    ugrade=ugrad(Te,:);
    ue = u_aux(Te,:);
    Conve = Conv(Te,:);
   
    [Me,Ke,Ce,Sme,Ske,Cse,Cle,Cre,Kse,Le,fe] = EleMat(Xe,Conve,ue,ugrade,nen,ngaus,wgp,N,Nxi,Neta,para);
        
    
    % Assembly
    M(Te,Te) = M(Te,Te) + Me; 
    K(Te,Te) = K(Te,Te) + Ke; 
    C(Te,Te) = C(Te,Te) + Ce;
    if strcmp(para.method,'SUPG')
    Sm(Te,Te)=Sm(Te,Te)+Sme;
    Sk(Te,Te)=Sk(Te,Te)+Ske;
    end
    f(Te)=f(Te)+fe;
    
    if strcmp(para.method,'OSS')
        Cs(Te,Te)= Cs(Te,Te)+Cse;
        Cl(Te,Te)= Cl(Te,Te)+Cle;
        Cr(Te,Te)= Cr(Te,Te)+Cre;
        Ks(Te,Te)= Ks(Te,Te)+Kse;
        L(Te,Te)= L(Te,Te)+Le;
    end
   
end

M = sparse(M);
K = sparse(K);
C = sparse(C);
Sm=sparse(Sm);

f = sparse(f);


if strcmp(para.method,'OSS')
    
    aux1=M^(-1)*L;
    
    Sk=Cs-Cl*aux1;% +aux1*Cr  -aux1*Ks*aux1;
    %Sk=Cl*(eye(nPt)-aux1);
    %Sk=2*Cs-Cl*aux1-aux1*Cr;
    
    Sk=sparse(Sk);
end

end


function [Me,Ke,Ce,Sme,Ske,Cse,Cle,Cre,Kse,Le,fe] = EleMat(Xe,Conve,ue,ugrade,nen,ngaus,wgp,N,Nxi,Neta,para)
%

Me = zeros(nen);
Ke = zeros(nen);
Ce = zeros(nen);
Sme = zeros(nen);
Ske = zeros(nen);
fe=zeros(nen,1);

Cse = zeros(nen);

Cle = zeros(nen);

Cre = zeros(nen);

Kse =zeros(nen);
Le = zeros(nen);

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
    aGradN = (ax*Nx) + (ay*Ny);
    
    aux=N_ig*Xe;
    fn=SourceTerm(aux);
    fn1=SourceTerm(aux);
    
    %shock capturing

    Projectedugrade=(N_ig*ugrade)';
    
    % interpolated gradient
    grtem=res*ue;
    k = shockcapturing(para,a,Projectedugrade,grtem);
    %k=0;

   if strcmp(para.method,'SUPG')
    
    tau=staticTau(para,a);
    P=aGradN';
    P=tau*P;
    
    
    Sme=Sme+P*N_ig*dvolu;
    Ske=Ske+P*(aGradN+para.sigma*N_ig)*dvolu;
    Me = Me + N_ig'*N_ig*dvolu; 
    
    aux_nu=para.nu+k;
    Ke = Ke + (aux_nu*(Nx'*Nx+Ny'*Ny))*dvolu;
    Ce = Ce + N_ig'*aGradN*dvolu;
    fe = fe + (P+N_ig')*(para.theta*fn1+(1-para.theta)*fn)*dvolu;
    
   elseif strcmp(para.method,'OSS')
       
    tau=staticTau(para,a);
  % tau=transientTau(tau,para);
    
    % standard mass matrix
    Me = Me + N_ig'*N_ig*dvolu; 
    
    % stiffness matrix with diffusion
    Ke = Ke + ((para.nu+k)*(Nx'*Nx+Ny'*Ny))*dvolu;
    
    % standard convection matrix
    Ce = Ce + N_ig'*aGradN*dvolu;
    
    % stabilization matrix
    % sysmmetric convection-stiff matrix
    Cse = Cse +tau*(aGradN'*aGradN)*dvolu;
    
    % left convection-stiff
    Cle = Cle + tau*aGradN'* N_ig*dvolu ;
    
    % righ convection stiff
    Cre = Cre +tau*(N_ig'*aGradN)*dvolu;
    
    % stiff-stiff
    Kse = Kse + tau*(N_ig'*N_ig)*dvolu;
    
    % 
    Le = Le + N_ig'*aGradN*dvolu;
    
    
   else
       disp('No Stabilization')
   end
   
end

end





function Tau=staticTau(para,a)
    
    Tau=(2*norm(a)/para.h+(4*para.nu/(para.h)^2)+para.sigma)^(-1);
    
    %Pe = norm(a)*para.h/(2*para.nu);
    %Tau = para.h*(1 + 9/Pe^2+((para.h*para.sigma)/(2*norm(a)))^2)^(-1/2)/(2*norm(a));
    
end


function transientTau=transientTau(staticTau,para)

    transientTau=(inv(staticTau)+inv(para.theta*para.dt))^(-1) ;   

end





