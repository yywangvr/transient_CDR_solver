function  u = Method(para,X,T,referenceElement,nStep)
 

Conv = ComputeVelocity(X,para.velo);
figure(1); hold on;
quiver(X(:,1),X(:,2),Conv(:,1),Conv(:,2))
plot(para.dom([1,2,2,1,1]), para.dom([3,3,4,4,3]), 'k')
axis equal


[ADir1,bDir1,nDir1] = BoundaryConditions(X,para.velo);
u=zeros(size(X,1),nStep+1);
u(:,1) = para.u0;

%switch method
%    case 1 % Crank-Nicolson + Galerkin
       
        
for n = 1:nStep
    for j=1:para.Picard
        
        % gradient of the solution
        
        if j==1
            u(:,n+1)=u(:,n);
     
        end
        
        ugrad= Solgradient(X,T,u(:,n+1),referenceElement);
        
        % matrix assembling
        [M,K,C,Sm,Sk,f] = FEM_matrices(X,T,Conv,u(:,n+1),referenceElement,ugrad,para);
        
        A = M+Sm + para.theta*(C+K+para.sigma*M+Sk)*para.dt;
        b = -(C+K+para.sigma*M+Sk)*u(:,n);
        
       %A = M + para.theta*(C+K+para.sigma*M)*para.dt;
        %b = -(C+K+para.sigma*M)*u(:,n);
        
        Ktot = sparse([A ADir1';ADir1 zeros(nDir1,nDir1)]);
        
        [L,U] = lu(Ktot);
        
        ftot=[(b*para.dt + f*para.dt);bDir1];
        
        
        sol = U\(L\ftot);
        Du = sol(1:length(f));
        
        
        u(:,n+1)=u(:,n) + Du;
        
        
       
        
        
    end
    
    
end
end        
%     case 2  % R22+Galerkin
%         
%         W=(1./24.)*[7 -1;13 5];
%         w=0.5*[1;1];
%         Id=[1,0;0,1];
%         n=size(W,1);
%         npt=size(K,1);
%         A=zeros(n*npt);
%         b=zeros(n*npt);
%         f=zeros(n*npt,1);
%         
%         ADir=zeros(n*nDir1,n*npt);
%         bDir=zeros(n*nDir1,1);
%         
%         for i=1:n
%             for j=1:n
%                 A((i-1)*npt+1:i*npt,(j-1)*npt+1:j*npt) = Id(i,j)*M/dt + W(i,j)*(C+K+sigma*M);
%                 b((i-1)*npt+1:i*npt,(j-1)*npt+1:j*npt) = -w(i)*(C+K+sigma*M);
%                 ADir((i-1)*nDir1+1:i*nDir1,(j-1)*npt+1:j*npt)=Id(i,j)*ADir1;
%                 
%                 
%             end
%             f((i-1)*npt+1:i*npt)=w(i)*f1;
%             bDir((i-1)*nDir1+1:i*nDir1)=bDir1;
%         end
%         Ktot = [A ADir';ADir zeros(n*nDir1,n*nDir1)];
%         
%         f=sparse(f);
%         [L,U] = lu(Ktot);
%         L=sparse(L);
%         U=sparse(U);
%         
%         for n = 1:nStep
%             
%             temp=[b(1:npt,1:npt)*u(:,n);b(npt+1:2*npt,npt+1:2*npt)*u(:,n)] + f;
%             
%             ftot=[temp;bDir];
%             sol =U\(L\ftot);
%             Du = sol(1:length(f));
%             Du=reshape(Du,npt,2);
%             u(:,n+1) = u(:,n) + sum(Du,2);
%             
%         end
%         
% end




