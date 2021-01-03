function fn=neumannBC(X,nodes,neumann)


% 1D quadrature
% ngaus = 3; 
% zgp = [-sqrt(15)/5; 0; sqrt(15)/5]; 
% wgp = [5/9   8/9   5/9];
% N   = [(zgp-1).*zgp/2   1-zgp.^2     (zgp+1).*zgp/2]; 
ngaus = 2; 
zgp = [-1/sqrt(3); 1/sqrt(3)]; 
wgp = [1   1];
 N =  [(1-zgp)/2     (1+zgp)/2];  
nPt= size(X,1);

%initialize the fn by a zero matrix
fn=zeros(nPt,1);

%form a new nodes coordinate
%X0=X(nodes,:);

%form a connectivity matrix corresponding to the nodes coodinate
T0=[nodes(1:1:end-1),nodes(2:1:end)];

nElem=size(T0,1);

%Loop on elements

for ielem=1:nElem
   Te=T0(ielem,:);
   Xe=X(Te,:);
   dd=Xe(end,:)-Xe(1,:);
   h=norm(dd);
   fe=zeros(size(Te,2),1);
   %loop on gauss points
   for ig=1:ngaus
       w_ig=wgp(ig)*h/2;
       N_ig=N(ig,:);
       fe=fe+w_ig*(N_ig')*neumann;
   end
   
   %assembly
   fn(Te)=fn(Te)+fe;
   
end

end




