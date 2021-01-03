function u0 = InitialCondition(X)

nPt = size(X,1);
u0 = zeros(nPt,1); 

sigma = 0.2; xref = 1/6; 
xdim = ( X(:,1) - xref) / sigma; 
ydim = ( X(:,2) - xref) / sigma; 
ind = find( (xdim.^2 + ydim.^2)<=1 ); 
u0(ind) = 0.25*(1 + cos(pi*xdim(ind))).*(1 + cos(pi*ydim(ind))); 
