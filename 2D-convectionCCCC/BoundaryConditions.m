function [A,b,nDir] = BoundaryConditions(X,velo)


x1 = min(X(:,1)); x2 = max(X(:,1)); xM = (x1+x2)/2;
y1 = min(X(:,2)); y2 = max(X(:,2)); yM = (y1+y2)/2; 

if velo == 1
    nodes_x1 = find( (abs(X(:,1) - x1)<1e-6) & (X(:,2)<=yM) ); 
    nodes_x2 = find( (abs(X(:,1) - x2)<1e-6) & (X(:,2)>=yM) ); 
    nodes_y1 = find( (abs(X(:,2) - y1)<1e-6) & (X(:,1)>=xM) ); 
    nodes_y2 = find( (abs(X(:,2) - y2)<1e-6) & (X(:,1)<=xM) ); 
    nodesDir = unique([nodes_x1; nodes_x2; nodes_y1; nodes_y2]); 
elseif velo == 2
    nodesDir = find( (abs(X(:,1) - x2)<1e-6)); 
elseif velo == 3
    nodes_x2 = find( (abs(X(:,1) - x2)<1e-6)); 
    nodes_y2 = find( (abs(X(:,2) - y2)<1e-6)); 
    nodesDir = unique([nodes_x2;nodes_y2]);
end

C = [nodesDir, zeros(size(nodesDir))]; 

nDir = size(C,1); 
nPt = size(X,1); 
A = zeros(nDir,nPt);
A(:,C(:,1)) = eye(nDir); 
b = C(:,2); 
