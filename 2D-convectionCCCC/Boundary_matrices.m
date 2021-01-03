function [Mout,Cout] = Boundary_matrices(X,T,Conv,referenceElement,velo)

T_boundary = OutflowBOundary(X,T,velo); 

elem = referenceElement.elem; 
nen = referenceElement.nen; 
p = referenceElement.degree; 

% 1D quadrature
ngaus = 3; 
zgp = [-sqrt(15)/5; 0; sqrt(15)/5]; 
wgp = [5/9   8/9   5/9];
% ngaus = 2; 
% zgp = [-1/sqrt(3); 1/sqrt(3)]; 
% wgp = [1   1];


% Number of nodes and elements
nPt = size(X,1); 
nElem = size(T_boundary,1); 

Mout = zeros(nPt);
Cout = zeros(nPt);

% Loop on elements
for i = 1:nElem
    ielem = T_boundary(i,1);
    Te = T(ielem,:); 
    Conve = Conv(Te,:);
    P1 = X(T_boundary(i,2),:); 
    P2 = X(T_boundary(i,3),:); 
    h = norm(P2 - P1); 
    n = T_boundary(i,4:5); 
    [N,Nxi,Neta] = ShapeFunc_aux(zgp,n,elem,p);
    
    Mout_e = zeros(nen); 
    Cout_e = zeros(nen);
    for ig = 1:ngaus
        N_ig    = N(ig,:);
        Nxi_ig  = Nxi(ig,:);
		Neta_ig = Neta(ig,:);
        %
        dvolu = wgp(ig)*h/2;
        % Derivatives of the shape functions on global coordinates
        Nx = Nxi_ig*2/h;
        Ny = Neta_ig*2/h;
        % velocity on the Gauss point
        a = N_ig*Conve; 
        % normal component of the velocity at the Gauss point
        an = n*a';
        % Contribution to element matrix
        Mout_e = Mout_e + an*(N_ig'*N_ig)*dvolu;
        Cout_e = Cout_e + an*(N_ig'*(a(1)*Nx+a(2)*Ny))*dvolu;
    end
    % Assembly
    Mout(Te,Te) = Mout(Te,Te) + Mout_e; 
    Cout(Te,Te) = Cout(Te,Te) + Cout_e; 
end



function T_boundary = OutflowBOundary(X,T,velo)

nElem = size(T,1); 
x1 = min(X(:,1)); x2 = max(X(:,1)); xM = (x1+x2)/2;
y1 = min(X(:,2)); y2 = max(X(:,2)); yM = (y1+y2)/2; 
T_boundary = zeros(nElem,5); ind = 1;  
if velo == 1
    for i = 1:nElem
        Te = T(i,:);
        Xe = X(Te,:); 
        xElem = Xe(:,1); 
        yElem = Xe(:,2); 
        xx1 = abs(xElem-x1); aux_x1 = find(xx1 < 1e-6);
        xx2 = abs(xElem-x2); aux_x2 = find(xx2 < 1e-6);
        yy1 = abs(yElem-y1); aux_y1 = find(yy1 < 1e-6); 
        yy2 = abs(yElem-y2); aux_y2 = find(yy2 < 1e-6);
        if length(aux_x1) == 2 && all(yElem(aux_x1) >= yM)
            T_boundary(ind,:) = [i, Te(aux_x1), -1, 0]; 
            ind = ind+1;
        end
        if length(aux_x2) == 2 && all(yElem(aux_x2) <= yM)
            T_boundary(ind,:) = [i, Te(aux_x2), 1, 0]; 
            ind = ind+1; 
        end
        if length(aux_y1) == 2 && all(xElem(aux_y1) <= xM)
            T_boundary(ind,:) = [i, Te(aux_y1), 0, -1]; 
            ind = ind+1; 
        end
        if length(aux_y2) == 2 && all(xElem(aux_y2) >= xM)
            T_boundary(ind,:) = [i, Te(aux_y2), 0, 1]; 
            ind = ind+1; 
        end
    end
    T_boundary = T_boundary(1:ind-1,:); 
elseif velo == 2
    for i = 1:nElem
        Te = T(i,:);
        Xe = X(Te,:); 
        xElem = Xe(:,1); 
        xx2 = abs(xElem-x2); aux_x2 = find(xx2 < 1e-6);
        if length(aux_x2) == 2 
            T_boundary(ind,:) = [i, Te(aux_x2), 1, 0]; 
            ind = ind+1; 
        end
    end
    T_boundary = T_boundary(1:ind-1,:); 
elseif velo == 3
    for i = 1:nElem
        Te = T(i,:);
        Xe = X(Te,:); 
        xElem = Xe(:,1); 
        yElem = Xe(:,2); 
        xx2 = abs(xElem-x2); aux_x2 = find(xx2 < 1e-6);
        yy2 = abs(yElem-y2); aux_y2 = find(yy2 < 1e-6);        
        if length(aux_x2) == 2 
            T_boundary(ind,:) = [i, Te(aux_x2), 1, 0]; 
            ind = ind+1; 
        end
        if length(aux_y2) == 2
            T_boundary(ind,:) = [i, Te(aux_y2), 0, 1]; 
            ind = ind+1; 
        end    
    end
    T_boundary = T_boundary(1:ind-1,:); 
end

% figure; PlotMesh(T,X,'b-');
% hold on 
% quiver(X(:,1),X(:,2),Conv(:,1),Conv(:,2))
% for i = 1:size(T_boundary,1)
%     ielem = T_boundary(i,1); 
%     nodes = T_boundary(i,2:3);
%     Te = T(ielem,:); 
%     Xe = X(Te,:); 
%     xm = mean(Xe(:,1)); ym = mean(Xe(:,2)); 
%     plot(xm,ym,'r*','LIneWidth',2)
%     plot(X(nodes,1), X(nodes,2), 'k*','LineWidth',2)
%     quiver(X(nodes(1),1), X(nodes(1),2), T_boundary(i,4), T_boundary(i,5),'k')
% end



function [N,Nxi,Neta] = ShapeFunc_aux(zgp_1D,n,elem,p)
ngaus = length(zgp_1D); 
if n(1) == 1
    zgp = [ ones(ngaus,1), zgp_1D]; 
elseif n(1) == -1
    zgp = [-ones(ngaus,1), zgp_1D]; 
elseif n(2) == 1
    zgp = [zgp_1D,  ones(ngaus,1)]; 
else
    zgp = [zgp_1D, -ones(ngaus,1)]; 
end
[N,Nxi,Neta] = ShapeFunc(elem,p,zgp);





