% This program solves a convection problem 
% in a square domain [-0.5,0.5]x[-0.5,0.5] 

clear; close all; clc


disp(' ')
disp('This program solves a pure convection equation in [-0.5,0.5]x[-0.5,0.5] ')


dom = [-0.5,0.5,-0.5,0.5]; 

% Element type
disp(' ')
elem = 0; p = 1; 
referenceElement = SetReferenceElement(elem,p); 

h = cinput('Spatial mesh size',0.05);

nx=(dom(2)-dom(1))/h;
ny=(dom(4)-dom(3))/h;

[X,T] = CreateMesh(dom,nx,ny,referenceElement);
figure; PlotMesh(T,X,'b-');

disp(' ')
disp('The following velocity fields can be considered');
disp('   [1] v(x,y) = (-x,-y)');
disp('   [2] v(x,y) = (1,0)');
disp('   [3] v(x,y) = (1,1)/sqrt(2)'); 
velo = cinput('Choose the convection field to be used on the computation',1); 
Conv = ComputeVelocity(X,velo);




% Method used for solving the problem
disp(' ')
disp ('The problem can be solved using one of the following methods: ');

disp ('    [1] Crank-Nicolson + Galerkin')
disp ('    [2] Lax-Wendroff + Galerkin')
disp ('    [3] Lax-Wendroff with lumped mass matrix + Galerkin')
method = cinput('Method = ',1);

figure(1); hold on; 
quiver(X(:,1),X(:,2),Conv(:,1),Conv(:,2))
plot(dom([1,2,2,1,1]), dom([3,3,4,4,3]), 'k')
axis equal

% Time discretization
tEnd = cinput('End time', 2*pi);
nStep = cinput('Numnber of time-steps', 120);
dt = tEnd / nStep; 
Courant = sqrt(Conv(1)^2+Conv(2)^2)*dt/h; 
disp(['Courant number: ', num2str(Courant)]); 

 
% Initial condition
u0 = InitialCondition(X); 

xx  = reshape(X(:,1), nx+1, ny+1)'; 
yy  = reshape(X(:,2), nx+1, ny+1)'; 
u0_dib = reshape(u0, nx+1, ny+1)';  
figure(2); clf; 
surface(xx,yy,u0_dib,'FaceColor','interp');
set(gca,'FontSize',16)
grid on
view(3)

figure(3); clf;
[C,h]=contour(xx,yy,u0_dib,[-0.1,-0.01,0.1:0.1:1.0]);
clabel(C,h,'FontSize',12);
axis equal; axis(dom);
set(gca,'FontSize',16,'XTick',-0.5:0.25:0.5,'YTick',-0.5:0.25:0.5)

[M,K,C] = FEM_matrices(X,T,Conv,referenceElement); 
[Mout,Cout] = Boundary_matrices(X,T,Conv,referenceElement,velo);
f = zeros(size(X,1),1);

% Boundary conditions (Lagrange multipliers)
[ADir,bDir,nDir] = BoundaryConditions(X,velo);


[A,B,methodname]=Method(M,C,K,dt,method);

Ktot = [A ADir';ADir zeros(nDir,nDir)];


u=zeros(size(X,1),nStep+1);
u(:,1) = u0;
for n = 1:nStep
    ftot=[(B*u(:,n) + f);bDir];
    sol = Ktot\ftot;
    Du = sol(1:length(f));
    u(:,n+1) = u(:,n) + Du; 
end



%xx  = reshape(X(:,1), ny+1, nx+1)'; 
%yy  = reshape(X(:,2), ny+1, nx+1)'; 
figure(4);  
for n=1:size(u,2)
    clf;
uu = reshape(u(:,n), nx+1, ny+1)';  
surface(xx,yy,uu,'FaceColor','interp');
set(gca,'FontSize',16)
grid on
view(3)
pause(0.05)
end
    




% Finish the code to be able to solve the problem
% ...
% ...
% ...

