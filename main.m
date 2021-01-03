% This program solves a convection-diffusion problem 
% in a square domain [0,2]x[0,3] 

clear; close all; clc

disp ('The following problem can be solved:');
disp ('     [0] Steady-convection-diffusion-reaction');


disp ('     [1] Unsteady-convection-diffusion-reaction');

addpath('FiniteElement')


Transient = cinput('Choose the problem type to be solved', 1);


disp(' ')
disp ('The following methods can be used:');
disp ('     [1] Galerkin');
disp ('     [2] SUPG');
disp ('     [3] OSS');
method = cinput ('Choose a method for solving the problem: ',3);

if method==2
    para.method = 'SUPG';
elseif method == 3
    para.method='OSS';
else
    para.method='Galerkin';
end
    



if Transient ==0
    disp(' ')
    disp('This program solves a convection-diffusion equation on [0,2]x[0,3]')
    disp(' ')
    disp('No source term is considered');
    
    % PDE coefficients
    a=cinput('Norm of Convective term',1);
    disp('Convection velocity is');
    velo = [a*cos(pi/6),a*sin(pi/6)];
    disp(velo);
    nu = cinput('Diffusion coefficient', 0.01);
    
    sigma=cinput('Reaction coefficient',1);
    neumann=cinput('Neumann boundary value',0);
    disp('')
    
    dom = [0,2,0,3];
    
    para.dom = dom;
    para.velo = velo;
    para.a = a;
    para.nu = nu;
    para.sigma=sigma;
    
    % Element type
    elem = 0; para.p = 1;
    referenceElement = SetReferenceElement(elem,para.p);
    disp('')
    para.h = cinput('Spatial mesh size',0.2);
    
    para.nx=(para.dom(2)-para.dom(1))/para.h;
    parany=(para.dom(4)-para.dom(3))/para.h;
    
    [X,T] = CreateMesh(dom,para.nx,para.ny,referenceElement);
    PlotMesh(T,X,'b-');
    
    
    Pe = a*para.h/(2*nu);
    disp(' ')
    disp(strcat('Peclet number: ',num2str(Pe)))
    
    
    
    % SYSTEM RESULTING OF DISCRETIZING THE WEAK FORM
    [K,f] = FEM_system(X,T,referenceElement,para);
    
    
    % BOUNDARY CONDITIONS
    % Boundary conditions are imposed using Lagrange multipliers
    nodes_y0 = [1:para.nx+1]';                           % Nodes on the boundary y=0
    nodes_x1 = [(para.nx+1):para.nx+1:(para.ny+1)*(nx+1)]' ;     % Nodes on the boundary x=1
    nodes_y1 = [para.ny*(para.nx+1)+para.nx:-1:para.ny*(para.nx+1)+1]' ;     % Nodes on the boundary y=1
    nodes_x0 = [(para.ny)*(para.nx+1)+1:-(para.nx+1):1]';     % Nodes on the boundary x=0
    
    
    % nodes on which solution is u=1
    nodesDir1 = nodes_x0( X(nodes_x0,2)>=1.5 );
    %nodesDir1 = nodes_x0  ;%( X(nodes_x0,2)>=1.5 );
    % nodes on which solution is u=0
    %nodesDir0 = [ nodes_x1 ];
    % Boundary condition matrix
    C = [nodesDir1, 2*ones(length(nodesDir1),1)];
    nDir = size(C,1);
    neq  = size(f,1);
    A = zeros(nDir,neq);
    A(:,C(:,1)) = eye(nDir);
    b = C(:,2);
    
    %apply neumann boundary condtion
    fn=neumannBC(X,nodes_x1,neumann);
    
    
    
    f=f+fn;
    
    % SOLUTION OF THE LINEAR SYSTEM
    % Entire matrix
    Ktot = [K A';A zeros(nDir,nDir)];
    ftot = [f;b];
    
    sol = Ktot\ftot;
    Temp = sol(1:neq);
    multip = sol(neq+1:end);
    
    
    % POSTPROCESS
    figure(2), clf
    xx  = reshape(X(:,1), para.nx+1, para.ny+1)';
    yy  = reshape(X(:,2), para.nx+1, para.ny+1)';
    sol = reshape(Temp, para.nx+1, para.ny+1)';
    surface(xx,yy,sol,'FaceColor','interp');
    set(gca, 'xTick',0:0.5:3, 'yTick',0:0.5:3, 'FontSize',12)
    xlabel('x','FontSize',14);
    ylabel('y','FontSize',14);
    zlabel('u','FontSize',14);
    
    grid on; view(3)
    
else
    
    
    
  
    
    disp('This program solves a convection-diffusion-reaction equation')
    para.dom = [-0.5,0.5,-0.5,0.5]; 
    
    % shock capturing
    
    para.alpha=0.8;
    
    % Element type
    disp(' ')
    elem = 0;
    
    para.p = cinput('Element Degress',1);
    referenceElement = SetReferenceElement(elem,para.p);
     
    para.h = cinput('Spatial mesh size',0.05);
    
    para.nx=(para.dom(2)-para.dom(1))/para.h;
    para.ny=(para.dom(4)-para.dom(3))/para.h;
    
    [X,T] = CreateMesh(para.dom,para.nx,para.ny,referenceElement);
    figure; PlotMesh(T,X,'b-');
    
    disp(' ')
    disp('The following velocity fields can be considered');
    disp('   [1] v(x,y) = (-y,x)');
    disp('   [2] v(x,y) = (1,0)');
    disp('   [3] v(x,y) = (1,1)/sqrt(2)');
    para.velo = cinput('Choose the convection field to be used on the computation',2);
   
    
    para.nu = cinput('Diffusion coefficient nu', 0.01);
    
    para.sigma = cinput('Reaction', 0);
    
    
      
    % Temporal discretization
    tEnd = cinput('End time', 2*pi);
    nStep = cinput('Numnber of time-steps', 110);
    para.dt = tEnd / nStep;
    para.theta=0.5;
    
    tic;%tic1
    t1=clock;
    % Initial condition
    
    para.u0 = InitialCondition(X);
    xx  = reshape(X(:,1), para.p*para.nx+1, para.p*para.ny+1)';
    yy  = reshape(X(:,2), para.p*para.nx+1, para.p*para.ny+1)';
    u0_dib = reshape(para.u0, para.p*para.nx+1, para.p*para.ny+1)';
    figure(2); clf;
    surface(xx,yy,u0_dib,'FaceColor','interp');
    set(gca,'FontSize',16)
    grid on
    view(3)
    
    % Internal iteration Picard method
    %para.Picard=cinput('Picard iteration', 1);
    para.Picard=1;
    para.tol=1e-6;
    
    

    u=Method(para,X,T,referenceElement,nStep);
    
    disp(['total time for the code:',num2str(etime(clock,t1))]);
   
%     figure(3); clf;
%     [C,h]=contour(xx,yy,u0_dib,[-0.1,-0.01,0.1:0.1:1.0]);
%     clabel(C,h,'FontSize',12);
%     axis equal; axis(para.dom);
%     set(gca,'FontSize',16,'XTick',0:0.25:2,'YTick',0:0.25:3)
    
    figure(4); 
    for n=1:size(u,2)
        clf;
        uu = reshape(u(:,n), para.p*para.nx+1, para.p*para.ny+1)';
        surface(xx,yy,uu,'FaceColor','interp');
        set(gca,'FontSize',16)
        grid on
        view(3)
        pause(0.1)
    end
end




