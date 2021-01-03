function [ A,B,methodName ] = Method(M,C,K,dt,method)

switch method
    case 1 % Crank-Nicolson + Galerkin 
        
    A = M + 0.5*dt*C;
    B = -dt*C; 
    methodName = 'CN';
    
    case 2  % Lax-Wendroff + Galerkin
     A = M;
     B = -dt*C- 0.5*dt^2*K; 
     methodName = 'LW';
    case 3 
     A = diag(sum(M));
     B = -dt*C- 0.5*dt^2*K; 
     methodName = 'LW';
        
        
end




end

