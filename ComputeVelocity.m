function Conv = ComputeVelocity(X,velo)
% 
% Velocity at the mesh nodes

if velo == 1
    Conv = [-X(:,2), X(:,1)];
elseif velo == 2
    n = size(X,1); 
    Conv = [ones(n,1), zeros(n,1)];
elseif velo == 3
    Conv = ones(size(X)) / sqrt(2); 
else
    error('not available velocity')
end