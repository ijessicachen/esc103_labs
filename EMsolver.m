
function [t,x,y] = EMsolver(A,x_0,y_0,T,N) 
%EMSolver Solves an initial value problem using Euler's Method
%   This function takes as input an A matrix that describes a specific 
%   differential equation problem of form Z'=AZ with Z a vector that
%   contains the variables (e.g., x(t), y(t)) in the differential equation.
%   It also takes initial values x_0=x(0) and y_0=y(0), the end time of the
%   function approximation T (with t in [0,T]), and the number of timesteps
%   N that should be used for the approximation.

    % compute a vector of timesteps of size dt=T/N
    dt = T/N; 
    t = 0:dt:T; 
     
    % initialize the vector that will contain the solution x(t) and y(t):
    % use a vector of NaN (Not a Number) entries, then set x(0) and y(0)
    SOL = NaN(2,length(t)); 
    SOL(1,1) = x_0; 
    SOL(2,1) = y_0; 
     
    % iteratively approximate the solution values at each timestep by
    % computing the slope using Z'=AZ and multiplying this by the timestep
    % size, then adding the result to the approximate value computed at the
    % previous timestep
    for n = 2:length(t)
        SOL(:,n) = SOL(:,n-1) + dt*A*SOL(:,n-1); 
    end 
     
    % return the x(t) and y(t) parts of the solution as different vectors
    x = SOL(1,:); 
    y = SOL(2,:); 
end