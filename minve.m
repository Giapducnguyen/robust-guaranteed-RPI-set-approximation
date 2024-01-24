function [X, x_e] = minve(V)
% MINVE function: [X, x_e] = minve(V)
% Function to find the minimum volume ellipsoid E that contains a polytopic
% convex set S defined such that
% 
% S = conv(x_1,...,x_m).
% 
% S is assumed to be bounded with non-empty interior.
% The input of the function is V, a matrix containing the coordinates of
% the points x_i on the row i, for all i in {1,...,m}.
% The outputs of the function are X and x_e respectivelly the ellipsoid
% matrix and center, such that the ellipsoid E is defined as follows
% 
% E = {x|(x-x_e)'*X^{-1}*(x-x_e)<=1}.

% R. Guicherd - December 2019
%% Sanity check

% Check for matrix type
if ~ismatrix(V)
    error('V is not a matrix!') 
end

% Check for real type
if ~isreal(V)
    error('V is not real!') 
end

% Problem dimensions
[m, n] = size(V); 

%% Optimization problem
% Optimization variables
X = sdpvar(n,n);
x_e = sdpvar(1,n,'full');

% Constraints
Cons = [];
for i = 1:1:m
    Cons = [Cons, [1 V(i,:)-x_e; (V(i,:)-x_e)' X]>= 0]; %#ok<AGROW>
end
clearvars i

% Options
Opts = sdpsettings('solver', 'sedumi', 'sedumi.eps', 1e-9, 'verbose', 0);

% Optimization
optimize(Cons, trace(X), Opts)

% Extract the values of X and x_e to output
X = value(X);
x_e = transpose(value(x_e));

end
%%%%% END OF MINVE FUNCTION %%%%%