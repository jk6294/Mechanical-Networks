function [Q, W, v0, err] = construct_conic_finite(X0, XF)
% Code for generating points along the solution space given specified nodes
% Able to construct conic in any dimension for any number of nodes
% 
% Inputs
% X0: dxN (dimensions x nodes) array of specified node initial positions
% XF: dxN (dimensions x nodes) array of specified node final positions
%
% Outputs
% Q: m+1 x m+1      conic matrix for the m-dimensional solution space
% W: 2d+1 x m       nullspace matrix of homogenous solutions
% v0: 2d+1 x 1      vector for the non-homogeneous solution
% err: scalar:      distance of particular solution to solution space

% Initialize Values and Matrices
d = size(X0,1);                          % Dimensions
N = size(X0,2);                          % Nodes
M = sym([2*X0', -2*XF', -ones(N,1)]);    % M from main text
O = [eye(d)   zeros(d) zeros(d,1);...    % O from supplement 
     zeros(d)  -eye(d) zeros(d,1);...
     zeros(1,2*d+1)];
b = sum(X0.^2 - XF.^2)';
p = [zeros([1, 2*d]), 1];                % Inidicator vector

% Find homogeneous and non-homogenous solutions to linear problem
W = null(M);           % Nullspace gives homogeneous solutions
v0 = pinv(M)*b;        % Particular Solution
err = norm(M*v0 - b);  % error

% Construct Conic Matrix by Applying Nonlinear Constraint
AP = W'*O*W;
BP = (2*v0'*O - p)*W/2;
CP = (v0'*O - p)*v0;
Q = [AP, BP';...
     BP, CP];
 
% Reconvert to Double
Q = double(Q);
W = double(W);
v0 = double(v0);
end