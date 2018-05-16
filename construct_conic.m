function [Q, W, v0, err] = construct_conic(X, U)
% Code for generating points along the solution space given specified nodes
% Able to construct conic in any dimension for any number of nodes
% 
% Inputs
% X: dxN (dimensions x nodes) array of specified node positions
% U: dxN (dimensions x nodes) array of specified node motions
%
% Outputs
% Q: m+1 x m+1      conic matrix for the m-dimensional solution space
% W: 2d+1 x m       nullspace matrix of homogenous solutions
% v0: 2d+1 x 1      vector for the non-homogeneous solution
% err: scalar:      distance of particular solution to solution space

% Initialize Values and Matrices
d = size(X,1);                          % Dimensions
N = size(X,2);                          % Nodes
At = [U' X' -ones(N,1)];                % \tilde{A} from main text
O = [zeros(d) eye(d) zeros(d,1);...     % O from supplement 
     eye(d)   zeros(d,d+1);...
     zeros(1,2*d+1)]/2;
b = sum(X.*U)';
p = [zeros([1, 2*d]), 1];               % Inidicator vector


% Find homogeneous and non-homogenous solutions to linear problem
W = null(At);           % Nullspace gives homogeneous solutions
v0 = pinv(At)*b;        % Particular Solution
err = norm(At*v0 - b);  % error

% Construct Conic Matrix by Applying Nonlinear Constraint
AP = W'*O*W;
BP = (2*v0'*O - p)*W/2;
CP = (v0'*O - p)*v0;
Q = [AP, BP';...
     BP, CP];
end