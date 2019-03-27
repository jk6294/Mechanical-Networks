function [Q, W, v0, err] = construct_conic(Xs, dXs, a)
% Code for generating points along the solution space given specified 
% nodes, and infinitesimal or finite motions. If z is not 0 or 1,
% default to infinitesimal
% 
% Inputs
% Xs:       d x n matrix of specified node initial positions
% dXs:      d x n matrix of specified node motions or final positions
% a:        1 x 1 scalar as 0 (infinitesimal) or 1 (finite)
%
% Outputs
% Q:    m+1 x m+1       conic matrix for the m-dimensional solution space
% W:    2d+1 x m        nullspace matrix of homogenous solutions
% v0:   2d+1 x 1        vector for the non-homogeneous solution
% err:  1 x 1           distance of particular solution to solution space

% Initialize Values and Matrices
d = size(Xs,1);                         % Dimensions
n = size(Xs,2);                         % Nodes
p = [zeros([2*d, 1]); 1];               % Inidicator vector

% Nonlinear constraint and constant based on infinitesimal or finite
if(a==1)
    M = sym([2*Xs', -2*dXs', -ones(n,1)]);  % M from main text
    O = [eye(d)   zeros(d) zeros(d,1);...   % Finite nonlinear constraint
         zeros(d)  -eye(d) zeros(d,1);...
         zeros(1,2*d+1)];
    b = sum(Xs.^2 - dXs.^2)';
else
    M = sym([dXs' Xs' -ones(n,1)]);     % Linearized matrix M from main
    O = [zeros(d) eye(d) zeros(d,1);... % Infinitesimal constraint
         eye(d)   zeros(d,d+1);...
         zeros(1,2*d+1)]/2;
    b = sum(sym(Xs.*dXs))';
end

% Find homogeneous and non-homogenous solutions to linear problem
W = null(M);           % Nullspace gives homogeneous solutions
v0 = pinv(M)*b;        % Particular solution v* in main text
err = norm(M*v0 - b);  % Error in existence of partincular solution
if(err > 1e-12)
    disp('Error in reconstructing particular solution');
end

% Construct Conic Matrix by Applying Nonlinear Constraint
AP = W'*O*W;
BP = W'*(2*O*v0 - p)/2;
CP = (v0'*O - p')*v0;
Q = [AP,  BP;...
     BP', CP];

% Reconvert to Double
Q = double(Q);
W = double(W);
v0 = double(v0);
end