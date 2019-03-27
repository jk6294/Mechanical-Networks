function [S, Xd] = rigidity(Xs, Xu, conn)
% Function to generate rigidity matrix as a function of node positions
% Inputs
%       Xs:     d x n matrix of coordinates for specified nodes
%       Xu:     d x m matrix of coordinates for unspecified nodes
%       conn:   E x 2 matrix of k edges connecting node i to node j
% Outputs
%       S:      k x k x s matrix of quadratic forms for s vectors in left
%               nullspace of rigidity matrix
%       Xd:     d(n+m) x k matrix of conformational motions in nullspace of
%               rigidity matrix
%
% Collect relevant parameters
n = size(Xs,2);         % Number specified nodes
m = size(Xu,2);         % Number unspecified nodes
N = n+m;                % Total nodes
X = [Xs Xu]';
d = size(Xs,1);         % Dimension of space
g = [];                 % Vector of constraints

% Compute Rigidity: 2 Dimensions
if(d==2)
XS = sym('x', [N, 1]); assume(XS, 'real');
YS = sym('y', [N, 1]); assume(YS, 'real');
I = zeros(size(conn,1), N);
for i = 1:size(conn,1)
    I(i,conn(i,1)) =  1;
    I(i,conn(i,2)) = -1;
end
J = [(XS(conn(:,1)) - XS(conn(:,2))) .* I,...
     (YS(conn(:,1)) - YS(conn(:,2))) .* I];
% Rigid Body Motions
r1 = [ones([N,1]); zeros([N,1])];           % x-translation
r2 = [zeros([N,1]); ones([N,1])];           % y-translation
r3 = [sqrt(X(:,1).^2 + X(:,2).^2); sqrt(X(:,1).^2 + X(:,2).^2)] .*...
     [-sin(atan2(X(:,2), X(:,1))); cos(atan2(X(:,2), X(:,1)))];
G = [r1 r2 r3];
R = matlabFunction(J, 'Optimize', false, 'Sparse', true, 'Vars', {[XS; YS]});

% Compute Rigidity: 3 Dimensions
elseif(d==3)
XS = sym('x', [N, 1]); assume(XS, 'real');
YS = sym('y', [N, 1]); assume(YS, 'real');
ZS = sym('z', [N, 1]); assume(ZS, 'real');
I = zeros(size(conn,1), N);
for i = 1:size(conn,1)
    I(i,conn(i,1)) =  1;
    I(i,conn(i,2)) = -1;
end
J = [(XS(conn(:,1)) - XS(conn(:,2))) .* I,...
     (YS(conn(:,1)) - YS(conn(:,2))) .* I,...
     (ZS(conn(:,1)) - ZS(conn(:,2))) .* I];
% Rigid Body Motions
r1 = [ones([N,1]); zeros([N,1]); zeros([N,1])];           % x-translation
r2 = [zeros([N,1]); ones([N,1]); zeros([N,1])];           % y-translation
r3 = [zeros([N,1]); zeros([N,1]); ones([N,1])];           % z-translation
r4 = [sqrt(X(:,1).^2 + X(:,2).^2); sqrt(X(:,1).^2 + X(:,2).^2); zeros([N,1])] .*...
     [-sin(atan2(X(:,2), X(:,1))); cos(atan2(X(:,2), X(:,1))); zeros([N,1])];
r5 = [zeros([N,1]); sqrt(X(:,2).^2 + X(:,3).^2); sqrt(X(:,2).^2 + X(:,3).^2)] .*...
     [zeros([N,1]); -sin(atan2(X(:,3), X(:,2))); cos(atan2(X(:,3), X(:,2)))];
r6 = [sqrt(X(:,3).^2 + X(:,1).^2); zeros([N,1]); sqrt(X(:,3).^2 + X(:,1).^2)] .*...
     [-sin(atan2(X(:,1), X(:,3))); zeros([N,1]); cos(atan2(X(:,1), X(:,3)))];
G = [r1 r2 r3 r4 r5 r6];
R = matlabFunction(J, 'Optimize', false, 'Sparse', true, 'Vars', {[XS;YS;ZS]});
end

% Compute non-rigid body motions in null(R)
RS = full(R(X(:)));
NS = null(RS);
W = null(RS');
Xd = double(NS*null(G'*NS));
k = size(Xd,2);
s = size(W,2);
S = zeros(k,k,s);
for i = 1:s
    w = W(:,i);
    for j = 1:k
        S(j,:,i) = w'*R(Xd(:,j))*Xd;
    end
end

end