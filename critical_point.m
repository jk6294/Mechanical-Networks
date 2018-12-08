function [JgS, HES, GS, Jg] = critical_point(Xs, Xu, conn)
% Variables
ns = size(Xs,2);
nu = size(Xu,2);
X = [Xs Xu];
N = ns + nu;
LVal = sqrt((Xs(1,conn(:,1)) - Xu(1,conn(:,2))).^2 +...
            (Xs(2,conn(:,1)) - Xu(2,conn(:,2))).^2);

% Initialize Symbolic
XS = sym('x', [N, 1]); assume(XS, 'real');
YS = sym('y', [N, 1]); assume(YS, 'real');

% Constraints
g = [];
for i = 1:size(conn,1)
    g = [g; (XS(conn(i,1)) - XS(conn(i,2)+ns))^2 + (YS(conn(i,1)) - YS(conn(i,2)+ns))^2];
end

% Energy
E = (LVal(1) - sqrt((XS(conn(1,1)) - XS(conn(1,2)+ns))^2 + (YS(conn(1,1)) - YS(conn(1,2)+ns))^2))^2;
for i = 2:size(conn,1)
    E = E + (LVal(i) - sqrt((XS(conn(i,1)) - XS(conn(i,2)+ns))^2 + (YS(conn(i,1)) - YS(conn(i,2)+ns))^2))^2;
end

% Derivatives
Jg = jacobian(g);
HE = hessian(E);
JgS = double(subs(subs(Jg,XS,X(1,:)'),YS,X(2,:)'));
HES = double(subs(subs(HE,XS,X(1,:)'),YS,X(2,:)'));
G = [];
for i = 1:length(g)
    G(:,:,i) = hessian(g(i), [XS;YS]);
end
GS = double(subs(subs(G,XS,X(1,:)'),YS,X(2,:)'));

end
