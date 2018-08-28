function [XMot, fC] = sim_motion3D_congrad(Xs, Xu, conn, LVal, delS, n, XN, XF, kL)

if ~exist('kL', 'var')
    kL = ones(length(LVal), 1);
end


%% Parameters
ns = size(Xs,2);
nu = size(Xu,2);
X = [Xs Xu];
N = ns + nu;
Us = (XF - Xs(:,XN))';


%% Symbolic
sInds = [XN, XN+N, XN+2*N];
nInds = setdiff(1:(3*N), sInds);
XNN = setdiff(1:N, XN);

% Create Symbolic Variables
XS = sym('x', [N, 1]); assume(XS, 'real');
YS = sym('y', [N, 1]); assume(YS, 'real');
ZS = sym('z', [N, 1]); assume(ZS, 'real');

XV = [XS(XNN); YS(XNN); ZS(XNN)];

% Energy
E = kL(1)*(LVal(1) - sqrt((XS(conn(1,1)) - XS(conn(1,2)+ns))^2 + (YS(conn(1,1)) - YS(conn(1,2)+ns))^2 + (ZS(conn(1,1)) - ZS(conn(1,2)+ns))^2))^2;
for i = 2:size(conn,1)
    E = E + kL(i)*(LVal(i) - sqrt((XS(conn(i,1)) - XS(conn(i,2)+ns))^2 +...
                                  (YS(conn(i,1)) - YS(conn(i,2)+ns))^2 +...
                                  (ZS(conn(i,1)) - ZS(conn(i,2)+ns))^2))^2;
end
EJ = jacobian(E, XV)';
EJ2 = hessian(E, XV);

% Ef = matlabFunction(E, 'Optimize', false, 'Vars', {XV});
% EJf = matlabFunction(EJ, 'Optimize', false, 'Vars', {XV});
Ef = matlabFunction(E, 'Optimize', false);
EJf = matlabFunction(EJ, 'Optimize', false);
EJ2f = matlabFunction(EJ2, 'Optimize', false);


%% Descent
% Initial Conditions
XP = zeros([3*N,n]);
XP(:,1) = [X(1,:) X(2,:) X(3,:)]';
XPS = num2cell(XP(:,1));
fC = zeros([1, n]); fC(1) = Ef(XPS{:});

fprintf([repmat('.',1,n-1) '\n\n']);
for i = 2:n
    fprintf('\b=\n');
    XP(:,i) = XP(:,i-1);
    XP(sInds,i) = XP(sInds,i) + delS*Us(:);
    
    EJfe = ones(3*N,1);
    while(norm(EJfe) > 1e-12)
        % Convert to cell
        XPS = num2cell(XP(:,i));
        % Evaluate Gradient
        EJfe = EJf(XPS{:});
%         XP(nInds,i) = XP(nInds,i) - EJfe(nInds)*80*delS;
%         disp(norm(EJfe(nInds)));

        % Newton Method
        EJ2fe = EJ2f(XPS{:});
        XP(nInds,i) = XP(nInds,i) - pinv(EJ2fe)*EJfe;
    end
    fC(i) = Ef(XPS{:});
end
xC = XP(1:N,:);
yC = XP([1:N]+N,:);
zC = XP([1:N]+2*N,:);
XC = zeros(3, N, n);
XC(1,:,:) = reshape(xC, [1, N, n]);
XC(2,:,:) = reshape(yC, [1, N, n]);
XC(3,:,:) = reshape(zC, [1, N, n]);

XMot = XC;

end