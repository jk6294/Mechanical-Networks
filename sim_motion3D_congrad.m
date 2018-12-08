function [XMot, fC] = sim_motion3D_congrad(Xs, Xu, conn, LVal, delS, n, XN, XF, pV, kL)

if ~exist('pV', 'var')
    pV = 0;
end
if ~exist('kL', 'var')
    kL = ones(length(LVal), 1);
end


%% Parameters
ns = size(Xs,2);
nu = size(Xu,2);
d = size(Xs,1);
X = [Xs Xu]';
N = ns + nu;
Us = (XF - Xs(:,XN))';


%% Symbolic
if(d==2)
    sInds = [XN, XN+N];
    nInds = setdiff(1:(2*N), sInds);
    % Create Symbolic Variables
    XS = sym('x', [N, 1]); assume(XS, 'real');
    YS = sym('y', [N, 1]); assume(YS, 'real');
    XNN = setdiff(1:N, XN);
    XV = [XS(XNN); YS(XNN)];
    % Energy
    E = kL(1)*(LVal(1) - sqrt((XS(conn(1,1)) - XS(conn(1,2)+ns))^2 + (YS(conn(1,1)) - YS(conn(1,2)+ns))^2))^2;
    for i = 2:size(conn,1)
        E = E + kL(i)*(LVal(i) - sqrt((XS(conn(i,1)) - XS(conn(i,2)+ns))^2 +...
                                      (YS(conn(i,1)) - YS(conn(i,2)+ns))^2))^2;
    end
elseif(d==3)
    sInds = [XN, XN+N, XN+2*N];
    nInds = setdiff(1:(3*N), sInds);
    % Create Symbolic Variables
    XS = sym('x', [N, 1]); assume(XS, 'real');
    YS = sym('y', [N, 1]); assume(YS, 'real');
    ZS = sym('z', [N, 1]); assume(ZS, 'real');
    XNN = setdiff(1:N, XN);
    XV = [XS(XNN); YS(XNN); ZS(XNN)];
    % Energy
    E = kL(1)*(LVal(1) - sqrt((XS(conn(1,1)) - XS(conn(1,2)+ns))^2 + (YS(conn(1,1)) - YS(conn(1,2)+ns))^2 + (ZS(conn(1,1)) - ZS(conn(1,2)+ns))^2))^2;
    for i = 2:size(conn,1)
        E = E + kL(i)*(LVal(i) - sqrt((XS(conn(i,1)) - XS(conn(i,2)+ns))^2 +...
                                      (YS(conn(i,1)) - YS(conn(i,2)+ns))^2 +...
                                      (ZS(conn(i,1)) - ZS(conn(i,2)+ns))^2))^2;
    end
end

EJ = jacobian(E, XV)';
EJ2 = hessian(E, XV);
Ef = matlabFunction(E, 'Optimize', false);
EJf = matlabFunction(EJ, 'Optimize', false);
EJ2f = matlabFunction(EJ2, 'Optimize', false);


%% Descent
% Initial Conditions
XP = zeros([d*N,n]);
XP(:,1) = X(:);
XPS = num2cell(XP(:,1));
fC = zeros([1, n]); fC(1) = Ef(XPS{:});

fprintf([repmat('.',1,n-1) '\n\n']);
for i = 2:n
    fprintf('\b=\n');
    XP(:,i) = XP(:,i-1);
    XP(sInds,i) = XP(sInds,i) + delS*Us(:);
    
    EJfe = ones(d*N,1);
    while(norm(EJfe) > 1e-12)
        % Convert to cell
        XPS = num2cell(XP(:,i));
        % Evaluate Gradient
        EJfe = EJf(XPS{:});
        % Newton Method
        EJ2fe = EJ2f(XPS{:});
        XP(nInds,i) = XP(nInds,i) - pinv(EJ2fe)*EJfe;
    end
    fC(i) = Ef(XPS{:});
end

XMot = zeros(d, N, n);
for i = 1:d
    XMot(i,:,:) = XP([1:N]+(i-1)*N,:);
end


%% Plot
% Plot Parameters
ms = 4;         % Marker Size
% Energy
C = parula(100);
EInd = discretize(fC,100);

if(pV ~= 0)
    hold on;
    if(d==2)
        % Traversal
        for i = 1:n
            plot(XMot(1,1:ns,i), XMot(2,1:ns,i), '.', 'linewidth', ms, 'markersize', ms, 'color', C(EInd(i),:))
            plot(XMot(1,(1:nu)+ns,i), XMot(2,(1:nu)+ns,i), '.', 'linewidth', ms, 'markersize', ms, 'color', C(EInd(i),:));
        end
    elseif(d==3)
        for i = 1:n
            plot3(XMot(1,1:ns,i), XMot(2,1:ns,i), XMot(3,1:ns,i),...
                '.', 'linewidth', ms/10, 'markersize', ms/10, 'color', C(EInd(i),:));
            plot3(XMot(1,(1:nu)+ns,i), XMot(2,(1:nu)+ns,i), XMot(3,(1:nu)+ns,i),...
                '.', 'linewidth', ms/10, 'markersize', ms/10, 'color', C(EInd(i),:));
        end
    end
    
    hold off;
end


end