function [XC, fC] = sim_motion3D_con(Xs, Xu, conn, delS, n, XN, XF)
% Function for simulating the motion of frames
% Inputs
%       Xs      3 X ns matrix of coordinates for specified nodes
%       Xu      3 X nu matrix of coordinates for unspecified nodes
%       conn:   k X 2 matrix of k edges connecting two nodes
%       delS:   Instantaneous length of motion
%       n:      Number of steps to move forward
%       XN      1 x 2 vector of specified nodes to move
%       XF      3 x 2 vector of final states for specified nodes
% Outputs
%       XC:     3 X N X n matrix of node motions across n time steps
%       fC:     1 X n matrix of motion errors


%% Lengths
ns = size(Xs,2);
nu = size(Xu,2);
X = [Xs Xu];
N = ns + nu;
LVal = sqrt((Xs(1,conn(:,1)) - Xu(1,conn(:,2))).^2 +...
            (Xs(2,conn(:,1)) - Xu(2,conn(:,2))).^2 +...
            (Xs(3,conn(:,1)) - Xu(3,conn(:,2))).^2);
Us = XF - Xs(:,XN);
        
        
%% Symbolic 
sInds = [XN, XN+N, XN+2*N];
nInds = setdiff(1:(3*N), sInds);
XNN = setdiff(1:N, XN);

XS = sym('x', [N, 1]); assume(XS, 'real');
YS = sym('y', [N, 1]); assume(YS, 'real');
ZS = sym('z', [N, 1]); assume(ZS, 'real');

XSC = XS(XN);
YSC = YS(XN);
ZSC = ZS(XN);

XSN = XS(XNN);
YSN = YS(XNN);
ZSN = ZS(XNN);



% Constraints
g = [];
for i = 1:size(conn,1)
    g = [g; (XS(conn(i,1)) - XS(conn(i,2)+ns))^2 + (YS(conn(i,1)) - YS(conn(i,2)+ns))^2 + (ZS(conn(i,1)) - ZS(conn(i,2)+ns))^2];
end
disp('1');
% J = jacobian(g, [XS; YS; ZS])/2;
disp('2');
% Jf = matlabFunction(J, 'Optimize', false, 'Sparse', true);
disp('3');

% Energy
E = (LVal(1) - sqrt((XS(conn(1,1)) - XS(conn(1,2)+ns))^2 + (YS(conn(1,1)) - YS(conn(1,2)+ns))^2 + (ZS(conn(1,1)) - ZS(conn(1,2)+ns))^2))^2;
for i = 2:size(conn,1)
    E = E + (LVal(i) - sqrt((XS(conn(i,1)) - XS(conn(i,2)+ns))^2 + (YS(conn(i,1)) - YS(conn(i,2)+ns))^2 + (ZS(conn(i,1)) - ZS(conn(i,2)+ns))^2))^2;
end
options = optimset('MaxFunEvals', 1e8, 'MaxIter', 1e6);
Ef = matlabFunction(E, 'Optimize', false, 'Vars', {[XSC; YSC; ZSC; XSN; YSN; ZSN]});
% Ef = matlabFunction(E, 'Optimize', false, 'Vars', {[XS; YS; ZS]});


% Rigid Body Motions
r1 = [ones([N,1]); zeros([N,1]); zeros([N,1])];           % x-translation
r2 = [zeros([N,1]); ones([N,1]); zeros([N,1])];           % y-translation
r3 = [zeros([N,1]); zeros([N,1]); ones([N,1])];           % z-translation
fRz = @(xP, yP) [sqrt(xP.^2 + yP.^2); sqrt(xP.^2 + yP.^2); zeros([length(xP),1])] .*...
                [-sin(atan2(yP, xP)); cos(atan2(yP, xP)); zeros([length(xP),1])];
fRx = @(xP, yP) [zeros([length(xP),1]); sqrt(xP.^2 + yP.^2); sqrt(xP.^2 + yP.^2)] .*...
                [zeros([length(xP),1]); -sin(atan2(yP, xP)); cos(atan2(yP, xP))];
fRy = @(xP, yP) [sqrt(xP.^2 + yP.^2); zeros([length(xP),1]); sqrt(xP.^2 + yP.^2)] .*...
                [-sin(atan2(yP, xP)); zeros([length(xP),1]); cos(atan2(yP, xP))];


%% RK4 Approximation
% Initial Conditions
xC = zeros([N, n]); xC(:,1) = X(1,:)';
yC = zeros([N, n]); yC(:,1) = X(2,:)';
zC = zeros([N, n]); zC(:,1) = X(3,:)';
fC = zeros([1, n]); fC(1) = 0;
Aeq = zeros([6,3*N]);
Aeq(:,[XN, XN+N, XN+2*N]) = eye(6);



fprintf([repmat('.',1,n-1) '\n\n']);
for i = 2:n
    XP1 = xC(:,i-1); YP1 = yC(:,i-1); ZP1 = zC(:,i-1);
    XP = [XP1; YP1; ZP1];
%     P = fmincon(Ef, XP, [], [], Aeq, [xC(XN,i-1)+delS*Us(1,XN)'; yC(XN,i-1)+delS*Us(2,XN)'; zC(XN,i-1)+delS*Us(3,XN)'], [], [], [], options);
    c = [xC(XN,i-1)+delS*Us(1,:)'; yC(XN,i-1)+delS*Us(2,:)'; zC(XN,i-1)+delS*Us(3,:)'];
    P = fminsearch(@(x) Ef([c; x]), XP(nInds'), options);
    fprintf('\b=\n');
    
    xC(XN,i) = xC(XN,i-1)+delS*Us(1,:)';
    yC(XN,i) = yC(XN,i-1)+delS*Us(2,:)';
    zC(XN,i) = zC(XN,i-1)+delS*Us(3,:)';
    
    xC(XNN,i) = P(1:length(XNN))';
    yC(XNN,i) = P([1:length(XNN)]+length(XNN))';
    zC(XNN,i) = P([1:length(XNN)]+2*length(XNN))';
    fC(i) = Ef([xC(XN,i); yC(XN,i); zC(XN,i); xC(XNN,i); yC(XNN,i); zC(XNN,i)]);

%     xC(:,i) = P(1:N)';
%     yC(:,i) = P([1:N]+N)';
%     zC(:,i) = P([1:N]+2*N)';
%     fC(i) = Ef([xC(:,i); yC(:,i); zC(:,i)]);
end

XC = zeros(3, N, n);
XC(1,:,:) = reshape(xC, [1, N, n]);
XC(2,:,:) = reshape(yC, [1, N, n]);
XC(3,:,:) = reshape(zC, [1, N, n]);


%% Plot
% Plot Parameters
ms = 8;        % Marker Size
lw = 2;         % Line Width
ea = .5;        % Edge Transparency
hold on
line([XC(1,conn(:,1),1); XC(1,conn(:,2)+ns,1)],...
     [XC(2,conn(:,1),1); XC(2,conn(:,2)+ns,1)],...
     [XC(3,conn(:,1),1); XC(3,conn(:,2)+ns,1)],...
     'linewidth', lw, 'color', [0 0 0 ea]);
% plot3(XC(1,1:ns,1), XC(2,1:ns,1), XC(3,1:ns,1),...
%     'ro', 'linewidth', ms)
% plot3(XC(1,(1:nu)+ns,1), XC(2,(1:nu)+ns,1), XC(3,(1:nu)+ns,1),...
%     'bo', 'linewidth', ms);
for i = 1:n
    plot3(XC(1,1:ns,i), XC(2,1:ns,i), XC(3,1:ns,i),...
        'r.', 'linewidth', ms/10)
    plot3(XC(1,(1:nu)+ns,i), XC(2,(1:nu)+ns,i), XC(3,(1:nu)+ns,i),...
        'b.', 'linewidth', ms/10);
end
hold off;
end