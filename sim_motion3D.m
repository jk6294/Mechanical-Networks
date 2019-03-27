function [XC, fC] = sim_motion3D(Xs, Xu, conn, delS, n, X0, pV)
% Function for simulating the motion of frames in 3-dimensions
% with one degree of freedom
%
% Inputs
% Xs:     3 x ns        matrix of coordinates for specified nodes
% Xu:     3 x nu        matrix of coordinates for unspecified nodes
% conn:   k x 2         matrix of k edges connecting node i to node j
% delS:   1 x 1         instantaneous length of motion
% n:      1 x 1         number of simulation time steps
% X0:     3 x N         vector of desired initial direction of motion
% pV:     1 x 1         scalar to determine plotting. 0 for no plot
%
% Outputs
% XC:     3 X N X n     matrix of node motions across n time steps
% fC:     1 X n         matrix of motion errors


%% Lengths
ns = size(Xs,2);        % Number of specified nodes
nu = size(Xu,2);        % Number of unspecified nodes
X = [Xs Xu];            % Combined node positions
N = ns + nu;            % Number of total nodes
% All edge lengths
LVal = sqrt((X(1,conn(:,1)) - X(1,conn(:,2))).^2 +...
            (X(2,conn(:,1)) - X(2,conn(:,2))).^2 +...
            (X(3,conn(:,1)) - X(3,conn(:,2))).^2)';
        
        
%% Symbolic 
XS = sym('x', [N, 1]); assume(XS, 'real');
YS = sym('y', [N, 1]); assume(YS, 'real');
ZS = sym('z', [N, 1]); assume(ZS, 'real');

% Construct rigidity matrix
I = zeros(size(conn,1), N);
for i = 1:size(conn,1)
    I(i,conn(i,1)) =  1;
    I(i,conn(i,2)) = -1;
end
J = [(XS(conn(:,1)) - XS(conn(:,2))) .* I,...
     (YS(conn(:,1)) - YS(conn(:,2))) .* I,...
     (ZS(conn(:,1)) - ZS(conn(:,2))) .* I];
Jf = matlabFunction(J, 'Optimize', false, 'Sparse', true);

% Energy
E = sum((LVal - sqrt((XS(conn(:,1)) - XS(conn(:,2))).^2 + (YS(conn(:,1)) - YS(conn(:,2))).^2 + (ZS(conn(:,1)) - ZS(conn(:,2))).^2)).^2);
Ef = matlabFunction(E, 'Optimize', false);

% Rigid Body Motions
r1 = [ones([N,1]); zeros([N,1]); zeros([N,1])];           % x-translation
r2 = [zeros([N,1]); ones([N,1]); zeros([N,1])];           % y-translation
r3 = [zeros([N,1]); zeros([N,1]); ones([N,1])];           % z-translation
fRz = @(xP, yP) [sqrt(xP.^2 + yP.^2); sqrt(xP.^2 + yP.^2); zeros([length(xP),1])] .*...
                [-sin(atan2(yP, xP)); cos(atan2(yP, xP)); zeros([length(xP),1])];
fRx = @(xP, yP) [zeros([length(xP),1]); sqrt(xP.^2 + yP.^2); sqrt(xP.^2 + yP.^2)] .*...
                [zeros([length(xP),1]); -sin(atan2(yP, xP)); cos(atan2(yP, xP))];
fRy = @(xP, yP) [sqrt(xP.^2 + yP.^2); zeros([length(xP),1]); sqrt(xP.^2 + yP.^2)] .*...
                [cos(atan2(yP, xP)); zeros([length(xP),1]); -sin(atan2(yP, xP))];


%% RK4 Approximation
% Initial Conditions
xC = zeros([N, n]); xC(:,1) = X(1,:)';
yC = zeros([N, n]); yC(:,1) = X(2,:)';
zC = zeros([N, n]); zC(:,1) = X(3,:)';
fC = zeros([1, n]); fC(1) = 0; X0 = X0';
delX = zeros([3*N, n]); delX(:,1) = X0(:);

% Print Total Progress
fprintf([repmat('.',1,n-1) '\n\n']);
for i = 2:n
    % Indicate Progress
    fprintf('\b=\n');
    
    % Step 1
    XP1 = xC(:,i-1); YP1 = yC(:,i-1); ZP1 = zC(:,i-1);
    XAC = num2cell([XP1; YP1; ZP1]);
    R = [r1 r2 r3 fRz(XP1, YP1) fRx(YP1, ZP1) fRy(ZP1, XP1)];
    K = full(Jf(XAC{:}));
    V = null(K);
    k1 = double(V*null(R'*V));
    k1 = sign(k1'*delX(:,i-1)) * k1 / sqrt(k1'*k1);
    
    % Step 2
    XP2 = xC(:,i-1) + delS*k1(1:N)/2; 
    YP2 = yC(:,i-1) + delS*k1((1:N)+N)/2;
    ZP2 = zC(:,i-1) + delS*k1((1:N)+2*N)/2;
    XAC = num2cell([XP2; YP2; ZP2]);
    R = [r1 r2 r3 fRz(XP2, YP2) fRx(YP2, ZP2) fRy(ZP2, XP2)];
    K = full(Jf(XAC{:}));
    V = null(K);
    k2 = double(V*null(R'*V));
    k2 = sign(k2'*delX(:,i-1)) * k2 / sqrt(k2'*k2);
    
    % Step 3
    XP3 = xC(:,i-1) + delS*k2(1:N)/2; 
    YP3 = yC(:,i-1) + delS*k2((1:N)+N)/2;
    ZP3 = zC(:,i-1) + delS*k2((1:N)+2*N)/2;
    XAC = num2cell([XP3; YP3; ZP3]);
    R = [r1 r2 r3 fRz(XP3, YP3) fRx(YP3, ZP3) fRy(ZP3, XP3)];
    K = full(Jf(XAC{:}));
    V = null(K);
    k3 = double(V*null(R'*V));
    k3 = sign(k3'*delX(:,i-1)) * k3 / sqrt(k3'*k3);
    
    % Step 4
    XP4 = xC(:,i-1) + delS*k3(1:N); 
    YP4 = yC(:,i-1) + delS*k3((1:N)+N);
    ZP4 = zC(:,i-1) + delS*k3((1:N)+2*N);
    XAC = num2cell([XP4; YP4; ZP4]);
    R = [r1 r2 r3 fRz(XP4, YP4) fRx(YP4, ZP4) fRy(ZP4, XP4)];
    K = full(Jf(XAC{:}));
    V = null(K);
    k4 = double(V*null(R'*V));
    k4 = sign(k4'*delX(:,i-1)) * k4 / sqrt(k4'*k4);
    
    % Update
    delX(:,i) = k1 + 2*k2 + 2*k3 + k4;
    xC(:,i) = xC(:,i-1) + delS * delX(1:N,i)/6;
    yC(:,i) = yC(:,i-1) + delS * delX((1:N)+N,i)/6;
    zC(:,i) = zC(:,i-1) + delS * delX((1:N)+2*N,i)/6;
    XAC = num2cell([xC(:,i); yC(:,i); zC(:,i)]);
    fC(i) = Ef(XAC{:});
end

XC = zeros(2, N, n);
XC(1,:,:) = reshape(xC, [1, N, n]);
XC(2,:,:) = reshape(yC, [1, N, n]);
XC(3,:,:) = reshape(zC, [1, N, n]);


%% Plot
% Plot Parameters
ms = 10;        % Marker Size
lw = 2;         % Line Width
ea = .5;        % Edge Transparency
if(pV ~= 0)
    hold on
    % Edges
    line([XC(1,conn(:,1),1); XC(1,conn(:,2),1)],...
         [XC(2,conn(:,1),1); XC(2,conn(:,2),1)],...
         [XC(3,conn(:,1),1); XC(3,conn(:,2),1)],...
         'linewidth', lw, 'color', [0 0 0 ea]);
    % Specified Nodes
    plot3(XC(1,1:ns,1), XC(2,1:ns,1), XC(3,1:ns,1),...
        'ro', 'linewidth', ms)
    % Unspecified Nodes
    plot3(XC(1,(1:nu)+ns,1), XC(2,(1:nu)+ns,1), XC(3,(1:nu)+ns,1),...
        'bo', 'linewidth', ms);
    % Node Movements Over Time
    for i = 1:n
        plot3(XC(1,1:ns,i), XC(2,1:ns,i), XC(3,1:ns,i),...
            'r.', 'linewidth', ms/10)
        plot3(XC(1,(1:nu)+ns,i), XC(2,(1:nu)+ns,i), XC(3,(1:nu)+ns,i),...
            'b.', 'linewidth', ms/10);
    end
    hold off;
end
end