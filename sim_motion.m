function [XC, fC] = sim_motion(Xs, Xu, conn, delS, n, X0)
% Function for simulating the motion of frames
% Inputs
%       Xs:     2 X ns matrix of coordinates for specified nodes
%       Xu:     2 x nu matrix of coordinates for unspecified nodes
%       conn:   k X 2 matrix of k edges connecting two nodes
%       delS:   Instantaneous length of motion
%       n:      Number of steps to move forward
%       X0      2 X N vector of desired initial direction of motion
% Outputs
%       XC:     2 X N X n matrix of node motions across n time steps
%       fC:     1 X N matrix of motion errors


%% Lengths
ns = size(Xs,2);
nu = size(Xu,2);
X = [Xs Xu];
N = ns + nu;
LVal = sqrt((Xs(1,conn(:,1)) - Xu(1,conn(:,2))).^2 +...
            (Xs(2,conn(:,1)) - Xu(2,conn(:,2))).^2);
        
        
%% Symbolic 
XS = sym('x', [N, 1]); assume(XS, 'real');
YS = sym('y', [N, 1]); assume(YS, 'real');

% Constraints
g = [];
for i = 1:size(conn,1)
    g = [g; (XS(conn(i,1)) - XS(conn(i,2)+ns))^2 + (YS(conn(i,1)) - YS(conn(i,2)+ns))^2];
end
disp('1');
J = jacobian(g, [XS; YS])/2;
disp('2');
Jf = matlabFunction(J, 'Optimize', false, 'Sparse', true);
disp('3');

% Energy
E = (LVal(1) - sqrt((XS(conn(1,1)) - XS(conn(1,2)+ns))^2 + (YS(conn(1,1)) - YS(conn(1,2)+ns))^2))^2;
for i = 2:size(conn,1)
    E = E + (LVal(i) - sqrt((XS(conn(i,1)) - XS(conn(i,2)+ns))^2 + (YS(conn(i,1)) - YS(conn(i,2)+ns))^2))^2;
end
Ef = matlabFunction(E, 'Optimize', false);


% Rigid Body Motions
r1 = [ones([N,1]); zeros([N,1])];           % x-translation
r2 = [zeros([N,1]); ones([N,1])];           % y-translation
fR3 = @(xP, yP) [sqrt(xP.^2 + yP.^2); sqrt(xP.^2 + yP.^2)] .*...
                [-sin(atan2(yP, xP)); cos(atan2(yP, xP))];
r3 = fR3(X(1,:)', X(2,:)');


%% RK4 Approximation
% Initial Conditions
xC = zeros([N, n]); xC(:,1) = X(1,:)';
yC = zeros([N, n]); yC(:,1) = X(2,:)';
fC = zeros([1, n]); fC(1) = 0; X0 = X0';
delX = zeros([2*N, n]); delX(:,1) = X0(:);


fprintf([repmat('.',1,n-1) '\n\n']);
for i = 2:n
    fprintf('\b=\n');
    % Step 1
    XP1 = xC(:,i-1); YP1 = yC(:,i-1);
    XAC = num2cell([XP1; YP1]);
    R = [r1 r2 fR3(XP1, YP1)];
    K = full(Jf(XAC{:}));
    V = null(K);
    k1 = double(V*null(R'*V));
    k1 = sign(k1'*delX(:,i-1)) * k1 / sqrt(k1'*k1);
    
    % Step 2
    XP2 = xC(:,i-1) + delS*k1(1:N)/2; YP2 = yC(:,i-1) + delS*k1([1:N]+N)/2;
    XAC = num2cell([XP2; YP2]);
    R = [r1 r2 fR3(XP2, YP2)];
    K = full(Jf(XAC{:}));
    V = null(K);
    k2 = double(V*null(R'*V));
    k2 = sign(k2'*delX(:,i-1)) * k2 / sqrt(k2'*k2);
    
    % Step 3
    XP3 = xC(:,i-1) + delS*k2(1:N)/2; YP3 = yC(:,i-1) + delS*k2([1:N]+N)/2;
    XAC = num2cell([XP3; YP3]);
    R = [r1 r2 fR3(XP3, YP3)];
    K = full(Jf(XAC{:}));
    V = null(K);
    k3 = double(V*null(R'*V));
    k3 = sign(k3'*delX(:,i-1)) * k3 / sqrt(k3'*k3);
    
    % Step 4
    XP4 = xC(:,i-1) + delS*k3(1:N); YP4 = yC(:,i-1) + delS*k3([1:N]+N);
    XAC = num2cell([XP4; YP4]);
    R = [r1 r2 fR3(XP4, YP4)];
    K = full(Jf(XAC{:}));
    V = null(K);
    k4 = double(V*null(R'*V));
    k4 = sign(k4'*delX(:,i-1)) * k4 / sqrt(k4'*k4);
    
    % Update
    delX(:,i) = k1 + 2*k2 + 2*k3 + k4;
    xC(:,i) = xC(:,i-1) + delS * delX(1:N,i)/6;
    yC(:,i) = yC(:,i-1) + delS * delX([1:N]+N,i)/6;
    XAC = num2cell([xC(:,i); yC(:,i)]);
    fC(i) = Ef(XAC{:});
end

XC = zeros(2, N, n);
XC(1,:,:) = reshape(xC, [1, N, n]);
XC(2,:,:) = reshape(yC, [1, N, n]);


%% Plot
% Plot Parameters
ms = 10;        % Marker Size
lw = 2;         % Line Width
ea = .5;        % Edge Transparency
hold on
line([XC(1,conn(:,1),1); XC(1,conn(:,2)+ns,1)],...
     [XC(2,conn(:,1),1); XC(2,conn(:,2)+ns,1)],...
     'linewidth', lw, 'color', [0 0 0 ea]);
plot(XC(1,1:ns,1), XC(2,1:ns,1), 'ro', 'linewidth', ms)
plot(XC(1,(1:nu)+ns,1), XC(2,(1:nu)+ns,1), 'bo', 'linewidth', ms);
for i = 1:n
    plot(XC(1,1:ns,i), XC(2,1:ns,i), 'r.', 'linewidth', ms)
    plot(XC(1,(1:nu)+ns,i), XC(2,(1:nu)+ns,i), 'b.', 'linewidth', ms);
end
hold off;

end