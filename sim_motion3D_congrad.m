function [XMot, fC] = sim_motion3D_congrad(Xs, Xu, conn, LVal, delS, n, XN, XF, pV, kL)
% Function for simulating the motion of elastic frames in 2 or 3 dimensions
% with energy calculation.
%
% Inputs
% Xs:       3 x ns      matrix of coordinates for specified nodes
% Xu:       3 x nu      matrix of coordinates for unspecified nodes
% conn:     k x 2       matrix of k edges connecting node i to node j
% LVal:     k x 1       vector of equilibrium edge lengths
% delS:     1 x 1       instantaneous length of motion
% n:        1 x 1       number of simulation time steps
% XN:       1 x s       vector of nodes to set position boundary conditions
% XF:       d x s       matrix of final positions for boundary nodes
% pV:       1 x 1       scalar to determine plotting. 0 for no plot
% kL:       k x 1       vector of spring constants for all edges

% Default to no plot, and spring constants of 1, if not specified
if ~exist('pV', 'var')
    pV = 0;
end
if ~exist('kL', 'var')
    kL = ones(length(LVal), 1);
end


%% Parameters
ns = size(Xs,2);        % Number of specified nodes
nu = size(Xu,2);        % Number of unspecified nodes
d = size(Xs,1);         % Dimension of embedding
X = [Xs Xu]';           % Combined node positions
N = ns + nu;            % Number of total nodes
Us = (XF - Xs(:,XN))';  % Displacement of boundary nodes


%% Symbolic
if(d==2)
    sInds = [XN, XN+N];               % Indices of controlled coordinates
    nInds = setdiff(1:(2*N), sInds);  % Indices of uncontrolled coordinates
    % Create Symbolic Variables
    XS = sym('x', [N, 1]); assume(XS, 'real');
    YS = sym('y', [N, 1]); assume(YS, 'real');
    XNN = setdiff(1:N, XN);
    XV = [XS(XNN); YS(XNN)];
    % Energy
    E = sum(kL.*(LVal - sqrt((XS(conn(:,1)) - XS(conn(:,2))).^2 + (YS(conn(:,1)) - YS(conn(:,2))).^2)).^2);
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
    E = sum(kL.*(LVal - sqrt((XS(conn(:,1)) - XS(conn(:,2))).^2 + (YS(conn(:,1)) - YS(conn(:,2))).^2 + (ZS(conn(:,1)) - ZS(conn(:,2))).^2)).^2);
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

% Print Total Progress
fprintf([repmat('.',1,n-1) '\n\n']);
for i = 2:n
    % Indicate Progress
    fprintf('\b=\n');
    
    % Initialize current step from past step
    XP(:,i) = XP(:,i-1);
    % Set displaced control nodes as boundary conditions
    XP(sInds,i) = XP(sInds,i) + delS*Us(:);
    
    % Initialize Gradient
    EJfe = ones(d*N,1);
    
    % Newton-Raphson until gradient norm below bound
    while(norm(EJfe) > 1e-12)
        % Convert to cell
        XPS = num2cell(XP(:,i));
        % Evaluate Gradient
        EJfe = EJf(XPS{:});
        % Newton Method
        EJ2fe = EJ2f(XPS{:});
        XP(nInds,i) = XP(nInds,i) - pinv(EJ2fe)*EJfe;
    end
    
    % Evaluate Energy
    fC(i) = Ef(XPS{:});
end

% Reshape simulationr esults
XMot = zeros(d, N, n);
for i = 1:d
    XMot(i,:,:) = XP((1:N)+(i-1)*N,:);
end


%% Plot
% Plot Parameters
% Energy
C = parula(100);
EInd = discretize(fC,100);

if(pV ~= 0)
    X = X';
    hold on;
    % Plot Parameters
    ms = 4;         % Marker Size
    lw = 1;         % Line Width
    ea = .5;        % Edge Transparency
    C_SN = [255 100 100]/255;       % Color of Specified Node
    C_UN = [100 100 255]/255;       % Color of Unspecified Node
    
    if(d==2)
        % Traversal through time
        for i = 1:n
            plot(XMot(1,1:ns,i), XMot(2,1:ns,i), '.', 'linewidth', ms, 'markersize', ms, 'color', C(EInd(i),:))
            plot(XMot(1,(1:nu)+ns,i), XMot(2,(1:nu)+ns,i), '.', 'linewidth', ms, 'markersize', ms, 'color', C(EInd(i),:));
        end
        
        % Edges
        line([X(1,conn(:,1)); X(1,conn(:,2))],...
             [X(2,conn(:,1)); X(2,conn(:,2))],...
             'linewidth', lw, 'color', [0 0 0 ea]);
        % Specified Nodes
        plot(Xs(1,:), Xs(2,:), 'o', 'linewidth', ms, 'markersize', ms, 'color', C_SN)
        % Unspecified Nodes
        plot(Xu(1,:), Xu(2,:), 'o', 'linewidth', ms, 'markersize', ms, 'color', C_UN(1,:));
        set(gca,'visible',0);
        set(gcf,'color','w');
    
    elseif(d==3)
        % Traversal through time
        for i = 1:n
            plot3(XMot(1,1:ns,i), XMot(2,1:ns,i), XMot(3,1:ns,i),...
                '.', 'linewidth', ms/10, 'markersize', ms/10, 'color', C(EInd(i),:));
            plot3(XMot(1,(1:nu)+ns,i), XMot(2,(1:nu)+ns,i), XMot(3,(1:nu)+ns,i),...
                '.', 'linewidth', ms/10, 'markersize', ms/10, 'color', C(EInd(i),:));
        end
        % Spherical point
        [xSp, ySp, zSp] = sphere(20);
        xSp = xSp/10; 
        ySp = ySp/10; 
        zSp = zSp/10; 
        % Edges
        line([X(1,conn(:,1)); X(1,conn(:,2))],...
             [X(2,conn(:,1)); X(2,conn(:,2))],...
             [X(3,conn(:,1)); X(3,conn(:,2))],...
             'linewidth', lw, 'color', [0 0 0 ea]);
        % Unspecified Nodes
        for i = 1:size(Xu,2)
            s = surf(xSp+Xu(1,i), ySp+Xu(2,i), zSp+Xu(3,i));
            s.FaceColor = C_UN(1,:);
            s.EdgeColor = 'none';
        end
        % Specified Nodes
        for i = 1:size(Xs,2)
            s = surf(xSp+Xs(1,i), ySp+Xs(2,i), zSp+Xs(3,i));
            s.FaceColor = C_SN;
            s.EdgeColor = 'none';
        end
        set(gca,'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[],...
                'XTick',[],'YTick',[],'ZTick',[],'box','on','boxstyle','back');
    end
    hold off;
end


end