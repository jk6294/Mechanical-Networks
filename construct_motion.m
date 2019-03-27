function [Uss, Uus, err] = construct_motion(Xs, Xsdot, Xu, conn, vS, vU)
% Function to generate instantaneous motions
%
% Inputs
% Xs:       d x ns  matrix of specified node positions
% Xsdot:    d x ns x z matrix of z sets of specified node motions
% Xu:       d x nu  matrix of unspecified node positions
% conn:     s x 2   matrix of connections between node i and node j
% vS:       1 x 1   If positive, v scales specified vectors. 0 to omit
% vU:       1 x 1   If positive, v scales unspecified vectors. 0 to omit
%
% Outputs
% Uss: d x ns x m matrix of m non-rigid degrees of freedom of sp. nodes
% Uus: d x nu x m matrix of m non-rigid degrees of freedom of unsp. nodes
% err: 1 x z: error between actual and reconstructed motion

% Initial Values
d = size(Xs,1);         % Dimension of embedding
ns = size(Xs,2);        % Number specified nodes
nu = size(Xu,2);        % Number unspecified nodes
N = ns + nu;            % Total number of nodes
X = [Xs Xu];            % Combined sp. and unsp. node positions
z = size(Xsdot,3);      % Number of motions as initial guesses

% Construct rigidity matrix
R = zeros(size(conn,1), d*N);
dInd = ([1:d]-1) * N;
for i = 1:size(conn,1)
    R(i,conn(i,1) + dInd) = (X(:,conn(i,1)) - X(:,conn(i,2)))';
    R(i,conn(i,2) + dInd) = (X(:,conn(i,2)) - X(:,conn(i,1)))';
end

% Find motions
D = null(R);                          % Dimension of infinitesimal motions
nSS = size(D,2) - d*N + size(conn,1); % Number of states of self-stress
if(d==2)
    % Rotation function
    fRot = @(xP, yP) [sqrt(xP.^2 + yP.^2); sqrt(xP.^2 + yP.^2)] .*...
                     [-sin(atan2(yP, xP)); cos(atan2(yP, xP))];
    % Matrix of rigid body motions
    T = [ones(1,N), zeros(1,N);...              % x-translation
         zeros(1,N), ones(1,N);...              % y-translation
         fRot(X(1,:)', X(2,:)')'];              % rotation
elseif(d==3)
    % Rotation Functions
    fRz = @(xP, yP) [sqrt(xP.^2 + yP.^2); sqrt(xP.^2 + yP.^2); zeros([length(xP),1])] .*...
                    [-sin(atan2(yP, xP)); cos(atan2(yP, xP)); zeros([length(xP),1])];
    fRx = @(xP, yP) [zeros([length(xP),1]); sqrt(xP.^2 + yP.^2); sqrt(xP.^2 + yP.^2)] .*...
                    [zeros([length(xP),1]); -sin(atan2(yP, xP)); cos(atan2(yP, xP))];
    fRy = @(xP, yP) [sqrt(xP.^2 + yP.^2); zeros([length(xP),1]); sqrt(xP.^2 + yP.^2)] .*...
                    [cos(atan2(yP, xP)); zeros([length(xP),1]); -sin(atan2(yP, xP))];
    % Matrix of rigid body motions
    T = [ones([1,N]), zeros([1,N]), zeros([1,N]);...    % x-translation
         zeros([1,N]), ones([1,N]), zeros([1,N]);...    % y-translation
         zeros([1,N]), zeros([1,N]), ones([1,N]);...    % z-translation
         fRz(X(1,:)', X(2,:)')';...                     % z-rotation
         fRx(X(2,:)', X(3,:)')';...                     % x-rotation
         fRy(X(3,:)', X(1,:)')'];                       % y-rotation
end

% Remove rigid body motions
dIndL = (repmat(dInd', [1,ns]) + [1:ns]);
Uss = zeros(d,ns,max(z,1));
Uus = zeros(d,nu,max(z,1));
% If no guess provided, raw reconstruction
if(z==0)
    % Project infinitesimal motions onto rigid body motions, and remove
    U = D*null(T*D);
    for j = 1:d
        Uss(j,:) = U((1:ns) + (j-1)*N)';
        Uus(j,:) = U((1:nu) + ns + (j-1)*N)';
    end
    err = 0;
% If guess provided, compute closest projection onto guess
else
    err = zeros(z,1);
    % Iterate across each guessed motion
    for i = 1:z
        % Minimum norm reconstruction of guessed motion using nullspace
        UsP = Xsdot(:,:,i);
        cU = pinv(D([dIndL(:)],:))*UsP(:);
        U = D*cU;
        for j = 1:d
            Uss(j,:,i) = U((1:ns) + (j-1)*N)';
            Uus(j,:,i) = U((1:nu) + ns + (j-1)*N)';
        end
        % Euclidean distance of reconstruction to guess
        err(i) = sum(sum((UsP - Uss(:,:,i)).^2));
    end
end

% Error warning for states of self-stress
if(nSS == 0)
    disp('No states of self-stress; motions should be finitely deformable');
else
    disp('States of self-stress; motions may not be finitely deformable');
end

% Visualize
% Plot Parameters
ms = 4;                         % Marker Size
lw = 1;                         % Line Width
ea = .5;                        % Edge Transparency
C_SN = [255 100 100]/255;       % Color of Specified Node
C_SA = [76 187 23;...           % Color of Specified Arrow
        50 255 50]/255;         
C_UN = [100 100 255]/255;       % Color of Solution Space   
if(vS~=0 && vU ~= 0)
    hold on;
    if(d==2)
        % Motions
        for i = 1:max(z,1)
            U = [Uss(:,:,i)*vS Uus(:,:,i)*vU];
            quiver(X(1,:), X(2,:),U(1,:), U(2,:), 0, 'linewidth', lw, 'color', C_SA(i,:));
        end
        % Edges
        line([X(1,conn(:,1)); X(1,conn(:,2))],...
             [X(2,conn(:,1)); X(2,conn(:,2))],...
             'linewidth', lw, 'color', [0 0 0 ea]);
        % Specified Nodes
        plot(Xs(1,:), Xs(2,:), 'o', 'linewidth', ms, 'markersize', ms, 'color', C_SN)
        % Unspecified Nodes
        plot(Xu(1,:), Xu(2,:), 'o', 'linewidth', ms, 'markersize', ms, 'color', C_UN);
        set(gca,'visible',0);
        set(gcf,'color','w');
    elseif(d==3)
        % Spherical point
        [xSp, ySp, zSp] = sphere(20);
        xSp = xSp/10; 
        ySp = ySp/10; 
        zSp = zSp/10; 
        % Motions
        for i = 1:max(z,1)
            U = [Uss(:,:,i)*vS Uus(:,:,i)*vU];
            quiver3(X(1,:), X(2,:), X(3,:), U(1,:), U(2,:), U(3,:),...
                   0, 'linewidth', lw, 'color', C_SA(i,:));
        end
        % Unspecified Nodes
        for i = 1:size(Xu,2)
            s = surf(xSp+Xu(1,i), ySp+Xu(2,i), zSp+Xu(3,i));
            s.FaceColor = C_UN;
            s.EdgeColor = 'none';
        end
        % Specified Nodes
        for i = 1:size(Xs,2)
            s = surf(xSp+Xs(1,i), ySp+Xs(2,i), zSp+Xs(3,i));
            s.FaceColor = C_SN;
            s.EdgeColor = 'none';
        end
        % Edges
        line([X(1,conn(:,1)); X(1,conn(:,2))],...
             [X(2,conn(:,1)); X(2,conn(:,2))],...
             [X(3,conn(:,1)); X(3,conn(:,2))],...
             'linewidth', lw, 'color', [0 0 0 ea]);
        set(gca,'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[],...
                'XTick',[],'YTick',[],'ZTick',[],'box','on','boxstyle','back');
    end
    hold off;
end
end