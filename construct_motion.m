function [Uss, Uus, err] = construct_motion(Xs, Us, Xu, conn, vS, vU, connS, connU)
% Function to generate instantaneous motions
%
% Inputs
% Xs: d x ns matrix of specified node positions
% Us: d x ns x z matrix of z sets of specified node motions as reference
% Xu: d x nu matrix of unspecified node positions
% conn: s x 2 matrix of connections between sp. node i and unsp. node j
% vS: scalar: If positive, v scales specified vectors. 0 to omit plot
% vU: scalar: If positive, v scales unspecified vectors. 0 to omit plot
% connS: k x 2 matrix of connections between specfied nodes i and j
%
% Outputs
% Us: d x ns x m matrix of m non-rigid degrees of freedom of sp. nodes
% Uu: d x nu x m matrix of m non-rigid degrees of freedom of unsp. nodes
% err: 1 x z: error between actual and reconstructed motion

if ~exist('connS', 'var')
    connS = [];
end
if ~exist('connU', 'var')
    connU = [];
end

% Initial Values
d = size(Xs,1);
ns = size(Xs,2);
nu = size(Xu,2);
L = ns + nu;
s1 = size(conn,1);
s2 = size(connS,1);
s3 = size(connU,1);
X = [Xs Xu];
z = size(Us,3);

% Construct rigidity matrix
R = zeros(s1, d*(L));
dInd = ([1:d]-1) * L;
for i = 1:s1
    R(i,conn(i,1) + dInd) = (Xs(:,conn(i,1)) - Xu(:,conn(i,2)))';
    R(i,conn(i,2) + ns + dInd) = (Xu(:,conn(i,2)) - Xs(:,conn(i,1)))';
end
for i = 1:s2
    R(s1+i,connS(i,1) + dInd) = (Xs(:,connS(i,1)) - Xs(:,connS(i,2)))';
    R(s1+i,connS(i,2) + dInd) = (Xs(:,connS(i,2)) - Xs(:,connS(i,1)))';
end
for i = 1:s3
    R(s1+s2+i,connU(i,1) + ns + dInd) = (Xu(:,connU(i,1)) - Xu(:,connU(i,2)))';
    R(s1+s2+i,connU(i,2) + ns + dInd) = (Xu(:,connU(i,2)) - Xu(:,connU(i,1)))';
end

% Find motions
D = null(R);
nSS = size(D,2) - (d*(L) - s1 - s2 - s3);
if(d==2)
    % Rotation function
    fRot = @(xP, yP) [sqrt(xP.^2 + yP.^2); sqrt(xP.^2 + yP.^2)] .*...
                    [-sin(atan2(yP, xP)); cos(atan2(yP, xP))];
    % Matrix of rigid body motions
    T = [ones(1,L), zeros(1,L);...              % x-translation
         zeros(1,L), ones(1,L);...              % y-translation
         fRot(X(1,:)', X(2,:)')'];              % rotation
elseif(d==3)
    fRz = @(xP, yP) [sqrt(xP.^2 + yP.^2); sqrt(xP.^2 + yP.^2); zeros([length(xP),1])] .*...
                    [-sin(atan2(yP, xP)); cos(atan2(yP, xP)); zeros([length(xP),1])];
    fRx = @(xP, yP) [zeros([length(xP),1]); sqrt(xP.^2 + yP.^2); sqrt(xP.^2 + yP.^2)] .*...
                    [zeros([length(xP),1]); -sin(atan2(yP, xP)); cos(atan2(yP, xP))];
    fRy = @(xP, yP) [sqrt(xP.^2 + yP.^2); zeros([length(xP),1]); sqrt(xP.^2 + yP.^2)] .*...
                    [cos(atan2(yP, xP)); zeros([length(xP),1]); -sin(atan2(yP, xP))];
    T = [ones([1,L]), zeros([1,L]), zeros([1,L]);...    % x-translation
         zeros([1,L]), ones([1,L]), zeros([1,L]);...    % y-translation
         zeros([1,L]), zeros([1,L]), ones([1,L]);...    % z-translation
         fRz(X(1,:)', X(2,:)')';...                     % z-rotation
         fRx(X(2,:)', X(3,:)')';...                     % x-rotation
         fRy(X(3,:)', X(1,:)')'];                       % y-rotation
end

% Remove rigid body motions
dIndL = (repmat(dInd', [1,ns]) + [1:ns]);
Uss = zeros(d,ns,max(z,1));
Uus = zeros(d,nu,max(z,1));
if(z==0)
    U = D*null(T*D);
    for j = 1:d
        Uss(j,:) = U((1:ns) + (j-1)*L)';
        Uus(j,:) = U((1:nu) + ns + (j-1)*L)';
    end
    err = 0;
else
    for i = 1:z
        UsP = Us(:,:,i);
        cU = pinv(D([dIndL(:)],:))*UsP(:);
        U = D*cU;
        for j = 1:d
            Uss(j,:,i) = U((1:ns) + (j-1)*L)';
            Uus(j,:,i) = U((1:nu) + ns + (j-1)*L)';
        end
        err(i) = sum(sum(abs(UsP - Uss(:,:,i))));
    end
end

% Error for states of self-stress
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
        for i = 1:max(z,1)
            X = [Xs Xu];
            U = [Uss(:,:,i)*vS Uus(:,:,i)*vU];
            quiver(X(1,:), X(2,:),U(1,:), U(2,:), 0, 'linewidth', lw, 'color', C_SA(i,:));
        end
        line([Xs(1,conn(:,1)); Xu(1,conn(:,2))],...
             [Xs(2,conn(:,1)); Xu(2,conn(:,2))],...
             'linewidth', lw, 'color', [0 0 0 ea]);
        if(s2 > 0)
            line([Xs(1,connS(:,1)); Xs(1,connS(:,2))],...
                 [Xs(2,connS(:,1)); Xs(2,connS(:,2))],...
                 'linewidth', lw, 'color', [0 0 0 ea]);
        end
        plot(Xs(1,:), Xs(2,:), 'o', 'linewidth', ms, 'markersize', ms, 'color', C_SN)
        plot(Xu(1,:), Xu(2,:), 'o', 'linewidth', ms, 'markersize', ms, 'color', C_UN);
        set(gca,'visible',0);
        set(gcf,'color','w');
    elseif(d==3)
        % Spherical point
        [xSp, ySp, zSp] = sphere(20);
        xSp = xSp/10; 
        ySp = ySp/10; 
        zSp = zSp/10; 
        for i = 1:max(z,1)
            X = [Xs Xu];
            U = [Uss(:,:,i)*vS Uus(:,:,i)*vU];
            quiver3(X(1,:), X(2,:), X(3,:), U(1,:), U(2,:), U(3,:),...
                   0, 'linewidth', lw, 'color', C_SA(i,:));
        end
        for i = 1:size(Xu,2)
            s = surf(xSp+Xu(1,i), ySp+Xu(2,i), zSp+Xu(3,i));
            s.FaceColor = C_UN;
            s.EdgeColor = 'none';
        end
        for i = 1:size(Xs,2)
            s = surf(xSp+Xs(1,i), ySp+Xs(2,i), zSp+Xs(3,i));
            s.FaceColor = C_SN;
            s.EdgeColor = 'none';
        end
        line([Xs(1,conn(:,1)); Xu(1,conn(:,2))],...
             [Xs(2,conn(:,1)); Xu(2,conn(:,2))],...
             [Xs(3,conn(:,1)); Xu(3,conn(:,2))],...
             'linewidth', lw, 'color', [0 0 0 ea]);
        if(s2 > 0)
            line([Xs(1,connS(:,1)); Xs(1,connS(:,2))],...
                 [Xs(2,connS(:,1)); Xs(2,connS(:,2))],...
                 [Xs(3,connS(:,1)); Xs(3,connS(:,2))],...
                 'linewidth', lw, 'color', [0 0 0 ea]);
        end
        if(s3 > 0)
            line([Xu(1,connU(:,1)); Xu(1,connU(:,2))],...
                 [Xu(2,connU(:,1)); Xu(2,connU(:,2))],...
                 [Xu(3,connU(:,1)); Xu(3,connU(:,2))],...
                 'linewidth', lw, 'color', [0 0 0 ea]);
        end
        set(gca,'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[],...
                'XTick',[],'YTick',[],'ZTick',[],'box','on','boxstyle','back');
    end
    hold off;
end

end