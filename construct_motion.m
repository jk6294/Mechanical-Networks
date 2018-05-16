function [Uss, Uus, err] = construct_motion(Xs, Us, Xu, conn, vS, vU)
% Function to generate instantaneous motions
%
% Inputs
% Xs: d x ns matrix of specified node positions
% Us: d x ns x z matrix of z sets of specified node motions as reference
% Xu: d x nu matrix of unspecified node positions
% conn: s x 2 matrix of connections between sp. node i and unsp. node j
% vS: scalar: If positive, v scales specified vectors. <0 to omit plot
% vU: scalar: If positive, v scales unspecified vectors. <0 to omit plot
%
% Outputs
% Us: d x ns x m matrix of m non-rigid degrees of freedom of sp. nodes
% Uu: d x nu x m matrix of m non-rigid degrees of freedom of unsp. nodes
% err: 1 x z: error between actual and reconstructed motion

% Initial Values
d = size(Xs,1);
ns = size(Xs,2);
nu = size(Xu,2);
L = ns + nu;
s = size(conn,1);
X = [Xs Xu];
z = size(Us,3);
cR = linspace(.5, 1, z);

% Construct rigidity matrix
R = zeros(s, d*(L));
dInd = ([1:d]-1) * L;
for i = 1:s
    R(i,conn(i,1) + dInd) = (Xs(:,conn(i,1)) - Xu(:,conn(i,2)))';
    R(i,conn(i,2) + ns + dInd) = (Xu(:,conn(i,2)) - Xs(:,conn(i,1)))';
end

% Find motions
D = null(R);
nSS = size(D,2) - (d*(L) - s);
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
                    [-sin(atan2(yP, xP)); zeros([length(xP),1]); cos(atan2(yP, xP))];
    T = [ones([1,L]), zeros([1,L]), zeros([1,L]);...    % x-translation
         zeros([1,L]), ones([1,L]), zeros([1,L]);...    % y-translation
         zeros([1,L]), zeros([1,L]), ones([1,L]);...    % z-translation
         fRz(X(1,:)', X(2,:)')';...                     % z-rotation
         fRx(X(2,:)', X(3,:)')';...                     % x-rotation
         fRy(X(3,:)', X(1,:)')'];                       % y-rotation
end

% Remove rigid body motions
% U = D * null(T*D);
dIndL = repmat(dInd', [1,ns]) + [1:ns];
Uss = zeros(d,ns,z);
Uus = zeros(d,nu,z);
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

% Error for states of self-stress
if(nSS == 0)
    disp('No states of self-stress; motions should be finitely deformable');
else
    disp('States of self-stress; motions may not be finitely deformable');
end

% Visualize
% Plot Parameters
ms = 10;        % Marker Size
lw = 2;         % Line Width
ea = .5;        % Edge Transparency
if(vS>=0 && vU >= 0)
    hold on;
    if(d==2)
        for i = 1:z
            quiver(Xs(1,:), Xs(2,:), Uss(1,:,i)*vS, Uss(2,:,i)*vS,...
                   0, 'linewidth', lw, 'color', [0, cR(i), 0]);
            quiver(Xu(1,:), Xu(2,:), Uus(1,:,i)*vU, Uus(2,:,i)*vU,...
                   0, 'linewidth', lw, 'color', [0, cR(i), 0]);
        end
        line([Xs(1,conn(:,1)); Xu(1,conn(:,2))],...
             [Xs(2,conn(:,1)); Xu(2,conn(:,2))],...
             'linewidth', lw, 'color', [0 0 0 ea]);
        plot(Xs(1,:), Xs(2,:), 'ro', 'linewidth', ms)
        plot(Xu(1,:), Xu(2,:), 'bo', 'linewidth', ms);
    elseif(d==3)
        for i = 1:z
            quiver3(Xs(1,:), Xs(2,:), Xs(3,:), Uss(1,:,i)*vS, Uss(2,:,i)*vS, Uss(3,:,i)*vS,...
                   0, 'linewidth', lw, 'color', [0, cR(i), 0]);
            quiver3(Xu(1,:), Xu(2,:), Xu(3,:), Uus(1,:,i)*vU, Uus(2,:,i)*vU, Uus(3,:,i)*vU,...
                   0, 'linewidth', lw, 'color', [0, cR(i), 0]);
        end
        plot3(Xs(1,:), Xs(2,:), Xs(3,:), 'ro', 'linewidth', ms)
        plot3(Xu(1,:), Xu(2,:), Xu(3,:), 'bo', 'linewidth', ms);
        line([Xs(1,conn(:,1)); Xu(1,conn(:,2))],...
             [Xs(2,conn(:,1)); Xu(2,conn(:,2))],...
             [Xs(3,conn(:,1)); Xu(3,conn(:,2))],...
             'linewidth', lw, 'color', [0 0 0 ea]);
    end
    hold off;
end

end