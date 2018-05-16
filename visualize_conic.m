function [] = visualize_conic(Xs, Us, R, nC, nP, vS, vU)
% Function for visualizing solution spaces in 2D and 3D spatial coordinates
% Solution space dimension must be between 1 and d
%
% Inputs
% Xs: d x n matrix of specified node spatial coordinates
% Us: d x n x z matrix of specified node instantaneous displacements
% R: m x 2 matrix of minimum and maximum spatial ranges to visualize
% nC: m x 1 vector of number of points to sample within range for curve
% nP: Number of points to sample to show displacements
% vS: Scalar: scales the specified arrows by this amount
% vU: Scalar: scales the unspecified arrows by this amount

z = size(Us,3);         % Total number of motions
for j = 1:z
    [Q, W, v0, err] = construct_conic(Xs, Us(:,:,j));
    d = size(Xs,1);         % Dimension of Space
    m = size(Q,1)-1;        % Number of homogeneous variables
    cR = linspace(.5, 1, z);
    
    % 2 Dimensional Space, 1D solution
    if(d==2)
        if(m==3)
            [xx, yy] = meshgrid(linspace(R(1,1),R(1,2),nP),...
                                linspace(R(2,1),R(2,2),nP));
            Cu = [xx(:) yy(:)]';
            Uu = zeros([size(Cu), z]);
            for i = 1:size(Cu,2)
                Uu(:,i,j) = (Cu(:,i)' - Xs')\[Us(:,:,j)'*Cu(:,i) - sum(Xs'.*Us(:,:,j)',2)];
            end
            hold on;
            quiver(Cu(1,:), Cu(2,:), Uu(1,:,j)*vU, Uu(2,:,j)*vU, 0, 'color', [0 0 cR(j)]);
            quiver(Xs(1,:), Xs(2,:), Us(1,:,j)*vS, Us(2,:,j)*vS, 0, 'filled', 'color', [0 cR(j) 0]);
            plot(Xs(1,:), Xs(2,:), 'o', 'linewidth', 10, 'color', 'r');
            hold off;
        elseif(m==2)
            % Change Coordinates
            P = [W([1:d],:) v0(1:d);...
                 zeros(1,d) 1]^-1;
            Q = P'*Q*P;
            % Create Mesh
            [xx, yy] = meshgrid(linspace(R(1,1),R(1,2),nC(1)),...
                                linspace(R(2,1),R(2,2),nC(2)));
            % Evaluate Conic Along Mesh
            F = Q(1,1)*xx.^2 + Q(2,2)*yy.^2 + Q(3,3) +...
                2*Q(2,1)*xx.*yy +...
                2*Q(3,1)*xx + 2*Q(3,2)*yy;
            % Generate contour along conic solution at 0
            hold on;
            C = contour(xx,yy,F,[0 0], 'color', [0 0 cR(j)]);
            Cu = C(:,floor(linspace(ceil(size(C,2)/nP),size(C,2)-1,nP)));
            % Generate displacements along curve
            Uu = zeros([2, nP, z]);
            for i = 1:nP
                Uu(:,i,j) = pinv(Cu(:,i)' - Xs')*[Us(:,:,j)'*Cu(:,i) - sum(Xs'.*Us(:,:,j)',2)];
            end
            quiver(Cu(1,:), Cu(2,:), Uu(1,:,j)*vU, Uu(2,:,j)*vU, 0, 'color', [0 0 cR(j)]);
            quiver(Xs(1,:), Xs(2,:), Us(1,:,j)*vS, Us(2,:,j)*vS, 0, 'filled', 'color', [0 cR(j) 0]);
            plot(Xs(1,:), Xs(2,:), 'o', 'linewidth', 10, 'color', 'r');
            hold off;
        end
    elseif(d==3)
        if(m==4)
            [xx,yy,zz] = meshgrid(linspace(R(1,1),R(1,2),nP),...
                                   linspace(R(2,1),R(2,2),nP),...
                                   linspace(R(3,1),R(3,2),nP));
            Cu = [xx(:) yy(:) zz(:)]';
            % Generate displacements along curve
            Uu = zeros([size(Cu), z]);
            for i = 1:size(Cu,2)
                Uu(:,i,j) = (Cu(:,i)' - Xs')\[Us(:,:,j)'*Cu(:,i) - sum(Xs'.*Us(:,:,j)',2)];
            end
            hold on;
            quiver3(Cu(1,:), Cu(2,:), Cu(3,:), Uu(1,:,j)*vU, Uu(2,:,j)*vU, Uu(3,:,j)*vU, 0, 'color', [0 0 cR(j)]);
            quiver3(Xs(1,:), Xs(2,:), Xs(3,:), Us(1,:,j)*vS, Us(2,:,j)*vS, Us(3,:,j)*vS, 0, 'filled', 'color', [0 cR(j) 0]);
            plot3(Xs(1,:), Xs(2,:), Xs(3,:), 'o', 'markersize', 6, 'linewidth', 10, 'color', 'r');
            hold off;
        elseif(m==3)
            % Change Coordinates
            P = [W([1:d],:) v0(1:d);...
                 zeros(1,d) 1]^-1;
            Q = P'*Q*P;
            % Create Mesh
            [xx, yy, zz] = meshgrid(linspace(R(1,1),R(1,2),nC(1)),...
                                   linspace(R(2,1),R(2,2),nC(2)),...
                                   linspace(R(3,1),R(3,2),nC(3)));
            % Evaluate Conic Along Mesh
            F = Q(1,1)*xx.^2 + Q(2,2)*yy.^2 + Q(3,3)*zz.^2 + Q(4,4)*1 +...
                2*Q(2,1)*xx.*yy + 2*Q(2,3)*yy.*zz + 2*Q(1,3)*xx.*zz +...
                2*Q(4,1)*xx + 2*Q(4,2)*yy + 2*Q(4,3)*zz;
            % Generate contour along conic solution at 0
            hold on;
            fI = isosurface(xx, yy, zz, F, 0);
            p = patch(fI);
            isonormals(xx,yy,zz,F,p)
            C = fI.vertices';
            Cu = C(:,floor(linspace(2,size(C,2)-1,nP)));
            % Generate displacements along curve
            Uu = zeros([3, nP, z]);
            for i = 1:nP
                Uu(:,i,j) = (Cu(:,i)' - Xs')\[Us(:,:,j)'*Cu(:,i) - sum(Xs'.*Us(:,:,j)',2)];
            end
            quiver3(Cu(1,:), Cu(2,:), Cu(3,:), Uu(1,:,j)*vU, Uu(2,:,j)*vU, Uu(3,:,j)*vU, 0, 'color', [0 0 cR(j)]);
            quiver3(Xs(1,:), Xs(2,:), Xs(3,:), Us(1,:,j)*vS, Us(2,:,j)*vS, Us(3,:,j)*vS, 0, 'filled', 'color', [0 cR(j) 0]);
            plot3(Xs(1,:), Xs(2,:), Xs(3,:), 'o', 'markersize', 6, 'linewidth', 10, 'color', 'r');
            hold off;
            p.FaceColor = [0 0 cR(j)];
            p.FaceAlpha = 0.35;
            p.EdgeColor = 'none';
        elseif(m==2)
            % First Coordinate Change to Cartesian Coordinates
            P = [W([1:d],:) v0(1:d)]^-1;
            % Define Plane of Intersection
            p = P(end,:)';
            % Second Coordinate Change to Plane z=1
            P2 = [null(p') p/(p'*p)];
            % Create Mesh
            [xx, yy] = meshgrid(linspace(R(1,1),R(1,2),nC(1)),...
                                linspace(R(2,1),R(2,2),nC(2)));
            XA = P2(1:2,1:2)^-1 * [xx(:) yy(:)]';
            xx = reshape(XA(1,:), size(xx));
            yy = reshape(XA(2,:), size(yy));
            % Evaluate Conic Along Mesh, with added constraint
            Q = P2'*P'*Q*P*P2;
            F = Q(1,1)*xx.^2 + Q(2,2)*yy.^2 + Q(3,3)*1 +...
                2*Q(2,1)*xx.*yy + 2*Q(3,2)*yy + 2*Q(3,1)*xx;
            % Generate contour along conic solution at 1 in second coordinates
            C = contour(xx,yy,F,[0 0], 'color', 'b', 'linestyle', 'none');
            C = [C(:,2:end); ones(1, size(C,2)-1)];
            % Convert second coordinates back to cartesian
            C = P2*C; C = C(:,sum(C > min(R,[],2) & C < max(R,[],2))==3);
            Cu = C(:,floor(linspace(2,size(C,2)-1,nP)));
            % Generate displacements along curve
            Uu = zeros([3, nP, z]);
            for i = 1:nP
                Uu(:,i,j) = (Cu(:,i)' - Xs')\[Us(:,:,j)'*Cu(:,i) - sum(Xs'.*Us(:,:,j)',2)];
            end
            hold on;
            plot3(C(1,:), C(2,:), C(3,:), '');
            quiver3(Cu(1,:), Cu(2,:), Cu(3,:), Uu(1,:,j)*vU, Uu(2,:,j)*vU, Uu(3,:,j)*vU, 0, 'color', [0 0 cR(j)]);
            quiver3(Xs(1,:), Xs(2,:), Xs(3,:), Us(1,:,j)*vS, Us(2,:,j)*vS, Us(3,:,j)*vS, 0, 'filled', 'color', [0 cR(j) 0]);
            plot3(Xs(1,:), Xs(2,:), Xs(3,:), 'o', 'markersize', 6, 'linewidth', 10, 'color', 'r');
            hold off;
        end
    end
end
end