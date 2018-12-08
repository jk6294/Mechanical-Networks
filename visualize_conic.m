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

% Visualization Parameters
LW_SA = 2;                      % Line Width of Specified Arrow
LW_UA = .7;                     % Line Width of Unspecified Arrows
LW_SS = 1;                      % Line Width of Solution Space
MS_SN = 4;                      % Marker Size of Specified Node
MS_UN = 2;                      % Marker Size of Unspecified Node
C_SN = [255 100 100]/255;       % Color of Specified Node
C_SA = [76 187 23;...           % Color of Specified Arrow
        50 255 50]/255;         
C_SS = [100 100 255;...         % Color of Solution Space
        100 200 255]/255;       

z = size(Us,3);                 % Total number of motions
for j = 1:z
    [Q, W, v0, err] = construct_conic(Xs, Us(:,:,j));
    d = size(Xs,1);             % Dimension of Space
    nDOF = size(Q,1)-1;         % Number of homogeneous variables
    
    % 2 Dimensional Space
    if(d==2)
        if(nDOF==3)
            [xx, yy] = meshgrid(linspace(R(1,1),R(1,2),nP),...
                                linspace(R(2,1),R(2,2),nP));
            Cu = [xx(:) yy(:)]';
            Uu = zeros([size(Cu), z]);
            for i = 1:size(Cu,2)
                Uu(:,i,j) = (Cu(:,i)' - Xs')\[Us(:,:,j)'*Cu(:,i) - sum(Xs'.*Us(:,:,j)',2)];
            end
        elseif(nDOF==2)
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
            C = contour(xx,yy,F,[0 0], 'linewidth', LW_SS, 'color', C_SS(j,:));
            Cu = C(:,floor(linspace(ceil(size(C,2)/nP),size(C,2)-1,nP)));
            % Generate displacements along curve
            Uu = zeros([2, nP, z]);
            for i = 1:nP
                Uu(:,i,j) = (Cu(:,i)' - Xs')\[Us(:,:,j)'*Cu(:,i) - sum(Xs'.*Us(:,:,j)',2)];
            end
        end
        % Plot
        pInd = ~(isnan(sum(Uu(:,:,j))) | isinf(sum(Uu(:,:,j))));
        hold on;
        quiver(Cu(1,pInd), Cu(2,pInd), Uu(1,pInd,j)*vU, Uu(2,pInd,j)*vU, 0,...
            'linewidth', LW_UA, 'color', C_SS(j,:));
        plot(Cu(1,pInd), Cu(2,pInd), 'o', 'linewidth', MS_UN, 'markersize', MS_UN, 'color', C_SS(j,:));
        quiver(Xs(1,:), Xs(2,:), Us(1,:,j)*vS, Us(2,:,j)*vS, 0, 'filled',...
            'linewidth', LW_SA, 'color', C_SA(j,:));
        quiver(Xs(1,:), Xs(2,:), Us(1,:,j)*vS, Us(2,:,j)*vS, 0, 'filled',...
            'linewidth', LW_SA/3, 'color', [1 1 1]);
        plot(Xs(1,:), Xs(2,:), 'o', 'linewidth', MS_SN, 'markersize', MS_SN, 'color', C_SN);
        hold off;
        set(gca,'visible',0);
        set(gcf,'color','w');
    elseif(d==3)
        % Spherical point
        [xSp, ySp, zSp] = sphere(20);
        xSp = xSp/10; 
        ySp = ySp/10; 
        zSp = zSp/10; 
        if(nDOF==4)
            [xx,yy,zz] = meshgrid(linspace(R(1,1),R(1,2),nP),...
                                  linspace(R(2,1),R(2,2),nP),...
                                  linspace(R(3,1),R(3,2),nP));
            Cu = [xx(:) yy(:) zz(:)]';
            % Generate displacements along curve
            Uu = zeros([size(Cu), z]);
            for i = 1:size(Cu,2)
                Uu(:,i,j) = (Cu(:,i)' - Xs')\[Us(:,:,j)'*Cu(:,i) - sum(Xs'.*Us(:,:,j)',2)];
            end
        elseif(nDOF==3)
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
            p.FaceColor = C_SS(j,:);
            p.FaceAlpha = 0.5;
            p.EdgeColor = 'none';
        elseif(nDOF==2)
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
            C = contour(xx,yy,F,[0 0], 'linestyle', 'none');
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
            plot3(C(1,:), C(2,:), C(3,:), '-', 'linewidth', LW_SS, 'color', C_SS(j,:));
        end
        % Plot
        quiver3(Cu(1,:), Cu(2,:), Cu(3,:), Uu(1,:,j)*vU, Uu(2,:,j)*vU, Uu(3,:,j)*vU, 0,...
            'filled', 'color', C_SS(j,:), 'linewidth', LW_UA);
        quiver3(Xs(1,:), Xs(2,:), Xs(3,:), Us(1,:,j)*vS, Us(2,:,j)*vS, Us(3,:,j)*vS, 0,...
            'filled', 'color', C_SA(j,:), 'linewidth', LW_SA);
        quiver3(Xs(1,:), Xs(2,:), Xs(3,:), Us(1,:,j)*vS, Us(2,:,j)*vS, Us(3,:,j)*vS, 0,...
            'filled', 'color', [1 1 1], 'linewidth', LW_SA/3);
        for i = 1:size(Cu,2)
            s = surf(xSp*MS_UN/MS_SN+Cu(1,i), ySp*MS_UN/MS_SN+Cu(2,i), zSp*MS_UN/MS_SN+Cu(3,i));
            s.FaceColor = C_SS(j,:);
            s.EdgeColor = 'none';
        end
        for i = 1:size(Xs,2)
            s = surf(xSp+Xs(1,i), ySp+Xs(2,i), zSp+Xs(3,i));
            s.FaceColor = C_SN;
            s.EdgeColor = 'none';
        end
        hold off;
        set(gca,'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[],...
                'XTick',[],'YTick',[],'ZTick',[],'box','on','boxstyle','back');
    end
end
end