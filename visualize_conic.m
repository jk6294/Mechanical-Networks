function [] = visualize_conic(Xs, Xsdot, R, nC, nP, vS, vU)
% Function for visualizing solution spaces in 2D and 3D spatial coordinates
% Solution space must be 1 or 2 dimensional
%
% Inputs
% Xs:    d x n matrix of specified node spatial coordinates
% Xsdot: d x n x z matrix of specified node motions
% R:     m x 2 matrix of minimum and maximum spatial ranges to visualize
% nC:    m x 1 vector of number of points to sample within range for curve
% nP:    1 x 1 Number of points to sample to show displacements
% vS:    1 x 1 scalar: scales the specified arrows by this amount
% vU:    1 x 1 scalar: scales the unspecified arrows by this amount

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
z = size(Xsdot,3);              % Total number of motions

% Draw solution spaces for each motion
for j = 1:z
    [Q, W, v0, err] = construct_conic(Xs, Xsdot(:,:,j), 0);
    d = size(Xs,1);             % Dimension of Space
    nDOF = size(Q,1)-1;         % Number of homogeneous variables
    
    % 2 Dimensional Embedding
    if(d==2)
        % Evenly sample points in 2D space and find motions
        if(nDOF==3)
            % Sample
            [xx, yy] = meshgrid(linspace(R(1,1),R(1,2),nP),...
                                linspace(R(2,1),R(2,2),nP));
            Cu = [xx(:) yy(:)]';                                
            Uu = zeros([size(Cu), z]);
            % Solve for motions at each point
            for i = 1:size(Cu,2)
                Uu(:,i,j) = (Cu(:,i)' - Xs') \ [Xsdot(:,:,j)'*Cu(:,i)...
                                              - sum(Xs'.*Xsdot(:,:,j)',2)];
            end
        % Sample conic in 2 dimensions and find zero-crossings
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
            % Generate motions along curve
            Uu = zeros([2, nP, z]);
            for i = 1:nP
                Uu(:,i,j) = (Cu(:,i)' - Xs')\[Xsdot(:,:,j)'*Cu(:,i)-sum(Xs'.*Xsdot(:,:,j)',2)];
            end
        end
        
        % Plot
        pInd = ~(isnan(sum(Uu(:,:,j))) | isinf(sum(Uu(:,:,j))));
        hold on;
        % Unspecified Arrows
        quiver(Cu(1,pInd), Cu(2,pInd), Uu(1,pInd,j)*vU, Uu(2,pInd,j)*vU, 0,...
            'linewidth', LW_UA, 'color', C_SS(j,:));
        % Unspecified Nodes
        plot(Cu(1,pInd), Cu(2,pInd), 'o', 'linewidth', MS_UN, 'markersize', MS_UN, 'color', C_SS(j,:));
        % Specified Arrows
        quiver(Xs(1,:), Xs(2,:), Xsdot(1,:,j)*vS, Xsdot(2,:,j)*vS, 0, 'filled',...
            'linewidth', LW_SA, 'color', C_SA(j,:));
        quiver(Xs(1,:), Xs(2,:), Xsdot(1,:,j)*vS, Xsdot(2,:,j)*vS, 0, 'filled',...
            'linewidth', LW_SA/3, 'color', [1 1 1]);
        % Specified Nodes
        plot(Xs(1,:), Xs(2,:), 'o', 'linewidth', MS_SN, 'markersize', MS_SN, 'color', C_SN);
        hold off;
        set(gca,'visible',0);
        set(gcf,'color','w');
        
    % 3 Dimensional Embedding
    elseif(d==3)
        % Spherical point for visualization
        [xSp, ySp, zSp] = sphere(20);
        xSp = xSp/10; 
        ySp = ySp/10; 
        zSp = zSp/10;
        if(nDOF==3)
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
                Uu(:,i,j) = (Cu(:,i)' - Xs') \ [Xsdot(:,:,j)'*Cu(:,i)-sum(Xs'.*Xsdot(:,:,j)',2)];
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
                Uu(:,i,j) = (Cu(:,i)' - Xs')\[Xsdot(:,:,j)'*Cu(:,i) - sum(Xs'.*Xsdot(:,:,j)',2)];
            end
            hold on;
            plot3(C(1,:), C(2,:), C(3,:), '-', 'linewidth', LW_SS, 'color', C_SS(j,:));
        end
        
        % Plot
        % Unspecified Arrows
        quiver3(Cu(1,:), Cu(2,:), Cu(3,:), Uu(1,:,j)*vU, Uu(2,:,j)*vU, Uu(3,:,j)*vU, 0,...
            'filled', 'color', C_SS(j,:), 'linewidth', LW_UA);
        % Specified Arrows
        quiver3(Xs(1,:), Xs(2,:), Xs(3,:), Xsdot(1,:,j)*vS, Xsdot(2,:,j)*vS, Xsdot(3,:,j)*vS, 0,...
            'filled', 'color', C_SA(j,:), 'linewidth', LW_SA);
        quiver3(Xs(1,:), Xs(2,:), Xs(3,:), Xsdot(1,:,j)*vS, Xsdot(2,:,j)*vS, Xsdot(3,:,j)*vS, 0,...
            'filled', 'color', [1 1 1], 'linewidth', LW_SA/3);
        % Unspecified Nodes
        for i = 1:size(Cu,2)
            s = surf(xSp*MS_UN/MS_SN+Cu(1,i), ySp*MS_UN/MS_SN+Cu(2,i), zSp*MS_UN/MS_SN+Cu(3,i));
            s.FaceColor = C_SS(j,:);
            s.EdgeColor = 'none';
        end
        % Specified Nodes
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