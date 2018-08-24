function [] = visualize_conic_finite(X0, XF, R, nC, nP, vS, vU)
% Function for visualizing solution spaces in 2D and 3D spatial coordinates
% Solution space dimension must be between 1 and d
%
% Inputs
% X0: d x n matrix of specified node initial spatial coordinates
% XF: d x n x j matrix of specified node final spatial coordinates
% R: m x 2 matrix of minimum and maximum spatial ranges to visualize
% nC: m x 1 vector of number of points to sample within range for curve
% nP: Number of points to sample to show displacements
% vS: Scalar: scales the specified arrows by this amount
% vU: Scalar: scales the unspecified arrows by this amount

% Visualization Parameters
LW_SA = 5;                      % Line Width of Specified Arrow
LW_UA = 1;                      % Line Width of Unspecified Arrows
LW_SS = 2;                      % Line Width of Solution Space
MS_SN = 8;                      % Marker Size of Specified Node
MS_UN = 3;                      % Marker Size of Unspecified Node
C_SN = [255 100 100]/255;       % Color of Specified Node
C_SA = [76 187 23;...           % Color of Specified Arrow
        50 255 50]/255;         
C_SS = [100 100 255;...         % Color of Solution Space
        100 200 255]/255;       
    
z = size(XF,3);                 % Total number of motions
Us = XF - X0;
for j = 1:z
    [Q, W, v0, err] = construct_conic_finite(X0, XF(:,:,j));
    d = size(X0,1);             % Dimension of Space
    nDOF = size(Q,1)-1;         % Number of homogeneous variables
    cR = linspace(.5, 1, z);
    
    % 2 Dimensional Space
    if(d==2)
        if(nDOF==3)
            [xx, yy] = meshgrid(linspace(R(1,1),R(1,2),nP),...
                                linspace(R(2,1),R(2,2),nP));
            Cu = [xx(:) yy(:)]';
            % Lengths
            l = [sum((X0(:,1)-Cu).^2); sum((X0(:,2)-Cu).^2)];
            b = [sum(XF(:,1,j).^2); sum(XF(:,2,j).^2)] - l;
            A = 2*XF(:,:,j)';
            vs = [A -ones(d,1)]\b;
            ns = null([A -ones(d,1)]);
            O = eye(d+1); O(d+1,d+1) = 0;
            p = [zeros(1,d) 1];
            % Coefficients
            alp = ns'*O*ns;
            bet = (2*vs'*O - p)*ns;
            gam = sum((vs'*O - p).*vs',2);
            a1 = (-bet + sqrt(bet.^2 - 4*alp.*gam))./(2*alp);
            a2 = (-bet - sqrt(bet.^2 - 4*alp.*gam))./(2*alp);
            Cu1 = (ns(1:d)*a1' + vs(1:d,:));
            Cu2 = (ns(1:d)*a2' + vs(1:d,:));
            
            pInds = find((sum(abs(imag(Cu1))) <= 1e-12) & (sum(abs(imag(Cu2)))<=1e-12));
            Cu1(abs(imag(Cu1)) > 1e-12) = nan;
            Cu2(abs(imag(Cu2)) > 1e-12) = nan;
            
            % Distances between solutions
            D1 = sum((Cu - Cu1).^2);
            D2 = sum((Cu - Cu2).^2);
            D = (D1 - D2 > 0) + 1;
            Cup = (D==1).*Cu1 + (D==2).*Cu2;
            Uu = Cup - Cu;
            hold on;
        elseif(nDOF==2)
            % Change Coordinates
            P1 = [W([1:d],:) v0(1:d);...
                 zeros(1,d) 1]^-1;
            P2 = [W([1:d]+d,:) v0([1:d]+d);...
                  zeros(1,d) 1]^-1;
            Q1 = P1'*Q*P1;
            Q2 = P2'*Q*P2;
            % Create Mesh
            [xx, yy] = meshgrid(linspace(R(1,1),R(1,2),nC(1)),...
                                linspace(R(2,1),R(2,2),nC(2)));
            % Evaluate Conic Along Mesh
            F1 = Q1(1,1)*xx.^2 + Q1(2,2)*yy.^2 + Q1(3,3) +...
                2*Q1(2,1)*xx.*yy +...
                2*Q1(3,1)*xx + 2*Q1(3,2)*yy;
            F2 = Q2(1,1)*xx.^2 + Q2(2,2)*yy.^2 + Q2(3,3) +...
                 2*Q2(2,1)*xx.*yy +...
                 2*Q2(3,1)*xx + 2*Q2(3,2)*yy;
            % Generate contour along conic solution at 0
            hold on;
            C = contour(xx,yy,F1,[0 0], 'linewidth', LW_SS, 'color', C_SS(j,:));
            contour(xx,yy,F2,[0 0], '--', 'linewidth', LW_SS, 'color', C_SS(j,:));
            Cu = C(:,floor(linspace(ceil(size(C,2)/nP),size(C,2)-1,nP)));
            % Generate displacements along curve
            S = [W v0] * P1 * [Cu; ones(1, size(Cu,2))];
            pInds = 1:size(Cu,2);
            Cup = S(3:4,:);
            Uu = Cup - Cu;
        end
        % Plot
        quiver(Cu(1,pInds), Cu(2,pInds), Uu(1,pInds)*vU, Uu(2,pInds)*vU, 0, 'linewidth', LW_UA, 'color', C_SS(j,:));
        quiver(X0(1,:), X0(2,:), Us(1,:,j)*vS, Us(2,:,j)*vS, 0, 'linewidth', LW_SA, 'color', C_SA(j,:));
        quiver(X0(1,:), X0(2,:), Us(1,:,j)*vS, Us(2,:,j)*vS, 0, 'linewidth', LW_SA/3, 'color', [1 1 1]);
        plot(Cu(1,pInds), Cu(2,pInds), 'o', 'markersize', MS_UN, 'linewidth', MS_UN, 'color', C_SS(j,:));
        plot(Cup(1,pInds), Cup(2,pInds), 'o', 'markersize', MS_UN, 'linewidth', MS_UN, 'color', [1 1 1]);
        plot(Cup(1,pInds), Cup(2,pInds), 'o', 'markersize', MS_UN+2, 'linewidth', MS_UN-2, 'color', C_SS(j,:));
        plot(X0(1,:), X0(2,:), 'o', 'linewidth', MS_SN, 'markersize', MS_SN, 'color', C_SN);
        plot(XF(1,:), XF(2,:), 'o', 'linewidth', MS_SN, 'markersize', MS_SN, 'color', [1 1 1]);
        plot(XF(1,:), XF(2,:), 'o', 'linewidth', MS_SN-6, 'markersize', MS_SN+6, 'color', C_SN);
        hold off;
        set(gca,'visible',0);
        set(gcf,'color','w');
    elseif(d==3)
        % Spherical point
        [xSp, ySp, zSp] = sphere(20);
        xSp = xSp/10; 
        ySp = ySp/10; 
        zSp = zSp/10; 
        if(nDOF==3)
            % Change Coordinates
            P1 = [W([1:d],:) v0(1:d);...
                 zeros(1,d) 1]^-1;
            P2 = [W([1:d]+d,:) v0([1:d]+d);...
                  zeros(1,d) 1]^-1;
            Q1 = P1'*Q*P1;
            Q2 = P2'*Q*P2;
            % Create Mesh
            [xx, yy, zz] = meshgrid(linspace(R(1,1),R(1,2),nC(1)),...
                                   linspace(R(2,1),R(2,2),nC(2)),...
                                   linspace(R(3,1),R(3,2),nC(3)));
            % Evaluate Conic Along Mesh
            F1 = Q1(1,1)*xx.^2 + Q1(2,2)*yy.^2 + Q1(3,3)*zz.^2 + Q1(4,4)*1 +...
                 2*Q1(2,1)*xx.*yy + 2*Q1(2,3)*yy.*zz + 2*Q1(1,3)*xx.*zz +...
                 2*Q1(4,1)*xx + 2*Q1(4,2)*yy + 2*Q1(4,3)*zz;
            F2 = Q2(1,1)*xx.^2 + Q2(2,2)*yy.^2 + Q2(3,3)*zz.^2 + Q2(4,4)*1 +...
                 2*Q2(2,1)*xx.*yy + 2*Q2(2,3)*yy.*zz + 2*Q2(1,3)*xx.*zz +...
                 2*Q2(4,1)*xx + 2*Q2(4,2)*yy + 2*Q2(4,3)*zz;
            % Generate contour along conic solution at 0
            hold on;
            fI = isosurface(xx, yy, zz, F1, 0);
            p1 = patch(fI);
            isonormals(xx,yy,zz,F1,p1)
            fI2 = isosurface(xx, yy, zz, F2, 0);
            p2 = patch(shrinkfaces(fI2,.6));
            isonormals(xx,yy,zz,F2,p2)
            C = fI.vertices';
            S = [W v0] * P1 * [C; ones(1, size(C,2))];
            Cup = S(4:6,:);
            vInds = find(prod((C < R(:,2)) & (C > R(:,1)) & (Cup < R(:,2)) & (Cup > R(:,1))));
            vSamp = floor(linspace(2,length(vInds)-1,nP));
            Cu = C(:,vInds(vSamp));
            Cup = Cup(:,vInds(vSamp));
            Uu = Cup - Cu;
            p1.FaceColor = C_SS(j,:);
            p1.FaceAlpha = .5;
            p1.EdgeColor = 'none';
            p2.FaceColor = C_SS(j,:);
            p2.FaceAlpha = .5;
            p2.EdgeColor = 'none';
        elseif(nDOF==2)
            % First Coordinate Change to Cartesian Coordinates
            P = [W([1:d],:) v0(1:d)]^-1;
            P2 = [W([1:d]+d,:) v0([1:d]+d)]^-1;
            % Define Plane of Intersection
            p = P(end,:)';
            p2 = P2(end,:)';
            % Second Coordinate Change to Plane z=1
            PS = [null(p') p/(p'*p)];
            P2S = [null(p2') p2/(p2'*p2)];
            % Create Mesh
            [xx, yy] = meshgrid(linspace(R(1,1),R(1,2),nC(1)),...
                                linspace(R(2,1),R(2,2),nC(2)));
            XA = PS(1:2,1:2)^-1 * [xx(:) yy(:)]';
            XA2 = P2S(1:2,1:2)^-1 * [xx(:) yy(:)]';
            xx1 = reshape(XA(1,:), size(xx));
            yy1 = reshape(XA(2,:), size(yy));
            xx2 = reshape(XA2(1,:), size(xx));
            yy2 = reshape(XA2(2,:), size(yy));
            % Evaluate Conic Along Mesh, with added constraint
            Q1 = PS'*P'*Q*P*PS;
            Q2 = P2S'*P2'*Q*P2*P2S;
            F1 = Q1(1,1)*xx1.^2 + Q1(2,2)*yy1.^2 + Q1(3,3)*1 +...
                 2*Q1(2,1)*xx1.*yy1 + 2*Q1(3,2)*yy1 + 2*Q1(3,1)*xx1;
            F2 = Q2(1,1)*xx2.^2 + Q2(2,2)*yy2.^2 + Q2(3,3)*1 +...
                 2*Q2(2,1)*xx2.*yy2 + 2*Q2(3,2)*yy2 + 2*Q2(3,1)*xx2;
            % Generate contour along conic solution at 1 in second coordinates
            C = contour(xx1,yy1,F1,[0 0], 'color', 'b', 'linestyle', 'none');
            C2 = contour(xx2,yy2,F2,[0 0], 'color', 'b', 'linestyle', 'none');
            C = [C(:,2:end); ones(1, size(C,2)-1)];
            C2 = [C2(:,2:end); ones(1, size(C2,2)-1)];
            % Convert second coordinates back to cartesian
            C = PS*C; C = C(:,sum(C > min(R,[],2) & C < max(R,[],2))==3);
            C2 = P2S*C2; C2 = C2(:,sum(C2 > min(R,[],2) & C2 < max(R,[],2))==3);
            % Map initial points to final points
            Cu = C(:,floor(linspace(2,size(C,2)-1,nP)));
            S = [W v0] * P * Cu;
            Cup = S([1:d]+d,:);
            Uu = Cup - Cu;
            % Plot
            hold on;
            plot3(C(1,:), C(2,:), C(3,:), 'b-', 'linewidth', LW);
            plot3(C2(1,:), C2(2,:), C2(3,:), 'b--', 'linewidth', LW);
            plot3(Cu(1,:), Cu(2,:), Cu(3,:), 'o', 'markersize', 4, 'linewidth', 4, 'color', 'b');
            plot3(Cup(1,:), Cup(2,:), Cup(3,:), 'o', 'markersize', 8, 'linewidth', 2, 'color', 'b');
            quiver3(Cu(1,:), Cu(2,:), Cu(3,:), Uu(1,:,j)*vU, Uu(2,:,j)*vU, Uu(3,:,j)*vU, 0, 'color', [0 0 cR(j)], 'linewidth', LW/2);
            quiver3(X0(1,:), X0(2,:), X0(3,:), Us(1,:,j)*vS, Us(2,:,j)*vS, Us(3,:,j)*vS, 0, 'filled', 'color', [cR(j) 0 0], 'linewidth', LW);
            plot3(X0(1,:), X0(2,:), X0(3,:), 'o', 'markersize', MS, 'linewidth', MW, 'color', 'r');
            plot3(X0(1,:)+Us(1,:), X0(2,:)+Us(2,:), X0(3,:)+Us(3,:), 'o', 'markersize', MSH, 'linewidth', MWH, 'color', 'r');
            hold off;
        end
        % Plot
        quiver3(Cu(1,:), Cu(2,:), Cu(3,:), Uu(1,:,j)*vU, Uu(2,:,j)*vU, Uu(3,:,j)*vU, 0,...
            'filled', 'color', C_SS(j,:), 'linewidth', LW_UA);
        quiver3(X0(1,:), X0(2,:), X0(3,:), Us(1,:,j)*vS, Us(2,:,j)*vS, Us(3,:,j)*vS, 0,...
            'filled', 'color', C_SA(j,:), 'linewidth', LW_SA);
        quiver3(X0(1,:), X0(2,:), X0(3,:), Us(1,:,j)*vS, Us(2,:,j)*vS, Us(3,:,j)*vS, 0,...
            'filled', 'color', [1 1 1], 'linewidth', LW_SA/3);
        for i = 1:size(Cu,2)
            s = surf(xSp*MS_UN/MS_SN+Cu(1,i), ySp*MS_UN/MS_SN+Cu(2,i), zSp*MS_UN/MS_SN+Cu(3,i));
            s.FaceColor = C_SS(j,:);
            s.EdgeColor = 'none';
        end
        for i = 1:size(X0,2)
            s = surf(xSp+X0(1,i), ySp+X0(2,i), zSp+X0(3,i));
            s.FaceColor = C_SN;
            s.EdgeColor = 'none';
        end
        hold off;
        set(gca,'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[],...
                'XTick',[],'YTick',[],'ZTick',[],'box','on','boxstyle','back');
    end
end
end