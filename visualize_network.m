function [] = visualize_network(Xs, Xu, conn)
% This function is to visualize a constructed network.
%
% Inputs
% Xs:       dxNs (dimensions x nodes) array of specified node positions
% Xu:       dxN (dimensions x nodes) array of unspecified node positions
% conn:     kx2 array of edges from specified node i to unspecified node j
%
% Outputs
% Figure. Assumes the desired figure placement is already selected


% Plot Parameters
C_SN = [255 100 100]/255;       % Color of Specified Node   
C_UN = [100 100 255]/255;       % Color of Unspecified Node
d = size(Xs,1);                 % Number of Dimensions
ms = 4;                         % Marker Size
lw = 1;                         % Line Width
ea = .5;                        % Edge Transparency

hold on
if(d==2)
    line([Xs(1,conn(:,1)); Xu(1,conn(:,2))],...
         [Xs(2,conn(:,1)); Xu(2,conn(:,2))],...
         'linewidth', lw, 'color', [0 0 0 ea]);
    plot(Xs(1,:), Xs(2,:), 'o', 'linewidth', ms, 'markersize', ms, 'color', C_SN)
    plot(Xu(1,:), Xu(2,:), 'o', 'linewidth', ms, 'markersize', ms, 'color', C_UN(1,:));
    set(gca,'visible',0);
    set(gcf,'color','w');
elseif(d==3)
    % Spherical point
    [xSp, ySp, zSp] = sphere(20);
    xSp = xSp/10; 
    ySp = ySp/10; 
    zSp = zSp/10; 
    line([Xs(1,conn(:,1)); Xu(1,conn(:,2))],...
         [Xs(2,conn(:,1)); Xu(2,conn(:,2))],...
         [Xs(3,conn(:,1)); Xu(3,conn(:,2))],...
         'linewidth', lw, 'color', [0 0 0 ea]);
    for i = 1:size(Xu,2)
        s = surf(xSp+Xu(1,i), ySp+Xu(2,i), zSp+Xu(3,i));
        s.FaceColor = C_UN(1,:);
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
hold off;

end