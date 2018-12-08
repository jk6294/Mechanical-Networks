function [Xu, fV] = construct_network_finite(X0, XF, Xui, conn, pV)
% Function for constructing networks given topology and initial guess
%
% Inputs
% X0: d x n matrix of specified node initial positions
% XF: d x n x z matrix of m specified node final positions
% Xui: d x k matrix of initial positions of unspecified nodes
% conn: s x 2 matrix of connections from sp. node i to unsp. node j
% pV: scalar: 1 to plot networks, 0 otherwise
%
% Outputs
% Xu: 2dxk matrix of unspecified node positions
% fV: 1xk vector of absolute value of conic at Xu (ideally 0)

% Visualization Parameters
C_SN = [255 100 100]/255;       % Color of Specified Node
C_UN = [100 100 255]/255;...    % Color of Unspecified Node
    
% Initial values
d = size(X0,1);
k = size(Xui,2);
z = size(XF, 3);

% Iterate across unspecified nodes and optimize position
Xu = zeros(2*d,k);
fV = zeros(1,k);

% Optimization Options
options = optimset('TolFun', 1e-32, 'TolX', 1e-32);

% Solve for positions
for i = 1:k
    sInds = conn(find(conn(:,2)==i),1);     % Specified nodes
    f = cell(1,z);
    for j = 1:z
        [Q, W, v, fV] = construct_conic_finite(X0(:,sInds), XF(:,sInds,j));
        if(z > 1)
            P = [W(1:d,:), v(1:d);...
                 zeros(1,d), 1]^-1;
            Q = P'*Q*P;
        end
        f{j} = @(x) ([x 1] * Q * [x 1]')^2;
    end
    fF = @(x) sum(cellfun(@(F)F(x), f));
    m = size(Q,1)-1;
    % Convert initial condition to coordinate representation
    if(z == 1)
        cP1 = double(W(1:d,:)\(Xui(:,i)-v(1:d)));
        x0 = cP1(1:m,:)';
    else
        x0 = Xui(:,i)';
    end
    [XP, fVal] = fminsearch(fF, x0, options);
    if(z > 1)
        Xu(1:d,i) = XP';
    else
        Xu(:,i) = [W(1:(2*d),:) v(1:(2*d))] * [XP 1]';
    end
    fV(i) = sqrt(fVal);
end

% Plot
if(pV == 1)
    % Plot Parameters
    ms = 4;         % Marker Size
    lw = 1;         % Line Width
    ea = .5;        % Edge Transparency

    hold on
    if(d==2)
        line([X0(1,conn(:,1)); Xu(1,conn(:,2))],...
             [X0(2,conn(:,1)); Xu(2,conn(:,2))],...
             'linewidth', lw, 'color', [0 0 0 ea]);
        plot(X0(1,:), X0(2,:), 'o', 'linewidth', ms, 'markersize', ms, 'color', C_SN)
        plot(Xu(1,:), Xu(2,:), 'o', 'linewidth', ms, 'markersize', ms, 'color', C_UN);
        set(gca,'visible',0);
        set(gcf,'color','w');
    elseif(d==3)
        % Spherical point
        [xSp, ySp, zSp] = sphere(20);
        xSp = xSp/10; 
        ySp = ySp/10; 
        zSp = zSp/10; 
        % Network
        line([X0(1,conn(:,1)); Xu(1,conn(:,2))],...
             [X0(2,conn(:,1)); Xu(2,conn(:,2))],...
             [X0(3,conn(:,1)); Xu(3,conn(:,2))],...
             'linewidth', lw, 'color', [0 0 0 ea]);
        for i = 1:size(X0,2)
            s = surf(xSp+X0(1,i), ySp+X0(2,i), zSp+X0(3,i));
            s.FaceColor = C_SN;
            s.EdgeColor = 'none';
        end
        for i = 1:size(Xu,2)
            s = surf(xSp+Xu(1,i), ySp+Xu(2,i), zSp+Xu(3,i));
            s.FaceColor = C_UN;
            s.EdgeColor = 'none';
        end
        set(gca,'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[],...
                'XTick',[],'YTick',[],'ZTick',[],'box','on','boxstyle','back');
    end
    hold off;
end
end