function [Xu fV] = construct_network(Xs, Us, Xui, conn, pV)
% Function for constructing networks given topology and initial guess
%
% Inputs
% Xs: d x n matrix of specified node positions
% Us: d x n x z matrix of m specified node motions
% Xui: d x k matrix of initial positions of unspecified nodes
% conn: s x 2 matrix of connections from sp. node i to unsp. node j
% pV: scalar: 1 to plot networks, 0 otherwise
%
% Outputs
% Xu: dxk matrix of unspecified node positions
% fV: 1xk vector of absolute value of conic at Xu (ideally 0)

% Initial values
d = size(Xs,1);
n = size(Xs,2);
k = size(Xui,2);
z = size(Us, 3);

% Iterate across unspecified nodes and optimize position
Xu = zeros(d,k);
fV = zeros(1,k);

% Optimization Options
options = optimset('TolFun', 1e-16);

% Solve for positions
for i = 1:k
    sInds = conn(find(conn(:,2)==i),1);     % Specified nodes
    f = cell(1,z);
    for j = 1:z
        [Q, W, v, fV] = construct_conic(Xs(:,sInds), Us(:,sInds,j));
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
        cP1 = double([W(1:d,:) v(1:d)]\Xui(:,i));
        x0 = cP1(1:m,:)';
    else
        x0 = Xui(:,i)';
    end
    [XP, fVal] = fminsearch(fF, x0, options);
    if(z > 1)
        Xu(:,i) = XP';
    else
        Xu(:,i) = [W(1:d,:) v(1:d)] * [XP 1]';
    end
    fV(i) = sqrt(fVal);
end

% Plot
if(pV == 1)
    % Plot Parameters
    ms = 10;        % Marker Size
    lw = 2;         % Line Width
    ea = .5;        % Edge Transparency

    hold on
    if(d==2)
        line([Xs(1,conn(:,1)); Xu(1,conn(:,2))],...
             [Xs(2,conn(:,1)); Xu(2,conn(:,2))],...
             'linewidth', lw, 'color', [0 0 0 ea]);
        plot(Xs(1,:), Xs(2,:), 'ro', 'linewidth', ms)
        plot(Xu(1,:), Xu(2,:), 'bo', 'linewidth', ms);
    elseif(d==3)
        line([Xs(1,conn(:,1)); Xu(1,conn(:,2))],...
             [Xs(2,conn(:,1)); Xu(2,conn(:,2))],...
             [Xs(3,conn(:,1)); Xu(3,conn(:,2))],...
             'linewidth', lw, 'color', [0 0 0 ea]);
        plot3(Xs(1,:), Xs(2,:), Xs(3,:), 'ro', 'linewidth', ms)
        plot3(Xu(1,:), Xu(2,:), Xu(3,:), 'bo', 'linewidth', ms);
        view(20, 10);
    end
    hold off;
end
end