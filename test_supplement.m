%% Prepare Space
clear; clc;


%% Figure 1
az = -175; el = 10;         % Views
Xs = [[-1 0 1 0];...
      [-1 1 -1 0];...
	  [0 0 0 1]]*.68;
Us1 = [[-1.0  0.0  1.0  0.0];...
       [-1.0  1.0 -1.0  0.0];...
       [0.0  0.0  0.0  1.0]]*.4;
Us2 = [[0.7 -0.4  1.2 -0.2];...
       [1.0 -0.2 -0.5  0.9];...
       [-1.0  0.0  0.0  0.0]]*.4;
Us3 = [[0.7 0.4  1.2 -0.2];...
       [1.0 -0.2  0.0  0.9];...
       [-0.0  0.5  0.2  0.0]]*.4;

fig = figure(1); clf;
subplot(6,3,[1 4 7]);
visualize_conic(Xs, cat(3,Us1,Us2), [-1 1; -1 1; -1 1]*1.8, [100 100 100], 0, 1, 0);
axis(1.2*[-1 1 -1 1 -1 1]);
view(az, el); camlight(50, 30); lighting gouraud; material([.4 1 0]);

subplot(6,3,[1 4 7]+1);
visualize_conic(Xs, cat(3,Us1,-.4*Us1+Us2), [-1 1; -1 1; -1 1]*1.8, [100 100 100], 0, 1, 0);
axis(1.2*[-1 1 -1 1 -1 1]);
view(az, el); camlight(50, 30); lighting gouraud; material([.4 1 0]);

subplot(6,3,[1 4 7]+2);
visualize_conic(Xs, cat(3,Us1,-.8*Us1+Us2), [-1 1; -1 1; -1 1]*1.8, [100 100 100], 0, 1, 0);
axis(1.2*[-1 1 -1 1 -1 1]);
view(az, el); camlight(50, 30); lighting gouraud; material([.4 1 0]);

subplot(6,3,[1 4 7]+9);
visualize_conic(Xs, cat(3,Us1,Us3), [-1 1; -1 1; -1 1]*1.8, [100 100 100], 0, 1, 0);
axis(1.2*[-1 1 -1 1 -1 1]);
view(az, el); camlight(50, 30); lighting gouraud; material([.4 1 0]);

subplot(6,3,[1 4 7]+10);
visualize_conic(Xs, cat(3,Us1,+.2*Us1+Us3), [-1 1; -1 1; -1 1]*1.8, [100 100 100], 0, 1, 0);
axis(1.2*[-1 1 -1 1 -1 1]);
view(az, el); camlight(50, 30); lighting gouraud; material([.4 1 0]);

subplot(6,3,[1 4 7]+11);
visualize_conic(Xs, cat(3,Us1,+.4*Us1+Us3), [-1 1; -1 1; -1 1]*1.8, [100 100 100], 0, 1, 0);
axis(1.2*[-1 1 -1 1 -1 1]);
view(az, el); camlight(50, 30); lighting gouraud; material([.4 1 0]);


% Size and Save Figure
fName = 'space_set';
set(gcf, 'Renderer', 'opengl'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 8.2 4.5];
fig.PaperSize = [8.2 4.5];
% print(['Figures\' fName], '-dpng','-r1500');


%% Figure 2
% Positions
XsSSS = [[-1 -1 1 1];...
         [-1 1 1 -1]];
% Motions
UsSSS = zeros(2,4,2);
UsSSS(:,:,1) = [[-1 -1 1 1];...
               [-1 1 1 -1]]/2;
UsSSS(:,:,2) = [[.2 1 -1 -.2];...
               [-1 .5 .5 -1]]/2;
% Module Connectivity
conn = [1 1; 2 1; 3 1];
% Initial Position Conditions
Xu0 = [-1 -1 1 1;...
       -1 0 0 -1];
   
% Visualize Intersections and Place Nodes
fig = figure(2); clf;
subplot(4,3,[1,4]);
visualize_conic(XsSSS(:,[2 3 4]), UsSSS(:,[2 3 4],:), [-1 1; -1 1]*2, [100 100], 0, 1, 1);
[XuSSS1, fV] = construct_network(XsSSS(:,[2 3 4]), UsSSS(:,[2 3 4],:), Xu0(:,1), conn, 1);
axis(2*[-1 1 -1 1]);
subplot(4,3,[1,4]+1);
visualize_conic(XsSSS(:,[1 3 4]), UsSSS(:,[1 3 4],:), [-1 1; -1 1]*2, [100 100], 0, 1, 1);
[XuSSS2, fV] = construct_network(XsSSS(:,[1 3 4]), UsSSS(:,[1 3 4],:), Xu0(:,2), conn, 1);
axis(2*[-1 1 -1 1]);
subplot(4,3,[1,4]+6);
visualize_conic(XsSSS(:,[1 2 4]), UsSSS(:,[1 2 4],:), [-1 1; -1 1]*2, [100 100], 0, 1, 1);
[XuSSS3, fV] = construct_network(XsSSS(:,[1 2 4]), UsSSS(:,[1 2 4],:), Xu0(:,3), conn, 1);
axis(2*[-1 1 -1 1]);
subplot(4,3,[1,4]+7);
visualize_conic(XsSSS(:,[1 2 3]), UsSSS(:,[1 2 3],:), [-1 1; -1 1]*2, [100 100], 0, 1, 1);
[XuSSS4, fV] = construct_network(XsSSS(:,[1 2 3]), UsSSS(:,[1 2 3],:), Xu0(:,4), conn, 1);
axis(2*[-1 1 -1 1]);
XuSSS = [XuSSS1 XuSSS2 XuSSS3 XuSSS4];
connSSS = [2 1; 3 1; 4 1; 1 2; 3 2; 4 2; 1 3; 2 3; 4 3; 1 4; 2 4; 3 4];
subplot(4,3,[1,4]+5)
visualize_network(XsSSS, XuSSS, connSSS);
axis(2*[-1 1 -1 1]);
[S, Xd] = rigidity(XsSSS, XuSSS, connSSS);


% Size and Save Figure
fName = 'pre_stress';
set(gcf, 'Renderer', 'opengl'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 8.2 4.1];
fig.PaperSize = [8.2 4.1];
% print(['Figures\' fName], '-dpng','-r1500');


%% Figure 3
% Positions
XsSSS = [[-1 -1 1 1];...
         [-1 1 1 -1]];
% Motions
UsSSS2 = zeros(2,4,2);
UsSSS2(:,:,1) = [[-1 -1 1 1];...
                 [-1 1 1 -1]]/2;
UsSSS2(:,:,2) = [[.2 -1 1 -.2];...
                 [-1 .5 .5 -1]]/2;
% Module Connectivity
conn = [1 1; 2 1; 3 1];
% Initial Position Conditions
Xu0 = [sqrt(2) 0 0 -sqrt(2);...
       0       1 1  0];
             
% Visualize Intersections and Place Nodes
fig = figure(3); clf;
subplot(4,3,[1,4]);
visualize_conic(XsSSS(:,[2 3 4]), UsSSS2(:,[2 3 4],:), [-1 1; -1 1]*2, [100 100], 0, 1, 1);
[XuSSS1, fV] = construct_network(XsSSS(:,[2 3 4]), UsSSS2(:,[2 3 4],:), Xu0(:,1), conn, 1);
axis(2*[-1 1 -1 1]);
subplot(4,3,[1,4]+1);
visualize_conic(XsSSS(:,[1 3 4]), UsSSS2(:,[1 3 4],:), [-1 1; -1 1]*2, [100 100], 0, 1, 1);
[XuSSS2, fV] = construct_network(XsSSS(:,[1 3 4]), UsSSS2(:,[1 3 4],:), Xu0(:,2), conn, 1);
axis(2*[-1 1 -1 1]);
subplot(4,3,[1,4]+6);
visualize_conic(XsSSS(:,[1 2 4]), UsSSS2(:,[1 2 4],:), [-1 1; -1 1]*2, [100 100], 0, 1, 1);
[XuSSS3, fV] = construct_network(XsSSS(:,[1 2 4]), UsSSS2(:,[1 2 4],:), Xu0(:,3), conn, 1);
axis(2*[-1 1 -1 1]);
subplot(4,3,[1,4]+7);
visualize_conic(XsSSS(:,[1 2 3]), UsSSS2(:,[1 2 3],:), [-1 1; -1 1]*2, [100 100], 0, 1, 1);
[XuSSS4, fV] = construct_network(XsSSS(:,[1 2 3]), UsSSS2(:,[1 2 3],:), Xu0(:,4), conn, 1);
axis(2*[-1 1 -1 1]);
XuSSS = [XuSSS1 XuSSS2 XuSSS3 XuSSS4];
connSSS = [2 1; 3 1; 4 1; 1 2; 3 2; 4 2; 1 3; 2 3; 4 3; 1 4; 2 4; 3 4];
subplot(4,3,[1,4]+5)
visualize_network(XsSSS, XuSSS, connSSS);
axis(2*[-1 1 -1 1]);

[S, Xd] = rigidity(XsSSS, XuSSS, connSSS);
[V, U] = eig(S);
xdotp = reshape(Xd*V*[sqrt(-U(2,2)/U(1,1)); 1],8,2);
xdotn = reshape(Xd*V*[-sqrt(-U(2,2)/U(1,1)); 1],8,2);

XsSSSp = XsSSS + 0.00001*xdotp(1:4,:)';
XuSSSp = XuSSS + 0.00001*xdotp(5:8,:)';
XsSSSn = XsSSS + 0.00001*xdotn(1:4,:)';
XuSSSn = XuSSS + 0.00001*xdotn(5:8,:)';

[XMp, fCp] = sim_motion(XsSSSp, XuSSSp, connSSS, .01, 200, xdotp,1);
[XMp, fCp] = sim_motion(XsSSSp, XuSSSp, connSSS, .01, 200, -xdotp,1);
[XMp, fCp] = sim_motion(XsSSSn, XuSSSn, connSSS, .01, 100, xdotn,1);
[XMp, fCp] = sim_motion(XsSSSn, XuSSSn, connSSS, .01, 100, -xdotn,1);

hold on;
quiver([XsSSS(1,:) XuSSS(1,:)], [XsSSS(2,:) XuSSS(2,:)],...
        -xdotp(:,1)', -xdotp(:,2)', 0, 'color', [76 187 23]/255);
quiver([XsSSS(1,:) XuSSS(1,:)], [XsSSS(2,:) XuSSS(2,:)],...
        -xdotn(:,1)', -xdotn(:,2)', 0, 'color', [50 255 50]/255);
hold off;

% Size and Save Figure
fName = 'branching';
set(gcf, 'Renderer', 'opengl'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 8.2 4.1];
fig.PaperSize = [8.2 4.1];
% print(['Figures\' fName], '-dpng','-r1500');


%% Figure 4
% Initial and Final Positions for 1 Module
Xs0 = 0.5*[-1 -1 1 1;...
           -1 1 1 -1];
XsT = Xs0/0.5;
% Connectivity for 1 Module
conn = [1 1; 2 1; 3 1; 4 1; 1 2; 2 2; 3 2; 4 2];
% Guess for Unspecified Node Position
x0 = [-sqrt(2) sqrt(2);...
       0       0];

fig = figure(4); clf;
% a: Finite Solution Space
subplot(4,5,[1 6])
visualize_conic_finite(Xs0, XsT, [-1 1; -1 1]*3, [100; 100],4, .75, .85);
axis(1.6*[-1 1 -1 1]);

% b: Finite Network and Motion
subplot(4,5,[1 6] + 10); cla;
[Xu, fV] = construct_network_finite(Xs0, XsT, x0, conn, 0);
[XMot, fC] = sim_motion(Xs0, Xu(1:2,:), conn, .01, 200, [Xs0 Xu(1:2,:)],1);
hold on;
plot(XsT(1,:), XsT(2,:), 'wo', 'markersize', 4, 'linewidth', 4);
plot(Xu(3,:), Xu(4,:), 'wo', 'markersize', 4, 'linewidth', 4);
plot(XsT(1,:), XsT(2,:), 'o', 'markersize', 7, 'linewidth', 1, 'color', [255 100 100]/255);
plot(Xu(3,:), Xu(4,:), 'o', 'markersize', 7, 'linewidth', 1, 'color', [100 100 255]/255);
hold off;
axis(1.6*[-1 1 -1 1]);

% c: Tesselate Network
subplot(4,5,[2 3 7 8 12 13 17 18]); cla;
[XsT, XuT, connT] = tesselate_network(Xs0, Xu(1:2,:), conn, [1 1]', [4 4]');
visualize_network(XsT, XuT, connT);
[XMot1, fC] = sim_motion(XsT, XuT, connT, .1, 170, [XsT XuT], 0);
axis(4*1.4*[-.5 1 -.5 1]);

% d: Expanded Network
subplot(4,5,[2 3 7 8 12 13 17 18]+2); cla;
visualize_network(XMot1(:,1:size(XsT,2),end), XMot1(:,[1:size(XuT,2)]+size(XsT,2),end), connT);
axis(4*1.4*[-.5 1 -.5 1]);

% Size and Save Figure
fName = 'tesselate_finite';
set(gcf, 'Renderer', 'opengl'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 8.2 2.6]*1;
fig.PaperSize = [8.2 2.6]*1;
% print(['Figures\' fName], '-dpng','-r1500');


%% Figure 5
Xs0 = [-sqrt(3)/2 0 sqrt(3)/2;...
       -1/2       1 -1/2];
XsT1 = 1.5*[-sqrt(3)/2 0 sqrt(3)/2;...
            -1/2       1 -1/2];
XsT2 = [-1/2 0 1/2;...
        -1/2 1.5 -1/2];
X0 = [-0.6  0.6  0.4;...
      -1.0 -1.0  1.6];
conn = [1 1; 1 2; 1 3; 2 1; 2 2; 2 3; 3 1; 3 2; 3 3];

fig = figure(5); clf;
subplot(1,8,[1 2]);
visualize_conic_finite(Xs0, XsT1, [-2 2; -2 2], [100;100], 0,1,1);
axis(1.5*[-1 1 -1 1]);

subplot(1,8,[3 4]);
visualize_conic_finite(Xs0, XsT2, [-2 2; -2 2], [100;100], 0,1,1);
axis(1.5*[-1 1 -1 1]);

subplot(1,8,[5 6]);
% Construct Network
[Xu, fV] = construct_network_finite(Xs0, cat(3,XsT1,XsT2), X0, conn, 1);
Xu = Xu(1:2,:);
LVal = sqrt((Xs0(1,conn(:,1))-Xu(1,conn(:,2))).^2 +...
            (Xs0(2,conn(:,1))-Xu(2,conn(:,2))).^2);
[M, E] = sim_motion3D_congrad(Xs0, Xu, conn, LVal, 0.001, 1001, [1 3], XsT1(:,[1 3]), 1);
axis(1.5*[-1 1 -1 1]);

subplot(1,8,[7 8]);
construct_network_finite(Xs0, cat(3,XsT1,XsT2), X0, conn, 1);
Xu = Xu(1:2,:);
[M2, E2] = sim_motion3D_congrad(Xs0, Xu, conn, LVal, 0.001, 1001, [1 3], XsT2(:,[1 3]), 1);
axis(1.5*[-1 1 -1 1]);

% Size and Save Figure
fName = 'tristable';
set(gcf, 'Renderer', 'opengl'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 8.2 1.8]*1;
fig.PaperSize = [8.2 1.8]*1;
% print(['Figures\' fName], '-dpng','-r1500');


%% Figure 6
az = -170; el = 7;         % Views
C_SS = [100 100 255;...         % Color of Solution Space
        100 200 255]/255;       
% Positions and Motions
d = 3;
Xs = [[-1 0 1 0];...
      [-1 1 -1 0];...
	  [0 0 0 1]]*.68;
Us = [[-1.0  0.0 -0.6  0.0];...
      [-1.0  0.5  1.0  0.0];...
      [0.0  0.0  0.0  -1]]*.4;
% Initial Positions
x0 = [-.7 -.5 .8]'*.6;
% Connectivity
conn = [1 1; 2 1; 3 1; 4 1];

fig = figure(6); clf;
subplot(1,3,1);
visualize_conic(Xs, Us, [-1 1; -1 1; -1 1]*1.8, [100 100 100], 0, 1, 0);
[Xu1, fV] = construct_network(Xs, Us, x0, conn, 1);
axis(1.2*[-1 1 -1 1 -1 1]);
view(az, el);
camlight(50,30); lighting gouraud; material([.4 1 0]);

subplot(1,3,2);
% Get Conic Parameters
[Q, W, v0, err] = construct_conic(Xs, Us);
% Change basis
P = [W(1:d,:), v0(1:d);...
     zeros(1,d), 1]^-1;
Q = P'*Q*P;
pInd = 1;
% Get Q*
A = Q(1:3,1:3);
Qs = [A -A*Xu1(:,pInd); -Xu1(:,pInd)'*A Xu1(:,pInd)'*A*Xu1(:,pInd)];
% Create Mesh
[xx, yy, zz] = meshgrid(linspace(-1.8,1.8,100),...
                        linspace(-1.8,1.8,100),...
                        linspace(-1.8,1.8,100));
F = Q(1,1)*xx.^2 + Q(2,2)*yy.^2 + Q(3,3)*zz.^2 + Q(4,4)*1 +...
    2*Q(2,1)*xx.*yy + 2*Q(2,3)*yy.*zz + 2*Q(1,3)*xx.*zz +...
    2*Q(4,1)*xx + 2*Q(4,2)*yy + 2*Q(4,3)*zz;
Fs = Qs(1,1)*xx.^2 + Qs(2,2)*yy.^2 + Qs(3,3)*zz.^2 + Qs(4,4)*1 +...
     2*Qs(2,1)*xx.*yy + 2*Qs(2,3)*yy.*zz + 2*Qs(1,3)*xx.*zz +...
     2*Qs(4,1)*xx + 2*Qs(4,2)*yy + 2*Qs(4,3)*zz;
hold on;
fI = isosurface(xx, yy, zz, F, 0);
p = patch(fI);
isonormals(xx,yy,zz,F,p)
            p.FaceColor = C_SS(1,:);
            p.FaceAlpha = 0.5;
            p.EdgeColor = 'none';
fIs = isosurface(xx, yy, zz, Fs, 0);
ps = patch(fIs);
isonormals(xx,yy,zz,Fs,ps)
            ps.FaceColor = C_SS(2,:);
            ps.FaceAlpha = 0.5;
            ps.EdgeColor = 'none';
hold off;
construct_network(Xs, Us, x0(:,pInd), conn(1:4,:), 1);
axis(1.2*[-1 1 -1 1 -1 1]);
view(az, el);
camlight(50,30); lighting gouraud; material([.4 1 0]);

% Construct Network
f = cell(1,2);
F = @(x) ([x 1] * Q * [x 1]')^2 + ([x 1] * Qs * [x 1]')^2;
options = optimset('TolFun', 1e-32, 'TolX', 1e-32);
[Xuu1, fVal] = fminsearch(F, [.5 0 -.2], options);
[Xuu2, fVal] = fminsearch(F, [-.8 -.5 .3], options);
Xu = [Xu1 Xuu1' Xuu2'];
conn = [conn; 1 2; 2 2; 3 2; 4 2; 1 3; 2 3; 3 3; 4 3];
subplot(1,3,3);
construct_motion(Xs, Us, Xu, conn, 1, .01, [], [pInd 2; pInd 3]);
axis(1.2*[-1 1 -1 1 -1 1]);
view(az, el);
camlight(50,30); lighting gouraud; material([.4 1 0]);

% Size and Save Figure
fName = 'nonbipartite';
set(gcf, 'Renderer', 'opengl'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 8.2 2.1]*1;
fig.PaperSize = [8.2 2.1]*1;
% print(['Figures\' fName], '-dpng','-r1500');



