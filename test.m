%% Prepare Space
clear; clc;


%% Figure 2
fig = figure(2); clf;
% Subplot Panel Sizes
nPan1 = [];
for i = 1:5; nPan1 = [nPan1, [1:7]+(i-1)*24+.5]; end
nPan2 = [];
for i = 1:6; nPan2 = [nPan2, [1:8]+(i-1)*24]; end
az = -160; el = 20;         % Views
azC = 50; elC = 30;         % Camera Position


% d: Positions and Motions
subplot(12,24,nPan1)
Xs = [[-1 1 0];...
      [0 0 1]];
Us = [[-1 1 -1];...
      [-1 -1 0]]/4;
visualize_conic(Xs, Us, [-1 1; -1 1]*2, [100 100], 14, 1, 1);
axis(1.3*[-1 1 -1 1]);


% e: Positions and Motions
subplot(12,24,nPan2+24*6)
Xs = [[-1 1 0 0 0];...
	  [-1 -1 1 0 0];...
	  [0 0 0 1 -1]];
Us = [[1 1 -1 0 0 ];...
	  [1 -1 0.1 0 0 ];...
	  [-1 -1 -1 1 1]]/4;
visualize_conic(Xs, Us, [-1 1; -1 1; -1 1]*4, [100 100 100], 40, 1, 1);
axis(1.5*[-1 1 -1 1 -1 1]);
view(az, el);
camlight(azC, elC); lighting gouraud; material([.4 1 0]);


% f: Positions and Motions
subplot(12,24,nPan1+8)
Xs = [[-1 1];...
      [0 0]];
Us = -[[-1 1];...
      [0 0]]/4;
visualize_conic(Xs, Us, [-1 1; -1 1]*1.7, [100 100], 7, 1, 1);
axis(1.8*[-1 1 -1 1]);


% g: Positions and Motions
subplot(12,24,nPan2+24*6+8)
Xs = [[-1 1 0 0];...
      [-1 -1 1 0];...
      [0 0 0 1]];
Us = [[-1 1 0 0]/2;...
      [-1 -1 1 0]/4;...
      [-1 -1 -1 -1.6]]/3;
visualize_conic(Xs, Us, [-1 1; -1 1; -1 1]*1.5, [100 100 100], 40, 1, 1);
axis(1.5*[-1 1 -1 1 -1 1]);
view(az, el);
camlight(azC, elC); lighting gouraud; material([.4 1 0]);


% h: Positions and Motions
subplot(12,24,nPan1+16)
Xs = [[-1 -1 1 1];...
      [-1 1 1 -1]];
Us = [[-1 -1 1 1];...
      [-1 1 1 -1]]/4;
visualize_conic(Xs, Us, [-1 1; -1 1]*2, [100 100], 16, 1, 1);
axis(1.5*[-1 1 -1 1]);


% i: Positions and Motions
subplot(12,24,nPan2+24*6+16)
Xs = [[-1 -1 -1 -1 1 1 1 1];...
      [-1 -1 1 1 -1 -1 1 1];...
	  [-1 1 -1 1 -1 1 -1 1]]*.8;
% Motions: 1D Manifold
Us = [[-1 -1 -1 -1 1 1 1 1];...
	  [-1 -1 1 1 -1 -1 1 1];...
	  [-1 1 -1 1 -1 1 -1 1]]/4;
visualize_conic(Xs, Us, [-1 1; -1 1; -1 1]*1.8, [100 100 100], 50, 1, 1);
axis(1.5*[-1 1 -1 1 -1 1]);
view(az, el);
camlight(azC, elC); lighting gouraud; material([.4 1 0]);


% Size and Save Figure
fName = 'bipartite';
set(gcf, 'Renderer', 'opengl'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 17.2 10.8];
fig.PaperSize = [17.2 10.8];
% print(['Figures\' fName], '-dpng','-r300');


%% Figure 3
fig = figure(3); clf;
% Subplot Panel Sizes
nPan1 = [];
for i = 1:5; nPan1 = [nPan1, [1:6]+(i-1)*24+1]; end
nPan2 = [];
for i = 1:6; nPan2 = [nPan2, [1:8]+(i-1)*24]; end
az = -160; el = 10;         % Views


% Positions and Motions
Xs = [[-1 -1 1 1];...
	  [-1 1 1 -1]];
Us = [[-1 -1 1 1];...
	  [-1 1 1 -1]]/2;
% Two different connectivities
conn1 = [1 1; 2 1; 3 1; 4 1;...
         1 2; 2 2; 3 2; 4 2];
conn2 = [1 2; 1 3; 2 1; 2 2; 2 3; 2 4;...
         3 1; 3 4; 4 1; 4 2; 4 3; 4 4];
% Initial guesses for unspecified node positions
xS1 = [[-sqrt(2) sqrt(2)];...
       [0 0]];
xS2 = [[-sqrt(2) 0 sqrt(2) 0];...
       [0 sqrt(2) 0 -sqrt(2)]];

   
% a: Visualize Conics
subplot(12,24,nPan1)
visualize_conic(Xs, Us, [-1 1; -1 1]*2, [100 100], 0, 1, 0);
axis(1.5*[-1 1 -1 1]);


% b: Find positions of unspecified nodes, and corresponding motions
subplot(12,24,nPan1+8)
[Xu1 fV1] = construct_network(Xs, Us, xS1, conn1, 0);
[Us1, Uu1, err1] = construct_motion(Xs, Us, Xu1, conn1, 1, 1);
axis(1.5*[-1 1 -1 1]);


% c:  Find positions of unspecified nodes, and corresponding motions
subplot(12,24,nPan1+16)
[Xu2 fV2] = construct_network(Xs, Us, xS2, conn2, 0);
[Us2, Uu2, err2] = construct_motion(Xs, Us, Xu2, conn2, 1, 1);
axis(1.5*[-1 1 -1 1]);


% Positions and Motions
Xs = [[-1 -1 -1 -1 1 1 1 1];...
	  [-1 -1 1 1 -1 -1 1 1];...
	  [-1 1 -1 1 -1 1 -1 1]]*.8;
Us = [[-1 -1 -1 -1 1 1 1 1];...
	  [-1 -1 1 1 -1 -1 1 1];...
	  [-1 1 -1 1 -1 1 -1 1]]/2;
% Two different connectivities
conn1 = [1 1; 2 1; 3 1; 4 1; 5 1; 6 1;...
         1 2; 3 2; 5 2; 6 2; 7 2; 8 2;...
         1 3; 2 3; 5 3; 6 3; 3 3; 8 3;...
         3 4; 4 4; 5 4; 2 4; 7 4; 8 4;...
         1 5; 3 5; 5 5; 4 5; 7 5; 8 5;...
         2 6; 4 6; 5 6; 6 6; 8 6];
conn2 = [1 1; 2 1; 3 1; 4 1; 5 1; 6 1; 7 1; 8 1;...
         1 2; 2 2; 3 2; 4 2; 5 2; 6 2; 7 2; 8 2;...
         1 3; 2 3; 3 3; 5 3; 6 3;...
         3 4; 4 4; 5 4; 7 4; 8 4;...
         1 5; 2 5; 3 5; 5 5; 7 5; 8 5];
% Initial guesses for unspecified node positions
xS1 = [[-sqrt(3)  sqrt(3)  0        0        0        0      ];...
       [0         0       -sqrt(3)  sqrt(3)  0        0      ];...
	   [0         0        0        0       -sqrt(3)  sqrt(3)]]*.8;
xS2 = [[-sqrt(3)  sqrt(3)  0        0        0      ];...
	   [0         0       -sqrt(3)  sqrt(3)  0      ];...
	   [0         0        0        0       -sqrt(3)]]*.8;
   

% d: Visualize conic solution space
subplot(12,24,nPan2+24*6)
visualize_conic(Xs, Us, [-1 1; -1 1; -1 1]*1.8, [100 100 100], 0, 1, 1);
axis(1.5*[-1 1 -1 1 -1 1]);
view(az, el);
camlight(azC, elC); lighting gouraud; material([.4 1 0]);


% e: Find positions of unspecified nodes, and corresponding motions
subplot(12,24,nPan2+24*6+8)
[Xu1 fV1] = construct_network(Xs, Us, xS1, conn1, 0);
[Us1, Uu1, err1] = construct_motion(Xs, Us, Xu1, conn1, 1, 1);
axis(1.5*[-1 1 -1 1 -1 1]);
view(az, el);
camlight(azC, elC); lighting gouraud; material([.4 1 0]);


% f: Find positions of unspecified nodes, and corresponding motions
subplot(12,24,nPan2+24*6+16)
[Xu2 fV2] = construct_network(Xs, Us, xS2, conn2, 0);
[Us2, Uu2, err2] = construct_motion(Xs, Us, Xu2, conn2, 1, 1);
axis(1.5*[-1 1 -1 1 -1 1]);
view(az, el);
camlight(azC, elC); lighting gouraud; material([.4 1 0]);


% Size and Save Figure
fName = 'construction';
set(gcf, 'Renderer', 'opengl'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 19.2 10.8];
fig.PaperSize = [19.2 10.8];
% print(['Figures\' fName], '-dpng','-r300');


%% Figure 4
fig = figure(4); clf;
% Subplot Panel Sizes
nPan1 = [];
for i = 1:5; nPan1 = [nPan1, [1:7]+(i-1)*24+.5]; end
nPan2 = [];
for i = 1:6; nPan2 = [nPan2, [1:7]+(i-1)*24+.5]; end
az = -175; el = 10;         % Views


% Positions and Motions
Xs = [[-1 0 1];...
	  [0 1 0]];
XsSSS = [[-1 -1 1 1];...
         [-1 1 1 -1]];
% Two separate motions
Us = zeros(2,3,2);
Us(:,:,1) = [[-1 0 1];...
             [0 1 0]]/2;
Us(:,:,2) = [[-1 1 1];...
             [.5 -1 .5]]/2;
% Two separate motions for singularity in rigidity matrix
UsSSS = zeros(2,4,2);
UsSSS(:,:,1) = [[-1 -1 1 1];...
               [-1 1 1 -1]]/2;
UsSSS(:,:,2) = [[.2 1 -1 -.2];...
               [-1 .5 .5 -1]]/2;
% Initial guess for unspecified nodes
x0 = [-.7 .7]';
% Connectivity
conn = [1 1; 2 1; 3 1];
connSSS = [1 1; 2 1; 3 1];


% a: Visualize solution space
subplot(12,24,nPan1)
visualize_conic(Xs, Us, [-1 1; -1 1]*2, [100 100], 0, 1, 0);
axis(1.5*[-1 1 -1 1]);


% c: Place unspecified nodes and construct motions
subplot(12,24,nPan1+8)
[Xu1, fV] = construct_network(Xs, Us, x0, conn, 0);
[Us1, Uu1, err] = construct_motion(Xs, Us, Xu1, conn, 1, 1);
axis(1.5*[-1 1 -1 1]);


% e: Construct states of self-stress
subplot(12,24,nPan1+16);
% subplot(2,2,1);
% visualize_conic(XsSSS(:,[2 3 4]), UsSSS(:,[2 3 4],:), [-1 1; -1 1]*2, [100 100], 0, 1, 1);
[XuSSS1, fV] = construct_network(XsSSS(:,[2 3 4]), UsSSS(:,[2 3 4],:), [-1 -1]', connSSS, 1);
% subplot(2,2,2);
% visualize_conic(XsSSS(:,[1 3 4]), UsSSS(:,[1 3 4],:), [-1 1; -1 1]*2, [100 100], 0, 1, 1);
[XuSSS2, fV] = construct_network(XsSSS(:,[1 3 4]), UsSSS(:,[1 3 4],:), [-1 0]', connSSS, 1);
% subplot(2,2,3);
% visualize_conic(XsSSS(:,[1 2 4]), UsSSS(:,[1 2 4],:), [-1 1; -1 1]*2, [100 100], 0, 1, 1);
[XuSSS3, fV] = construct_network(XsSSS(:,[1 2 4]), UsSSS(:,[1 2 4],:), [1 0]', connSSS, 1);
% subplot(2,2,4);
% visualize_conic(XsSSS(:,[1 2 3]), UsSSS(:,[1 2 3],:), [-1 1; -1 1]*2, [100 100], 0, 1, 1);
[XuSSS4, fV] = construct_network(XsSSS(:,[1 2 3]), UsSSS(:,[1 2 3],:), [1 -1]', connSSS, 1);
XuSSS = [XuSSS1 XuSSS2 XuSSS3 XuSSS4];
connSSS = [2 1; 3 1; 4 1; 1 2; 3 2; 4 2; 1 3; 2 3; 4 3; 1 4; 2 4; 3 4];
[Us1, Uu1, err] = construct_motion(XsSSS, UsSSS, XuSSS, connSSS, 1, 0);
axis(1.5*[-1 1 -1 1]);
% [XM, fC] = sim_motion(XsSSS, XuSSS, connSSS, .1, 110, -[XsSSS XuSSS],1);


% Positions and Motions
Xs = [[-1 0 1 0];...
      [-1 1 -1 0];...
	  [0 0 0 1]]*.8;
% Two separate motions
Us = zeros(3,4,2);
Us(:,:,1) = [[-1.0  0.0  1.0  0.0];...
             [-1.0  1.0 -1.0  0.0];...
             [0.0  0.0  0.0  1.0]]/2;
Us(:,:,2) = [[0.7 -0.4  1.2 -0.2];...
             [1.0 -0.2 -0.5  0.9];...
             [-1.0  0.0  0.0  0.0]]/2;
% Initial guess for unspecified nodes
x0 = [-.7 -.5 .8;...
      1.2 -.4 -.75;...
      -.4 -1.5 -.3;...
      .2 .6 .6]'*.8;
x0SSS = [-.7 -.5 .8;...
         1.2 -.4 -.75;...
         -.4 -1.5 -.3;...
         .2 .6 .6;...
         -.05 .58 -.85]'*.8;...
         
% Connectivity
conn = [1 1; 2 1; 3 1; 4 1;...
        1 2; 2 2; 3 2; 4 2;...
        1 3; 2 3; 3 3; 4 3;...
        1 4; 2 4; 3 4; 4 4];
connSSS = [1 1; 2 1; 3 1; 4 1;...
           1 2; 2 2; 3 2; 4 2;...
           1 3; 2 3; 3 3; 4 3;...
           1 4; 2 4; 3 4; 4 4;...
           1 5; 2 5; 3 5; 4 5];
    
% b: Visualize solution space
subplot(12,24,nPan2+24*6)
visualize_conic(Xs, Us, [-1 1; -1 1; -1 1]*1.8, [100 100 100], 0, 1, 0);
axis(1.5*[-1 1 -1 1 -1 1]);
view(az, el);
camlight(azC, elC); lighting gouraud; material([.4 1 0]);


% d: Place unspecified nodes and construct motion
subplot(12,24,nPan2+24*6+8)
[Xu1, fV] = construct_network(Xs, Us, x0, conn, 0);
[Us1, Uu1, err] = construct_motion(Xs, Us, Xu1, conn, 1, 1);
axis(1.5*[-1 1 -1 1 -1 1]);
view(az, el);
camlight(azC, elC); lighting gouraud; material([.4 1 0]);


% f: Place unspecified nodes and construct motion
subplot(12,24,nPan2+24*6+16)
[Xu1, fV] = construct_network(Xs, Us, x0SSS, connSSS, 0);
% Lengths
LVal = sqrt((Xs(1,connSSS(:,1)) - Xu1(1,connSSS(:,2))).^2 +...
            (Xs(2,connSSS(:,1)) - Xu1(2,connSSS(:,2))).^2 +...
            (Xs(3,connSSS(:,1)) - Xu1(3,connSSS(:,2))).^2);
[Us1, Uu1, err] = construct_motion(Xs, Us, Xu1, connSSS, .01, .01);
axis(1.5*[-1 1 -1 1 -1 1]);
view(az, el);
camlight(azC, elC); lighting gouraud; material([.4 1 0]);


% Size and Save Figure
fName = 'multi';
set(gcf, 'Renderer', 'opengl'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 19.2 10.8];
fig.PaperSize = [19.2 10.8];
% print(['Figures\' fName], '-dpng','-r300');


%% Figure 5b-d
fig = figure(5); clf;
% Subplot Panel Sizes
nPan1 = [];
for i = 1:6; nPan1 = [nPan1, [1:4]+(i-1)*24]; end
nPan2 = [];
for i = 1:5; nPan2 = [nPan2, [1:24]+(i-1)*24]; end


% Positions and Motions
Xs = [-1 -1 1 1;...
      -1 1 1 -1];
Us = [-1 -1 1 1;...
      -1 1 1 -1];
% Initial guess for unspecified nodes
x0 = [-sqrt(2) 0;...
       sqrt(2) 0]';
% Connectivity
conn = [1 1; 2 1; 3 1; 4 1; 1 2; 2 2; 3 2; 4 2];


% a: Independent modules
subplot(12,24,nPan1)
construct_network(Xs + [0;1.2], Us, x0 + [0;1.2], conn, 1);
construct_network(Xs - [0;1.2], Us, x0 - [0;1.2], conn, 1);
axis([-1.8 1.8 -2.8 2.8]);


% b: Combined modules with motion
subplot(12,24,nPan1 + 4)
[Xu, fV] = construct_network(Xs, Us, x0, conn, 0);
[XsT, XuT, connT] = tesselate_network(Xs, Xu, conn, [2 2]', [1 2]');
[Us, Uu, err] = construct_motion(XsT, zeros(0,0,0), XuT, connT, -2, -2);
axis([-1.8 1.8 -1.8 3.8]);


% c: Tesselation expanded
subplot(12,24,[nPan1 nPan1+4]+8)
[XsT, XuT, connT] = tesselate_network(Xs, Xu, conn, [2 2]', [5 5]');
[XMot1, fC] = sim_motion(XsT, XuT, connT, .1, 130, [XsT XuT],0);
construct_motion(XMot1(:,1:size(XsT,2),end), zeros(0,0,0), XMot1(:,[1:size(XuT,2)]+size(XsT,2),end), connT, .1, .1);
axis([-4 12 -2 10]);


% d: Tesselation contracted
subplot(12,24,[nPan1 nPan1+4]+16);
[XMot2, fC] = sim_motion(XsT, XuT, connT, .1, 130, -[XsT XuT],0);
construct_motion(XMot2(:,1:size(XsT,2),end), zeros(0,0,0), XMot2(:,[1:size(XuT,2)]+size(XsT,2),end), connT, .1, .1);
axis([-4 12 -2 10]);


% 3D
subplot(12,24,nPan2+24*6); cla;

% Positions
Xs1 = [[2 1 3 3.5 3.5];...
       [0 0 0 -1 1];...
       [1 2 3 2.5 2.5]];
Xs2 = [[6 7 5 4.5 4.5];...
       [0 0 0 -1 1];....
       [1 2 3 2.5 2.5]];
% Motions
Us1 = [[1 -1 1 0 0];...
       [0 0 0 1 -1];...
       [-1 1 0 0 0]]*.4;
Us2 = [[1 -1 1 0 0];...
       [0 0 0 1 -1];...
       [1 -1 0 0 0]]*.4;
% Initial guess for unspecified nodes
xS1 = [[3.092 0.808 2.183];...
       [2.758 -0.930 1.516];...
       [2.758 0.930 1.516];...
       [2 -.77 2.061];...
       [2 .77 2.061]]';
xS2 = [[5.081 -0.648 3.313];...
       [5.081 0.648 3.313];...
       [6.184 -1.000 1.737];...
       [5.121 -1.287 2.626];...
       [5.121 1.287 2.626]]';
% Connectivity
conn = [1 1; 2 1; 3 1; 4 1; 5 1;...
        1 2; 2 2; 3 2; 4 2; 5 2;...
        1 3; 2 3; 3 3; 4 3; 5 3;...
        1 4; 2 4; 3 4; 4 4;...
        1 5; 2 5; 3 5; 5 5];

% Simulate Motion
[Xu1 fV] = construct_network(Xs1, Us1, xS1, conn, 0);
[Usp, Uup, err] = construct_motion(Xs1, -Us1, Xu1, conn, 1, 1);

[Xu2 fV] = construct_network(Xs2, Us2, xS2, conn, 0);
[Usp, Uup, err] = construct_motion(Xs2, -Us2, Xu2, conn, 1, 1);

axis([.9 19.1 -3 1.4 .5 3.5]);
view(-8, 10);


[XsT, XuT, connT] = tesselate_network([Xs1 Xs2-[1;0;0]], [Xu1 Xu2-[1;0;0]], [conn; conn+[size(Xs1,2), size(Xu1,2)]], [0 0 0]', [1 1 1]');
[XMot, fC] = sim_motion3D(XsT+[7.2 0 0]', XuT+[7.2 0 0]', connT, .01, 70, ([XsT XuT]+[7.2 0 0]'), 0, [3 8]);
construct_motion(XMot(:,1:size(XsT,2),end), zeros(0,0,0), XMot(:,[1:size(XuT,2)]+size(XsT,2),end), connT, .01, .01, [3 8]);

[XsT, XuT, connT] = tesselate_network([Xs1 Xs2-[1;0;0]], [Xu1 Xu2-[1;0;0]], [conn; conn+[size(Xs1,2), size(Xu1,2)]], [0 0 0]', [1 1 1]');
[XMot, fC] = sim_motion3D(XsT+[13 0 0]', XuT+[13 0 0]', connT, .01, 110, -([XsT XuT]+[13 0 0]'), 0, [3 8]);
construct_motion(XMot(:,1:size(XsT,2),end), zeros(0,0,0), XMot(:,[1:size(XuT,2)]+size(XsT,2),end), connT, .01, .01, [3 8]);



% Effects
lightangle(0, 40);
lighting gouraud; 
material([.4 1 0]);


% Size and Save Figure
fName = 'combine';
set(gcf, 'Renderer', 'opengl'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 19.2 10.8];
fig.PaperSize = [19.2 10.8];
% print(['Figures\' fName], '-dpng','-r300');


%% Cooperativity Design
figure(6); clf;
nPan1 = [];
for i = 1:5; nPan1 = [nPan1, [1:4]+(i-1)*24]; end
nPan2 = [];
for i = 1:6; nPan2 = [nPan2, [1:24]+(i-1)*24]; end
az = 30; el = 25;         % Views
azC = 50; elC = 30;         % Camera Position
     
xS1 = [-2 -2 -1 -1;...
       -1  1  0  0;...
        0  0 -1  1]/2;
xF1 = xS1 - [ 0     0     0     0;...
             -0.2   0.2   0     0;...
              0     0    -0.06  0.06];
U = [-0.4  -0.4   0     0;...
      0.2  -0.2   0     0;...
      0     0     0.06 -0.06];
% subplot(2,3,1); cla;
visualize_conic_finite(xS1, xF1, [-2 -.5; -1 1; -1 1]*1, [100; 100; 100], 2, 1, 1);
visualize_conic(xS1, U, [-2 -.5; -1 1; -1 1]*1, [100; 100; 100], 0, 1, 1);
view(az, el);


%% Construct Test
figure(6); clf; 
x01 = [-0.879  0.175  0.657;...
       -0.803 -0.222  0.172;...
       -0.894 -0.100 -0.859;...
       -0.833  0.201 -0.455;...
       -0.561 -0.100 -0.273;...
       -0.542  0.111  0.172]'; 
% x01 = [-1.232  0.390  0.576;...
%        -1.232 -0.390  0.576;...
%        -1.183  0.333  0.172;...
%        -1.183 -0.333  0.172;...
%        -1.192  0.359 -0.434;...
%        -1.192 -0.359 -0.434]';
% x01 = [-1.354  0.440  0.677;...
%        -1.354 -0.440  0.677;...
%        -1.475  0.423  0.091;...
%        -1.475 -0.423  0.091;...
%        -1.394  0.455 -0.700;...
%        -1.394 -0.455 -0.700]';
% x01 = x01 + (rand(size(x01))-.5)*.1;
conn1 = [1 1; 2 1; 3 1; 4 1;...
         1 2; 2 2; 3 2; 4 2;...
         1 3; 2 3; 3 3; 4 3;...
         1 4; 2 4; 3 4; 4 4;...
         1 5; 2 5; 3 5; 4 5;...
         1 6; 2 6; 3 6; 4 6];
cla;
[Xu1, fV] = construct_network_finite(xS1, xF1, x01, conn1, 1);


%% Simulate Test
LVal1 = sqrt((xS1(1,conn1(:,1)) - Xu1(1,conn1(:,2))).^2 +...
             (xS1(2,conn1(:,1)) - Xu1(2,conn1(:,2))).^2 +...
             (xS1(3,conn1(:,1)) - Xu1(3,conn1(:,2))).^2);
[MSi, ESi] = sim_motion3D_congrad(xS1, Xu1(1:3,:), conn1, LVal1, 0.001, 1001, [1 2], xF1(:,[1 2]));


%% Video
figure(7); clf;
XMot = MSi;
EA = ESi;
connA = conn1;
XsA = xS1;
XuA = Xu1;

ns = size(XsA,2);
nu = size(XuA,2);

% kTC = num2cell(kT*5);
b = squeeze(sqrt(sum((XMot(:,3,:) - XMot(:,4,:)).^2)));
b = b / max(b);

for i = 1:10:size(XMot,3)
    subplot(6,1,1:4);
    cla;
    hold on;
    L = line([XMot(1,connA(:,1),i); XMot(1,connA(:,2)+ns,i)],...
         [XMot(2,connA(:,1),i); XMot(2,connA(:,2)+ns,i)],...
         [XMot(3,connA(:,1),i); XMot(3,connA(:,2)+ns,i)],...
         'color', [0 0 0 .5]);
%     set(L, {'LineWidth'}, kTC(1:size(connA,1)));
    plot3(XMot(1,1:ns,i), XMot(2,1:ns,i), XMot(3,1:ns,i),...
        'ro', 'linewidth', 15)
    plot3(XMot(1,(1:nu)+ns,i), XMot(2,(1:nu)+ns,i), XMot(3,(1:nu)+ns,i),...
        'bo', 'linewidth', 15);
    hold off;
    set(gca, 'visible', 'off');
    axis([-2.0 2.2 -1.5 1.5 -1 1]);
%     view(-5, 20);
%     view(0, 90);
%     view(0,0);
    
    subplot(6,1,5:6);
    plot(EA(1:i), 'linewidth', 4);
    hold on;
    plot(b*max(EA));
%     plot(EK);
    hold off;
%     axis([0, length(EA) 0 max([EA EK])]);
    xlabel('t'); ylabel('Energy');
    set(gca, 'xtick', [], 'ytick', [], 'fontsize', 24);
    drawnow;
end



%% Figure 6
figure(6); clf;
nPan1 = [];
for i = 1:5; nPan1 = [nPan1, [1:4]+(i-1)*24]; end
nPan2 = [];
for i = 1:6; nPan2 = [nPan2, [1:24]+(i-1)*24]; end
az = 30; el = 25;         % Views
azC = 50; elC = 30;         % Camera Position

% Test
x0 = [-1 -1 1 1;...
       -1 1 1 -1]/2;
xF = [-1 -1 1 1;...
       -1 1 1 -1]*0.8;
subplot(2,3,1);
visualize_conic_finite(x0, xF, [-2 2; -2 2], [100; 100], 8, 1, 1);


xS1 = [-2 -2 -1 -1;...
       -1  1  0  0;...
        0  0 -1  1]/2;
xF1 = xS1 - [ 0     0     0     0;...
             -0.2   0.2   0     0;...
              0     0    -0.06  0.06];
        
xS2 = [ 2  2  1  1;...
       -1  1  0  0;...
        0  0 -1  1]/2;
xF2 = xS2 - [ 0     0     0     0;...
             -0.2   0.2   0     0;...
              0     0    -0.06  0.06];

          

% Nodes internal mid causes initial downward deflections
% x01 = [-1.273  0.445  0.838;...
%        -1.273 -0.445  0.838;...
%        -1.273  0.449 -0.859;...
%        -1.273 -0.449 -0.859;...
%        -0.431  0.000  0.859;...
%        -0.431  0.000 -0.859]';
% x01 = x01 + (rand(size(x01))-.5)*.1;

% Works but is bulky
% x01 = [-1.561 -0.876  0.415;...
%        -1.561  0.876  0.315;...
%        -1.390 -0.955  0.10;...
%        -1.390  0.955  0.20;...
%        -1.476 -0.700 -0.233;...
%        -1.476  0.700 -0.333]';

x02 = [-x01(1,:); x01(2:3,:)];
conn1 = [1 1; 2 1; 3 1; 4 1;...
         1 2; 2 2; 3 2; 4 2;...
         1 3; 2 3; 3 3; 4 3;...
         1 4; 2 4; 3 4; 4 4;...
         1 5; 2 5; 3 5; 4 5;...
         1 6; 2 6; 3 6; 4 6];
     
xC = [-1 -1  1  1  0      ;...
       0  0  0  0  sqrt(2);...
      -1  1 -1  1  0]/2   ;
x0C = [ 0        sqrt(2)  sqrt(2);...
        0        sqrt(2) -sqrt(2);...
        sqrt(2)  sqrt(2)  0.0    ;...
       -sqrt(2)  sqrt(2)  0.0     ]'/4;
uC = xC/2;
connC = [1 1; 2 1; 3 1; 4 1; 5 1;...
         1 2; 2 2; 3 2; 4 2; 5 2;...
         1 3; 2 3; 3 3; 4 3; 5 3;...
         1 4; 2 4; 3 4; 4 4; 5 4];

     
     
subplot(2,3,1); cla;
visualize_conic_finite(xS1, xF1, [-4 2; -1 1; -1 1]*1, [100; 100; 100], 100, 1, 1);
visualize_conic(xS1, xF1-xS1 + [1 1 0 0; 0 0 0 0; 0 0 0 0]/4, [-4 2; -1 1; -1 1]*1, [100; 100; 100], 0, 1, 1);



subplot(2,3,1);
visualize_conic_finite(xS1, xF1, [-2 2; -1 1; -1 1]*1, [100; 100; 100], 10, 1, 1);
view(az, el);
camlight(azC, elC); lighting gouraud; material([.4 1 0]);
subplot(2,3,2);
visualize_conic_finite(xS2, xF2, [-2 2; -1 1; -1 1]*1, [100; 100; 100], 10, 1, 1);
view(az, el);
camlight(azC, elC); lighting gouraud; material([.4 1 0]);
subplot(2,3,3);
[Xu1, fV] = construct_network_finite(xS1, xF1, x01, conn1, 1);
subplot(2,3,4);
[Xu2, fV] = construct_network_finite(xS2, xF2, x02, conn1, 1);
subplot(2,3,5);
visualize_conic(xC, uC, [-1 1; -1 1; -1 1]*2, [100; 100; 100], 10, 1, 1);
[xCu, fV] = construct_network(xC, uC, x0C, connC, 1);
view(az, el);
camlight(azC, elC); lighting gouraud; material([.4 1 0]);
subplot(2,3,6);
[XsT, XuT, connT] = tesselate_network([xS1 xS2 xC], [Xu1(1:3,:) Xu2(1:3,:) xCu], [conn1; conn1+[4 6]; connC+[8 12]], [0 0 0]', [1 1 1]');
% [XV, FV] = sim_motion3D_con(XsT, XuT, connT, 0.001, 800, [1 2], xF1(:,[1 2]));
LVal = sqrt((XsT(1,connT(:,1)) - XuT(1,connT(:,2))).^2 +...
            (XsT(2,connT(:,1)) - XuT(2,connT(:,2))).^2 +...
            (XsT(3,connT(:,1)) - XuT(3,connT(:,2))).^2);

LVal1 = sqrt((xS1(1,conn1(:,1)) - Xu1(1,conn1(:,2))).^2 +...
             (xS1(2,conn1(:,1)) - Xu1(2,conn1(:,2))).^2 +...
             (xS1(3,conn1(:,1)) - Xu1(3,conn1(:,2))).^2);
         
LVal2 = sqrt((xS2(1,conn1(:,1)) - Xu2(1,conn1(:,2))).^2 +...
             (xS2(2,conn1(:,1)) - Xu2(2,conn1(:,2))).^2 +...
             (xS2(3,conn1(:,1)) - Xu2(3,conn1(:,2))).^2);
drawnow;
         
         
% Simulate Single Module
kT = [ones(48,1); 1*ones(23,1)];
[MSi, ESi] = sim_motion3D_congrad(xS2, Xu2(1:3,:), conn1, LVal2, 0.001, 1000, [1 2], xF2(:,[1 2]));
[MSSi, ESSi] = sim_motion3D_congrad(MSi(:,1:4,end), MSi(:,5:10,end), conn1, LVal1, 0.001, 20, [1], MSi(:,[1]));
D = zeros(size(xS1,2), size(xS1,2), length(ESi));
for i = 1:length(ESi)
    D(:,:,i) = squareform(pdist(MSi(:,1:4,i)'));
end
DS = squareform(pdist(xF1'));
DSS = squareform(pdist(MSSi(:,:,end)'));
plot(squeeze(sum(sum((D - DS).^2)))');

% Simulate Original ==> Second Deformation
[MK, EK] = sim_motion3D_congrad(XsT, XuT, connT, LVal, 0.001, 1000, [1 2], xF1(:,[1 2]), kT);
% Simulate Original ==> First Deformation
[M, E] = sim_motion3D_congrad(XsT, XuT, connT, LVal, 0.001, 1000, [5 6], xF2(:,[1 2]), kT);
XsT2 = M(:,1:size(XsT,2),end);
XuT2 = M(:,[1:size(XuT,2)]+size(XsT,2),end);
% Simulate First Deformation ==> Stablize
[MS, ES] = sim_motion3D_congrad(XsT2, XuT2, connT, LVal, 0.001, 2, [1], XsT2(:,[1]), kT);
XsT3 = MS(:,1:size(XsT,2),end);
XuT3 = MS(:,[1:size(XuT,2)]+size(XsT,2),end);
% Simulate First Deformation ==> Second Deformation
[MS2, ES2] = sim_motion3D_congrad(XsT3, XuT3, connT, LVal, 0.001, 1000, [1 2], xF1(:,[1 2]), kT);

figure(1); clf;
subplot(3,1,1);
plot(E);
subplot(3,1,2);
plot(ES);
subplot(3,1,3);
plot(EK);
hold on;
plot(ES2);
hold off;


%% Video
XMot = zeros(3,25,1);
XMot(:,:,1:size(M,3)) = M;
XMot(:,:,[1:size(MS2,3)]+size(M,3)) = MS2;
connA = connT;
EA = [E ES2];
XsA = XsT;
XuA = XuT;

% XMot = MK;
% EA = EK;
% connA = connT;
% XsA = XsT;
% XuA = XsT;

% XMot = MSi;
% EA = ESi;
% connA = conn1;
% XsA = xS1;
% XuA = Xu1;

h = figure(7); clf;
filename = 'cooperativity3D.gif';
ns = size(XsA,2);
nu = size(XuA,2);

kTC = num2cell(kT*5);

for i = 1:20:size(XMot,3)
    subplot(6,1,1:4);
    cla;
    hold on;
    L = line([XMot(1,connA(:,1),i); XMot(1,connA(:,2)+ns,i)],...
         [XMot(2,connA(:,1),i); XMot(2,connA(:,2)+ns,i)],...
         [XMot(3,connA(:,1),i); XMot(3,connA(:,2)+ns,i)],...
         'color', [0 0 0 .5]);
    set(L, {'LineWidth'}, kTC(1:size(connA,1)));
    plot3(XMot(1,1:ns,i), XMot(2,1:ns,i), XMot(3,1:ns,i),...
        'ro', 'linewidth', 15)
    plot3(XMot(1,(1:nu)+ns,i), XMot(2,(1:nu)+ns,i), XMot(3,(1:nu)+ns,i),...
        'bo', 'linewidth', 15);
%         plot3(Xs(1,:), Xs(2,:), Xs(3,:), 'ro', 'markersize', 14, 'linewidth', 2);
%     plot3(XsT(1,:)+Us(1,:), Xs(2,:)+Us(2,:), Xs(3,:)+Us(3,:), 'ro', 'markersize', 14, 'linewidth', 2);
    hold off;
    set(gca, 'visible', 'off');
    axis([-2.0 2.2 -1.5 1.5 -1 1]);
%     axis([-2 1 -1 1 -1 1]);
%     view(-5, 20);
%     view(0, 90);
%     view(0,0);
    
    subplot(6,1,5:6);
    plot(EA(1:i), 'linewidth', 4);
    hold on;
    plot(squeeze(sqrt(sum((XMot(:,3,:) - XMot(:,4,:)).^2)))/10000)
    plot(EK);
    hold off;
    axis([0, length(EA) 0 max([EA EK])]);
    xlabel('t'); ylabel('Energy');
    set(gca, 'xtick', [], 'ytick', [], 'fontsize', 24);
    drawnow;
    
%     frame = getframe(h);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%     
%     if i == 1
%         imwrite(imind,cm,filename,'gif','Loopcount',inf,'DelayTime',.03);
%     else
%         imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',.03);
%     end
    
end









