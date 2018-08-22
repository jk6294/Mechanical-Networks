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
Us = [[-1 1];...
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
az = -160; el = 20;         % Views


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
x01 = [[-sqrt(2) sqrt(2)];...
       [0 0]];
x02 = [[-sqrt(2) 0 sqrt(2) 0];...
       [0 sqrt(2) 0 -sqrt(2)]];

   
% a: Visualize Conics
subplot(12,24,nPan1)
visualize_conic(Xs, Us, [-1 1; -1 1]*2, [100 100], 0, 1, 0);
axis(1.5*[-1 1 -1 1]);


% b: Find positions of unspecified nodes, and corresponding motions
subplot(12,24,nPan1+8)
[Xu1 fV1] = construct_network(Xs, Us, x01, conn1, 0);
[Us1, Uu1, err1] = construct_motion(Xs, Us, Xu1, conn1, 1, 1);
axis(1.5*[-1 1 -1 1]);


% c:  Find positions of unspecified nodes, and corresponding motions
subplot(12,24,nPan1+16)
[Xu2 fV2] = construct_network(Xs, Us, x02, conn2, 0);
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
x01 = [[-sqrt(3)  sqrt(3)  0        0        0        0      ];...
       [0         0       -sqrt(3)  sqrt(3)  0        0      ];...
	   [0         0        0        0       -sqrt(3)  sqrt(3)]]*.8;
x02 = [[-sqrt(3)  sqrt(3)  0        0        0      ];...
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
[Xu1 fV1] = construct_network(Xs, Us, x01, conn1, 0);
[Us1, Uu1, err1] = construct_motion(Xs, Us, Xu1, conn1, 1, 1);
axis(1.5*[-1 1 -1 1 -1 1]);
view(az, el);
camlight(azC, elC); lighting gouraud; material([.4 1 0]);


% f: Find positions of unspecified nodes, and corresponding motions
subplot(12,24,nPan2+24*6+16)
[Xu2 fV2] = construct_network(Xs, Us, x02, conn2, 0);
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
for i = 1:6; nPan2 = [nPan2, [1:8]+(i-1)*24]; end
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
UsSSS(:,:,2) = [[.2 -1 1 -.2];...
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


% b: Place unspecified nodes and construct motions
subplot(12,24,nPan1+8)
[Xu1, fV] = construct_network(Xs, Us, x0, conn, 0);
[Us1, Uu1, err] = construct_motion(Xs, Us, Xu1, conn, 1, 1);
axis(1.5*[-1 1 -1 1]);


% c: Construct states of self-stress
subplot(12,24,nPan1+16);
[XuSSS1, fV] = construct_network(XsSSS(:,[1 2 3]), UsSSS(:,[1 2 3],:), [-sqrt(2); 0], connSSS, 0);
[XuSSS2, fV] = construct_network(XsSSS(:,[1 2 4]), UsSSS(:,[1 2 4],:), [-.5; 1.3], connSSS, 0);
[XuSSS4, fV] = construct_network(XsSSS(:,[2 3 4]), UsSSS(:,[2 3 4],:), [sqrt(2); 0], connSSS, 0);
[XuSSS3, fV] = construct_network(XsSSS(:,[1 3 4]), UsSSS(:,[1 3 4],:), [.5; 1.3], connSSS, 0);
XuSSS = [XuSSS1 XuSSS2 XuSSS3 XuSSS4];
connSSS = [1 1; 2 1; 3 1; 1 2; 2 2; 4 2; 1 3; 3 3; 4 3; 2 4; 3 4; 4 4];
[Us1, Uu1, err] = construct_motion(XsSSS, UsSSS, XuSSS, connSSS, 1, 1);


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
% Connectivity
conn = [1 1; 2 1; 3 1; 4 1; 1 2; 2 2; 3 2; 4 2;...
        1 3; 2 3; 3 3; 4 3; 1 4; 2 4; 3 4; 4 4];
    
% d: Visualize solution space
subplot(12,24,nPan2+24*6)
visualize_conic(Xs, Us, [-1 1; -1 1; -1 1]*1.8, [100 100 100], 0, 1, 0);
axis(1.5*[-1 1 -1 1 -1 1]);
view(az, el);
camlight(azC, elC); lighting gouraud; material([.4 1 0]);


% e: Place unspecified nodes and construct motion
subplot(12,24,nPan2+24*6+8)
[Xu1, fV] = construct_network(Xs, Us, x0, conn, 0);
[Us1, Uu1, err] = construct_motion(Xs, Us, Xu1, conn, 1, 1);
axis(1.5*[-1 1 -1 1 -1 1]);
view(az, el);
camlight(azC, elC); lighting gouraud; material([.4 1 0]);


%% Figure 5b-d
figure(5); clf;
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


% b: Combined modules with motion
subplot(12,24,nPan1 + 4)
[Xu, fV] = construct_network(Xs, Us, x0, conn, 0);
[XsT, XuT, connT] = tesselate_network(Xs, Xu, conn, [2 2]', [1 2]');
[Us, Uu, err] = construct_motion(XsT, zeros(0,0,0), XuT, connT, -2, -2);


% c: Tesselation expanded
subplot(12,24,[nPan1 nPan1+4]+8)
[XsT, XuT, connT] = tesselate_network(Xs, Xu, conn, [2 2]', [5 5]');
[XMot1, fC] = sim_motion(XsT, XuT, connT, .1, 130, [XsT XuT],0);
construct_motion(XMot1(:,1:size(XsT,2),end), zeros(0,0,0), XMot1(:,[1:size(XuT,2)]+size(XsT,2),end), connT, .1, .1);
axis([-4 12 -2 10]);


% d: Tesselation contracted
subplot(12,24,[nPan1 nPan1+4]+16)
[XMot2, fC] = sim_motion(XsT, XuT, connT, .1, 110, -[XsT XuT],0);
construct_motion(XMot2(:,1:size(XsT,2),end), zeros(0,0,0), XMot2(:,[1:size(XuT,2)]+size(XsT,2),end), connT, .1, .1);
axis([-4 12 -2 10]);


%% 3D
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
x01 = [[3.092 0.808 2.183];...
       [2.758 -0.930 1.516];...
       [2.758 0.930 1.516];...
       [2 -.77 2.061];...
       [2 .77 2.061]]';
x02 = [[5.081 -0.648 3.313];...
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
[Xu1 fV] = construct_network(Xs1, Us1, x01, conn, 0);
[Usp, Uup, err] = construct_motion(Xs1, -Us1, Xu1, conn, 1, 1);

[Xu2 fV] = construct_network(Xs2, Us2, x02, conn, 0);
[Usp, Uup, err] = construct_motion(Xs2, -Us2, Xu2, conn, 1, 1);

axis([.9 19.1 -3 1.4 .5 3.5]);
view(-10, 10);


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


%% Figure 6
