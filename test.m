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
subplot(11,24,nPan1)
Xs = [[-1 1 0];...
      [ 0 0 1]];
Us = [[-1 1 -1];...
      [-1 -1 0]]/4;
visualize_conic(Xs, Us, [-1 1; -1 1]*2, [100; 100], 14, 1, 1);
axis(1.4*[-1 1 -1 1]);


% e: Positions and Motions
subplot(11,24,nPan2+24*5)
Xs = [[-1 1 0 0 0];...
	  [-1 -1 1 0 0];...
	  [0 0 0 1 -1]]*.6;
Us = [[1 1 -1 0 0 ];...
	  [1 -1 0.1 0 0 ];...
	  [-1 -1 -1 1 1]]*.25;
visualize_conic(Xs, Us, [-1 1; -1 1; -1 1]*4, [100 100 100], 40, 1, 1);
axis(1.1*[-1 1 -1 1 -1 1]);
view(az, el);
camlight(azC, elC); lighting gouraud; material([.4 1 0]);


% f: Positions and Motions
subplot(11,24,nPan1+8)
Xs = [[-1 1];...
      [0 0]];
Us = -[[-1 1];...
      [0 0]]/4;
visualize_conic(Xs, Us, [-1 1; -1 1]*1.7, [100 100], 7, 1, 1);
axis(1.8*[-1 1 -1 1]);


% g: Positions and Motions
subplot(11,24,nPan2+24*5+8)
Xs = [[-1 1 0 0];...
      [-1 -1 1 0];...
      [0 0 0 1]]*.6;
Us = [[-1 1 0 0]/2;...
      [-1 -1 1 0]/4;...
      [-1 -1 -1 -1.6]]*.2;
visualize_conic(Xs, Us, [-1 1; -1 1; -1 1]*1.1, [100 100 100], 32, 1, 1);
axis(1.1*[-1 1 -1 1 -1 1]);
view(az, el);
camlight(azC, elC); lighting gouraud; material([.4 1 0]);


% h: Positions and Motions
subplot(11,24,nPan1+16)
Xs = [[-1 -1 1 1];...
      [-1 1 1 -1]];
Us = [[-1 -1 1 1];...
      [-1 1 1 -1]]/4;
visualize_conic(Xs, Us, [-1 1; -1 1]*2, [100 100], 16, 1, 1);
axis(1.6*[-1 1 -1 1]);


% i: Positions and Motions
subplot(11,24,nPan2+24*5+16)
Xs = [[-1 -1 -1 -1 1 1 1 1];...
      [-1 -1 1 1 -1 -1 1 1];...
	  [-1 1 -1 1 -1 1 -1 1]]*.6;
% Motions: 1D Manifold
Us = [[-1 -1 -1 -1 1 1 1 1];...
	  [-1 -1 1 1 -1 -1 1 1];...
	  [-1 1 -1 1 -1 1 -1 1]]*.2;
visualize_conic(Xs, Us, [-1 1; -1 1; -1 1]*1.1, [100 100 100], 40, 1, 1);
axis(1.1*[-1 1 -1 1 -1 1]);
view(az, el);
camlight(azC, elC); lighting gouraud; material([.4 1 0]);


% Size and Save Figure
fName = 'bipartite';
set(gcf, 'Renderer', 'opengl'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6.6 3.8];
fig.PaperSize = [6.6 3.8];
% print(['Figures\' fName], '-dpng','-r1500');


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
         1 2; 2 2; 3 2; 4 2] + [0 4];
conn2 = [1 2; 1 3; 2 1; 2 2; 2 3; 2 4;...
         3 1; 3 4; 4 1; 4 2; 4 3; 4 4] + [0 4];
% Initial guesses for unspecified node positions
xS1 = [[-sqrt(2) sqrt(2)];...
       [0 0]];
xS2 = [[-sqrt(2) 0 sqrt(2) 0];...
       [0 sqrt(2) 0 -sqrt(2)]];

   
% a: Visualize Conics
subplot(11,24,nPan1)
visualize_conic(Xs, Us, [-1 1; -1 1]*2, [100 100], 0, 1, 0);
axis(1.5*[-1 1 -1 1]);


% b: Find positions of unspecified nodes, and corresponding motions
subplot(11,24,nPan1+8)
[Xu1] = construct_network(Xs, Us, xS1, conn1, 0, 0);
construct_motion(Xs, Us, Xu1(1:2,:), conn1, 1, 1);
axis(1.5*[-1 1 -1 1]);


% c:  Find positions of unspecified nodes, and corresponding motions
subplot(11,24,nPan1+16)
[Xu2] = construct_network(Xs, Us, xS2, conn2, 0, 0);
construct_motion(Xs, Us, Xu2(1:2,:), conn2, 1, 1);
axis(1.5*[-1 1 -1 1]);


% Positions and Motions
Xs = [[-1 -1 -1 -1 1 1 1 1];...
	  [-1 -1 1 1 -1 -1 1 1];...
	  [-1 1 -1 1 -1 1 -1 1]]*.55;
Us = [[-1 -1 -1 -1 1 1 1 1];...
	  [-1 -1 1 1 -1 -1 1 1];...
	  [-1 1 -1 1 -1 1 -1 1]]*.3;
% Two different connectivities
conn1 = [1 1; 2 1; 3 1; 4 1; 5 1; 6 1;...
         1 2; 3 2; 5 2; 6 2; 7 2; 8 2;...
         1 3; 2 3; 5 3; 6 3; 3 3; 8 3;...
         3 4; 4 4; 5 4; 2 4; 7 4; 8 4;...
         1 5; 3 5; 5 5; 4 5; 7 5; 8 5;...
         2 6; 4 6; 5 6; 6 6; 8 6] + [0 8];
conn2 = [1 1; 2 1; 3 1; 4 1; 5 1; 6 1; 7 1; 8 1;...
         1 2; 2 2; 3 2; 4 2; 5 2; 6 2; 7 2; 8 2;...
         1 3; 2 3; 3 3; 5 3; 6 3;...
         3 4; 4 4; 5 4; 7 4; 8 4;...
         1 5; 2 5; 3 5; 5 5; 7 5; 8 5] + [0 8];
% Initial guesses for unspecified node positions
xS1 = [[-sqrt(3)  sqrt(3)  0        0        0        0      ];...
       [0         0       -sqrt(3)  sqrt(3)  0        0      ];...
	   [0         0        0        0       -sqrt(3)  sqrt(3)]]*.8;
xS2 = [[-sqrt(3)  sqrt(3)  0        0        0      ];...
	   [0         0       -sqrt(3)  sqrt(3)  0      ];...
	   [0         0        0        0       -sqrt(3)]]*.8;
   

% d: Visualize conic solution space
subplot(11,24,nPan2+24*5)
visualize_conic(Xs, Us, [-1 1; -1 1; -1 1]*1.8, [100 100 100], 0, 1, 1);
axis(1.1*[-1 1 -1 1 -1 1]);
view(az, el);
camlight(azC, elC); lighting gouraud; material([.4 1 0]);


% e: Find positions of unspecified nodes, and corresponding motions
subplot(11,24,nPan2+24*5+8)
[Xu1] = construct_network(Xs, Us, xS1, conn1, 0, 0);
construct_motion(Xs, Us, Xu1(1:3,:), conn1, 1, 1);
axis(1.1*[-1 1 -1 1 -1 1]);
view(az, el);
camlight(azC, elC); lighting gouraud; material([.4 1 0]);


% f: Find positions of unspecified nodes, and corresponding motions
subplot(11,24,nPan2+24*5+16)
[Xu2] = construct_network(Xs, Us, xS2, conn2, 0, 0);
[~, Uu2, err2] = construct_motion(Xs, Us, Xu2(1:3,:), conn2, 1, 1);
axis(1.1*[-1 1 -1 1 -1 1]);
view(az, el);
camlight(azC, elC); lighting gouraud; material([.4 1 0]);


% Size and Save Figure
fName = 'construction';
set(gcf, 'Renderer', 'opengl'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 8.2 4.1];
fig.PaperSize = [8.2 4.1];
% print(['Figures\' fName], '-dpng','-r1500');


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
conn = [1 1; 2 1; 3 1] + [0 3];
connSSS = [1 1; 2 1; 3 1] + [0 3];


% a: Visualize solution space
subplot(11,24,nPan1)
visualize_conic(Xs, Us, [-1 1; -1 1]*2, [100 100], 0, 1, 0);
axis(1.5*[-1 1 -.8 1]);


% c: Place unspecified nodes and construct motions
subplot(11,24,nPan1+8)
[Xu1] = construct_network(Xs, Us, x0, conn, 0, 0);
construct_motion(Xs, Us, Xu1(1:2,:), conn, 1, 1);
axis(1.5*[-1 1 -.8 1]);


% e: Construct states of self-stress
subplot(11,24,nPan1+16); cla;
[XuSSS1] = construct_network(XsSSS(:,[2 3 4]), UsSSS(:,[2 3 4],:), [-1 -1]', connSSS, 0, 0);
[XuSSS2] = construct_network(XsSSS(:,[1 3 4]), UsSSS(:,[1 3 4],:), [-1 0]', connSSS, 0, 0);
[XuSSS3] = construct_network(XsSSS(:,[1 2 4]), UsSSS(:,[1 2 4],:), [1 0]', connSSS, 0, 0);
[XuSSS4] = construct_network(XsSSS(:,[1 2 3]), UsSSS(:,[1 2 3],:), [1 -1]', connSSS, 0, 0);
XuSSS = [XuSSS1 XuSSS2 XuSSS3 XuSSS4];
connSSS = [2 1; 3 1; 4 1; 1 2; 3 2; 4 2; 1 3; 2 3; 4 3; 1 4; 2 4; 3 4] + [0 4];
construct_motion(XsSSS, UsSSS, XuSSS(1:2,:), connSSS, .01, .01);
axis(1.5*[-1 1 -1 1]);


% Positions and Motions
Xs = [[-1 0 1 0];...
      [-1 1 -1 0];...
	  [0 0 0 1]]*.68;
% Two separate motions
Us = zeros(3,4,2);
Us(:,:,1) = [[-1.0  0.0  1.0  0.0];...
             [-1.0  1.0 -1.0  0.0];...
             [0.0  0.0  0.0  1.0]]*.4;
Us(:,:,2) = [[0.7 -0.4  1.2 -0.2];...
             [1.0 -0.2 -0.5  0.9];...
             [-1.0  0.0  0.0  0.0]]*.4;
% Initial guess for unspecified nodes
x0 = [-.7 -.5 .8;...
      1.2 -.4 -.75;...
      -.4 -1.5 -.3;...
      .2 .6 .6]'*.6;
x0SSS = [-.7 -.5 .8;...
         1.2 -.4 -.75;...
         -.4 -1.5 -.3;...
         .2 .6 .6;...
         -.05 .58 -.85]'*.6;...
         
% Connectivity
conn = [1 1; 2 1; 3 1; 4 1;...
        1 2; 2 2; 3 2; 4 2;...
        1 3; 2 3; 3 3; 4 3;...
        1 4; 2 4; 3 4; 4 4] + [0 4];
connSSS = [1 1; 2 1; 3 1; 4 1;...
           1 2; 2 2; 3 2; 4 2;...
           1 3; 2 3; 3 3; 4 3;...
           1 4; 2 4; 3 4; 4 4;...
           1 5; 2 5; 3 5; 4 5] + [0 4];

    
% b: Visualize solution space
subplot(11,24,nPan2+24*5)
visualize_conic(Xs, Us, [-1 1; -1 1; -1 1]*1.8, [100 100 100], 0, 1, 0);
axis(1.2*[-1 1 -1 1 -1 1]);
view(az, el);
camlight(azC, elC); lighting gouraud; material([.4 1 0]);


% d: Place unspecified nodes and construct motion
subplot(11,24,nPan2+24*5+8)
[Xu1] = construct_network(Xs, Us, x0, conn, 0, 0);
construct_motion(Xs, Us, Xu1(1:3,:), conn, 1, 1);
axis(1.2*[-1 1 -1 1 -1 1]);
view(az, el);
camlight(azC, elC); lighting gouraud; material([.4 1 0]);


% f: Place unspecified nodes and construct motion
subplot(11,24,nPan2+24*5+16)
[Xu1] = construct_network(Xs, Us, x0SSS, connSSS, 0, 0);
% Lengths
[~, Uu1, ~] = construct_motion(Xs, Us, Xu1(1:3,:), connSSS, .01, .01);
axis(1.2*[-1 1 -1 1 -1 1]);
view(az, el);
camlight(azC, elC); lighting gouraud; material([.4 1 0]);


% Size and Save Figure
fName = 'multi';
set(gcf, 'Renderer', 'opengl'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 8.2 4.1];
fig.PaperSize = [8.2 4.1];
% print(['Figures\' fName], '-dpng','-r1500');


%% Figure 5
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
conn = [1 1; 2 1; 3 1; 4 1; 1 2; 2 2; 3 2; 4 2] + [0 4];


% a: Independent modules
subplot(11,24,nPan1)
construct_network(Xs + [0;1.2], Us, x0 + [0;1.2], conn, 1, 0);
construct_network(Xs - [0;1.2], Us, x0 - [0;1.2], conn, 1, 0);
axis([-1.8 1.8 -2.8 2.8]);


% b: Combined modules with motion
subplot(11,24,nPan1 + 4)
[Xu] = construct_network(Xs, Us, x0, conn, 0, 0);
[XsT, XuT, connT] = tesselate_network(Xs, Xu(1:2,:), conn, [2 2]', [1 2]');
[Us, Uu] = construct_motion(XsT, zeros(0,0,0), XuT, connT, -2, -2);
axis([-1.8 1.8 -1.8 3.8]);


% c: Tesselation expanded
subplot(11,24,[nPan1 nPan1+4]+8)
[XsT, XuT, connT] = tesselate_network(Xs, Xu(1:2,:), conn, [2 2]', [5 5]');
[XMot1] = sim_motion(XsT, XuT, connT, .1, 130, [XsT XuT],0);
construct_motion(XMot1(:,1:size(XsT,2),end), zeros(0,0,0), XMot1(:,(1:size(XuT,2))+size(XsT,2),end), connT, .1, .1);
axis([-4 12 -2 10]);


% d: Tesselation contracted
subplot(11,24,[nPan1 nPan1+4]+16);
[XMot2] = sim_motion(XsT, XuT(1:2,:), connT, .1, 130, -[XsT XuT],0);
construct_motion(XMot2(:,1:size(XsT,2),end), zeros(0,0,0), XMot2(:,(1:size(XuT,2))+size(XsT,2),end), connT, .1, .1);
axis([-4 12 -2 10]);


% 3D
subplot(11,24,nPan2+24*6); cla;

% Positions
Xs1 = [[2 1 3 3.5 3.5];...
       [0 0 0 -1 1];...
       [1 2 3 2.5 2.5]]*.8;
Xs2 = [[5 6 4 3.5 3.5];...
       [0 0 0 -1 1];....
       [1 2 3 2.5 2.5]]*.8;
% Motions
Us1 = [[1 -1 1 0 0];...
       [0 0 0 1 -1];...
       [-1 1 0 0 0]]*.3;
Us2 = [[1 -1 1 0 0];...
       [0 0 0 1 -1];...
       [1 -1 0 0 0]]*.3;
% Initial guess for unspecified nodes
xS1 = [[3.092 0.808 2.183];...
       [2.758 -0.930 1.516];...
       [2.758 0.930 1.516];...
       [2 -.77 2.061];...
       [2 .77 2.061]]'*.8;
xS2 = [[4.081 -0.648 3.313];...
       [4.081 0.648 3.313];...
       [5.184 -1.000 1.737];...
       [4.121 -1.287 2.626];...
       [4.121 1.287 2.626]]'*.8;
% Connectivity
conn = [1 1; 2 1; 3 1; 4 1; 5 1;...
        1 2; 2 2; 3 2; 4 2; 5 2;...
        1 3; 2 3; 3 3; 4 3; 5 3;...
        1 4; 2 4; 3 4; 4 4;...
        1 5; 2 5; 3 5; 5 5] + [0 5];

% Simulate Motion
[Xu1] = construct_network(Xs1, Us1, xS1, conn, 0, 0);
construct_motion(Xs1, -Us1, Xu1(1:3,:), conn, 1, 1);

[Xu2] = construct_network(Xs2, Us2, xS2, conn, 0, 0);
[Usp, Uup, err] = construct_motion(Xs2+[.8 0 0]', -Us2, Xu2(1:3,:)+[.8 0 0]', conn, 1, 1);

[XsT, XuT, connT] = tesselate_network([Xs1 Xs2], [Xu1(1:3,:) Xu2(1:3,:)], [conn; conn+[size(Xs1,2), size(Xu1,2)]], [0 0 0]', [1 1 1]');
[XMot] = sim_motion3D(XsT+[5.7 0 0]', XuT(1:3,:)+[5.75 0 0]', [connT; 3 8], .01, 70, ([XsT XuT]+[5.7 0 0]'), 0);
construct_motion(XMot(:,1:size(XsT,2),end), zeros(0,0,0), XMot(:,(1:size(XuT,2))+size(XsT,2),end), [connT; 3 8], .01, .01);

[XsT, XuT, connT] = tesselate_network([Xs1 Xs2], [Xu1(1:3,:) Xu2(1:3,:)], [conn; conn+[size(Xs1,2), size(Xu1,2)]], [0 0 0]', [1 1 1]');
[XMot] = sim_motion3D(XsT+[10.4 0 0]', XuT(1:3,:)+[10.4 0 0]', [connT; 3 8], .01, 110, -([XsT XuT]+[10 0 0]'), 0);
construct_motion(XMot(:,1:size(XsT,2),end), zeros(0,0,0), XMot(:,(1:size(XuT,2))+size(XsT,2),end), [connT; 3 8], .01, .01);

% Effects
axis([.5 15.4 -3.1 1.1 .5 2.6]);
view(-7, 20);
delete(findall(gcf,'Type','light'))
h = light;
h.Position = [2, -10, 10];
lighting gouraud; 
material([.4 1 0]);

% Size and Save Figure
fName = 'combine';
set(gcf, 'Renderer', 'opengl'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 8.2 3.5]*1.4;
fig.PaperSize = [8.2 3.5]*1.4;
% print(['Figures\' fName], '-dpng','-r1500');


%% Figure 6
fig = figure(6); clf;
% Subplot Panel Sizes
nPan1 = [];
for i = 1:4; nPan1 = [nPan1, [1:4]+(i-1)*24]; end
nPan2 = [];
for i = 1:3; nPan2 = [nPan2, [1:8]+(i-1)*24]; end
nPan3 = [];
for i = 1:4; nPan3 = [nPan3, [1:10]+(i-1)*24]; end
az = 7; el = 14;                % Views
C_SA = [76 187 23;...           % Color of Different Motions
        50 255 50]/255;         


% a: 2 Dimensional Solution Space
subplot(12,24,nPan1); cla;
Xs0 = [-1 -1 1 1;...
       -1 1 1 -1];
XsF = Xs0 + [-1 -1 1 1;...
             -1 1 1 -1]*1;
visualize_conic_finite(Xs0, XsF, [-1 1; -1 1]*3, [100; 100],4, .75, .85);
axis(3.0*[-1 1 -1 1]);


% b: 2 Dimensional Motion Example
subplot(12,24,nPan1+24*4); cla;
% Initial guess for unspecified nodes
x0 = [-sqrt(2) 0;...
       sqrt(2) 0]'*2;
% Connectivity
conn = [1 1; 2 1; 3 1; 4 1; 1 2; 2 2; 3 2; 4 2] + [0 4];
% Construct and Simulate Motion
[Xu] = construct_network(Xs0, XsF, x0, conn, 0, 1);
[XMot, fC] = sim_motion(Xs0, Xu(1:2,:), conn, .01, 400, [Xs0 Xu(1:2,:)],1);
hold on;
plot(XsF(1,:), XsF(2,:), 'wo', 'markersize', 4, 'linewidth', 4);
plot(Xu(3,:), Xu(4,:), 'wo', 'markersize', 4, 'linewidth', 4);
plot(XsF(1,:), XsF(2,:), 'o', 'markersize', 7, 'linewidth', 1, 'color', [255 100 100]/255);
plot(Xu(3,:), Xu(4,:), 'o', 'markersize', 7, 'linewidth', 1, 'color', [100 100 255]/255);
hold off;
axis(3.0*[-1 1 -1 1]);


% c: 3 Dimensional Solution Space
subplot(12,24,nPan1+24*8); cla;
xS1 = [-3   -3   -1   -1;...
       -1    1    0    0;...
        0    0   -1    1]/2;
xF1 = xS1 + [ 0.2   0.2  -0.15  -0.15;...
             -0.2   0.2   0     0;...
              0     0    -0.15   0.15];
visualize_conic_finite(xS1, xF1, [-1.8 -0.0; -.3 .8; -.75 .75]*1, [100; 100; 100], 8, 1, 1);
axis(1*[-1.8 -0.0 -.9 .9 -.8 .8]);
view(az, el);
lighting gouraud; 
lightangle(45, 45);
material([.4 1 0]);


% d: 3 Dimensional initial and final positions, top
subplot(12,24,nPan2+5); cla;
x01 = [-0.8483  0.2323  0.3333;...
       -0.8483 -0.2323  0.3333;...
       -0.8111  0.2564 -0.3778;...
       -0.8111 -0.2564 -0.3778;...
       -1.1489  0.4317  0.2590;...
       -1.1489 -0.4317  0.2590]';
x01 = x01 + [0.01  0.00  0.02  0.00  0.01  0.00;...
            -0.01  0.00  0.00  0.00  0.00  0.00;...
             0.01  0.00  0.00  0.00  0.06  0.00]/2;
conn1 = [1 1; 2 1; 3 1; 4 1;...
         1 2; 2 2; 3 2; 4 2;...
         1 3; 2 3; 3 3; 4 3;...
         1 4; 2 4; 3 4; 4 4;...
         1 5; 2 5; 3 5; 4 5;...
         1 6; 2 6; 3 6; 4 6] + [0 4];
[Xu1] = construct_network(xS1, xF1, x01, conn1, 0, 1);
visualize_network(xS1-[1.5 0 0]', Xu1(1:3,:)-[1.5 0 0]', conn1);
visualize_network(xF1+[1.5 0 0]', Xu1(4:6,:)+[1.5 0 0]', conn1);
hold on;
quiver3(xS1(1,[1,2])-1.5,xS1(2,[1,2]),xS1(3,[1,2]),...
        xF1(1,[1,2])-xS1(1,[1,2]),xF1(2,[1,2])-xS1(2,[1,2]),xF1(3,[1,2])-xS1(3,[1,2]),0,...
        'linewidth', 1, 'color', C_SA(1,:));
hold off;
axis([-3.1 1.0 -0.8 0.8 -0.8 0.8]);
view(0, 90);
lighting gouraud;
lightangle(45, 45);
material([.4 1 0]);
set(gca, 'visible', 0);


% e: 3 Dimensional Single Module Energy
subplot(12,24,nPan2+24*3+5); cla;
LVal1 = sqrt((xS1(1,conn1(:,1)) - Xu1(1,conn1(:,2)-4)).^2 +...
             (xS1(2,conn1(:,1)) - Xu1(2,conn1(:,2)-4)).^2 +...
             (xS1(3,conn1(:,1)) - Xu1(3,conn1(:,2)-4)).^2)';
hold on;
[M1, E1] = sim_motion3D_congrad(xS1, Xu1(1:3,:), conn1, LVal1, 0.001, 1001, [1 2], xF1(:,[1 2]));
plot(E1, 'linewidth',1);
hold off;
set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);
box off;
axis([0 length(E1) 0 max(E1)]);


% f: 3 Dimensional initial and final positions, side
subplot(12,24,nPan2+24*6+5); cla;
visualize_network(xS1-[1.5 0 0]', Xu1(1:3,:)-[1.5 0 0]', conn1);
visualize_network(xF1+[1.5 0 0]', Xu1(4:6,:)+[1.5 0 0]', conn1);
hold on;
h = quiver3(xS1(1,[3,4])-1.5,xS1(3,[3,4]),xS1(2,[3,4]),...
            xF1(1,[3,4])-xS1(1,[3,4]),xF1(3,[3,4])-xS1(3,[3,4]),xF1(2,[3,4])-xS1(2,[3,4]),0,...
            'linewidth', 1, 'color', C_SA(1,:));
t = hgtransform('Parent',gca);
R = makehgtform('xrotate',pi/2);
set(t,'Matrix',R);
set(h,'Parent',t);      % The arrows should point correctly now
hold off;
axis([-3.1 1.0 -0.8 0.8 -0.8 0.8]);
view(0, 0);
lighting gouraud;
lightangle(45, 45);
material([.4 1 0]);
set(gca, 'visible', 0);


% g: 3 Dimensional Single Module Energy
subplot(12,24,nPan2+24*9+5); cla;
[M2, E2] = sim_motion3D_congrad(xS1, Xu1(1:3,:), conn1, LVal1, 0.001, 1001, [3 4], xF1(:,[3 4]));
plot(E2, 'linewidth', 1);
set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);
box off;
axis([0 length(E1) 0 max(E2)]);


% h: Coupling Modules
subplot(12,24,nPan3+14); cla;
% Position Mirror Image
xS2 = [-xS1(1,:); xS1(2:3,:)];
xF2 = [-xF1(1,:); xF1(2:3,:)];
Xu2 = [-Xu1(1,:); Xu1(2:3,:)];
% Coupling Components
xC = [-1 -1  1  1  0      ;...
       0  0  0  0  sqrt(2);...
      -1  1 -1  1  0]/2   ;
x0C = [ 0       -sqrt(2)  sqrt(2);...
        0       -sqrt(2) -sqrt(2);...
        sqrt(2) -sqrt(2)  0.0    ;...
       -sqrt(2) -sqrt(2)  0.0     ]'/4;
uC = xC/2;
connC = [1 1; 2 1; 3 1; 4 1; 5 1;...
         1 2; 2 2; 3 2; 4 2; 5 2;...
         1 3; 2 3; 3 3; 4 3; 5 3;...
         1 4; 2 4; 3 4; 4 4; 5 4] + [0 5];
[xCu, fV] = construct_network(xC, uC, x0C, connC, 1, 0);
[XsT, XuT, connT] = tesselate_network([xS1 xS2 xC], [Xu1(1:3,:) Xu2(1:3,:) xCu(1:3,:)], [conn1+[0 9]; conn1+[4 15]; connC+[8 20]], [0 0 0]', [1 1 1]');
visualize_network(xS1-[.6 0 0]', Xu1(1:3,:)-[.6 0 0]', conn1);
visualize_network(xS2+[.6 0 0]', Xu2(1:3,:)+[.6 0 0]', conn1);
axis([-2.4 2.4 -.8 .8 -.8 .8]);
view(az, el);
lighting gouraud;
lightangle(50, 60);
material([.4 1 0]);


% i: Simulate and Energy
subplot(12,24,nPan3+24*4+14); cla;
% Simulate Original ==> First Deformation
kT = ones(68,1);
LVal = sqrt((XsT(1,connT(:,1)) - XuT(1,connT(:,2)-9)).^2 +...
            (XsT(2,connT(:,1)) - XuT(2,connT(:,2)-9)).^2 +...
            (XsT(3,connT(:,1)) - XuT(3,connT(:,2)-9)).^2)';
[M, ES] = sim_motion3D_congrad(XsT, XuT, connT, LVal, 0.0004, 2513, [5 6], xF2(:,[1 2]), 0, kT);
XsT2 = M(:,1:size(XsT,2),end);
XuT2 = M(:,(1:size(XuT,2))+size(XsT,2),end);
% Simulate First Deformation ==> Second Deformation
[MS2, ES2] = sim_motion3D_congrad(XsT2, XuT2, connT, LVal, 0.0004, 2500, [1 2], xF1(:,[1 2]), 0, kT);
MA = cat(3,M,MS2);
EA = [ES ES2];
% Plot Networks
visualize_network(MA(:,1:9,1)+[-4 0 0]', MA(:,(1:16)+9,1)+[-4 0 0]', connT);
visualize_network(MA(:,1:9,size(M,3)), MA(:,(1:16)+9,size(M,3)), connT);
visualize_network(MA(:,1:9,end)+[4 0 0]', MA(:,(1:16)+9,end)+[4 0 0]', connT);
hold on;
quiver3(MA(1,[5,6],1)-4,MA(2,[5,6],1),MA(3,[5,6],1),...
        xF2(1,[1,2])-xS2(1,[1,2]),xF2(2,[1,2])-xS2(2,[1,2]),xF2(3,[1,2])-xS2(3,[1,2]),0,...
        'linewidth', 1, 'color', C_SA(1,:));
quiver3(MA(1,[1,2],size(M,3)),MA(2,[1,2],size(M,3)),MA(3,[1,2],size(M,3)),...
        xF1(1,[1,2])-xS1(1,[1,2]),xF1(2,[1,2])-xS1(2,[1,2]),xF1(3,[1,2])-xS1(3,[1,2]),0,...
        'linewidth', 1, 'color', C_SA(1,:));
hold off;
axis([-5.6 5.6 -1.5 2 -1 1]);
view(0,90);
set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[],'visible',0);
lighting gouraud;
lightangle(50, 60);
material([.4 1 0]);


% j: Energy
subplot(12,24,nPan3+24*8+14); cla;
hold on;
plot(EA, 'linewidth', 1);
axis([0 length(EA) 0 max(EA)]);
set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);
hold off;


% Size and Save Figure
fName = 'finite';
set(gcf, 'Renderer', 'opengl'); 
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 8.2 4.0]*1;
fig.PaperSize = [8.2 4.0]*1;
% print(['Figures\' fName], '-dpng','-r1500');

