%% Prepare Space
clear; clc;


%% Figure 2d
% Positions and Motions
Xs = [[-1 1 0];...
      [0 0 1]];
Us = [[-1 1 -1];...
      [-1 -1 0]]/2;
figure(1); clf;
visualize_conic(Xs, Us, [-1 1; -1 1]*2, [100 100], 40, .5, .5);
axis(1.8*[-1 1 -1 1]);


%% Figure 2e
% Positions and Motions
Xs = [[-1 1 0 0 0];...
	  [-1 -1 1 0 0];...
	  [0 0 0 1 -1]];
Us = [[1 1 -1 0 0 ];...
	  [1 -1 0.1 0 0 ];...
	  [-1 -1 -1 1 1]];
figure(1); clf;
visualize_conic(Xs, Us, [-1 1; -1 1; -1 1]*4, [100 100 100], 60, .2, .2);
axis(1.5*[-1 1 -1 1 -1 1]);
view(-160, 20);


%% Figure 2f
% Positions and Motions
Xs = [[-1 1];...
      [0 0]];
Us = [[-1 1];...
      [0 0]];
figure(1); clf;
visualize_conic(Xs, Us, [-1 1; -1 1]*2, [100 100], 9, .2, .2);
axis(1.8*[-1 1 -1 1]);


%% Figure 2g
% Positions and Motions
Xs = [[-1 1 0 0];...
      [-1 -1 1 0];...
      [0 0 0 1]];
Us = [[-1 1 0 0]/2;...
      [-1 -1 1 0]/4;...
      [-1 -1 -1 -1.6]];
figure(1); clf;
visualize_conic(Xs, Us, [-1 1; -1 1; -1 1]*1.5, [100 100 100], 120, .2, .2);
axis(1.5*[-1 1 -1 1 -1 1]);
view(-160, 20);


%% Figure 2h
% Positions and Motions
Xs = [[-1 -1 1 1];...
      [-1 1 1 -1]];
Us = [[-1 -1 1 1]/2;...
      [-1 1 1 -1]/2];
figure(1); clf;
visualize_conic(Xs, Us, [-1 1; -1 1]*2, [100 100], 21, .5, .5);
axis(1.8*[-1 1 -1 1]);


%% Figure 2i
% Positions and Motions
Xs = [[-1 -1 -1 -1 1 1 1 1];...
      [-1 -1 1 1 -1 -1 1 1];...
	  [-1 1 -1 1 -1 1 -1 1]];
% Motions: 1D Manifold
Us = [[-1 -1 -1 -1 1 1 1 1];...
	  [-1 -1 1 1 -1 -1 1 1];...
	  [-1 1 -1 1 -1 1 -1 1]];
figure(1); clf;
visualize_conic(Xs, Us, [-1 1; -1 1; -1 1]*1.8, [100 100 100], 200, .2, .2);
axis(1.8*[-1 1 -1 1 -1 1]);
view(-160, 20);


%% Figure 3a-c
% Positions and Motions
Xs = [[-1 -1 1 1];...
	  [-1 1 1 -1]];
Us = [[-1 -1 1 1];...
	  [-1 1 1 -1]];
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

figure(1); clf;
% Visualize Conics
subplot(1,3,1);
visualize_conic(Xs, Us, [-1 1; -1 1]*2, [100 100], 21, .5, 0);
axis(1.8*[-1 1 -1 1]);
% Find positions of unspecified nodes, and corresponding motions
subplot(1,3,2);
[Xu1 fV1] = construct_network(Xs, Us, x01, conn1, 0);
[Us1, Uu1, err1] = construct_motion(Xs, Us, Xu1, conn1, .5, .5);
axis(1.8*[-1 1 -1 1]);
subplot(1,3,3);
[Xu2 fV2] = construct_network(Xs, Us, x02, conn2, 0);
[Us2, Uu2, err2] = construct_motion(Xs, Us, Xu2, conn2, .5, .5);
axis(1.8*[-1 1 -1 1]);
disp('Conic evaluation at Xu')
disp([fV1 fV2]);
disp('Motion reconstruction error');
disp([err1, err2]);


%% Figure 3d-f
% Positions and Motions
Xs = [[-1 -1 -1 -1 1 1 1 1];...
	  [-1 -1 1 1 -1 -1 1 1];...
	  [-1 1 -1 1 -1 1 -1 1]];
Us = [[-1 -1 -1 -1 1 1 1 1];...
	  [-1 -1 1 1 -1 -1 1 1];...
	  [-1 1 -1 1 -1 1 -1 1]];
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
	   [0         0        0        0       -sqrt(3)  sqrt(3)]];
x02 = [[-sqrt(3)  sqrt(3)  0        0        0      ];...
	   [0         0       -sqrt(3)  sqrt(3)  0      ];...
	   [0         0        0        0       -sqrt(3)]];

figure(1); clf;
% Visualize conic solution space
subplot(1,3,1);
visualize_conic(Xs, Us, [-1 1; -1 1; -1 1]*1.8, [100 100 100], 200, .4, 0);
axis(1.8*[-1 1 -1 1 -1 1]);
view(-160, 20);
% Find positions of unspecified nodes, and corresponding motions
subplot(1,3,2);
[Xu1 fV1] = construct_network(Xs, Us, x01, conn1, 0);
[Us1, Uu1, err1] = construct_motion(Xs, Us, Xu1, conn1, .5, .5);
axis(1.8*[-1 1 -1 1 -1 1]);
view(-160, 20);
subplot(1,3,3);
[Xu2 fV2] = construct_network(Xs, Us, x02, conn2, 0);
[Us2, Uu2, err2] = construct_motion(Xs, Us, Xu2, conn2, .5, .5);
axis(1.8*[-1 1 -1 1 -1 1]);
view(-160, 20);
disp('Conic evaluation at Xu')
disp([fV1 fV2]);
disp('Motion reconstruction error');
disp([err1, err2]);


%% Figure 4 b-c
% Positions and Motions
Xs = [[-1 0 1];...
	  [0 1 0]];
% Two separate motions
Us = zeros(2,3,2);
Us(:,:,1) = [[-1 0 1];...
             [0 1 0]];
Us(:,:,2) = [[-1 1 1];...
             [.5 -1 .5]];
% Initial guess for unspecified nodes
x0 = [-.7 .7]';
% Connectivity
conn = [1 1; 2 1; 3 1];

% Plot
figure(1); clf;
subplot(1,2,1);
visualize_conic(Xs, Us, [-1 1; -1 1]*2, [100 100], 21, .5, 0);
axis(1.8*[-1 1 -1 1]);
subplot(1,2,2);
[Xu1, fV] = construct_network(Xs, Us, x0, conn, 0);
[Us1, Uu1, err] = construct_motion(Xs, Us, Xu1, conn, .5, .5);
axis(1.8*[-1 1 -1 1]);
disp('Conic evaluation at Xu')
disp([fV]);
disp('Motion reconstruction error');
disp([err]);


%% Figure 4 e-f
% Positions and Motions
Xs = [[-1 0 1 0];...
      [-1 1 -1 0];...
	  [0 0 0 1]];
% Two separate motions
Us = zeros(3,4,2);
Us(:,:,1) = [[-1.0  0.0  1.0  0.0];...
             [-1.0  1.0 -1.0  0.0];...
             [0.0  0.0  0.0  1.0]];
Us(:,:,2) = [[0.7 -0.4  1.2 -0.2];...
             [1.0 -0.2 -0.5  0.9];...
             [-1.0  0.0  0.0  0.0]];
% Initial guess for unspecified nodes
x0 = [-.7 -.5 .8;...
      1.2 -.4 -.75;...
      -.4 -1.5 -.3;...
      .2 .6 .6]';
% Connectivity
conn = [1 1; 2 1; 3 1; 4 1; 1 2; 2 2; 3 2; 4 2;...
        1 3; 2 3; 3 3; 4 3; 1 4; 2 4; 3 4; 4 4];
    
% Plot
figure(1); clf;
subplot(1,2,1);
visualize_conic(Xs, Us, [-1 1; -1 1; -1 1]*1.8, [100 100 100], 200, .4, 0);
axis(1.8*[-1 1 -1 1 -1 1]);
view(-175, 10);
subplot(1,2,2);
[Xu1, fV] = construct_network(Xs, Us, x0, conn, 0);
[Us1, Uu1, err] = construct_motion(Xs, Us, Xu1, conn, .5, .5);
axis(1.8*[-1 1 -1 1 -1 1]);
view(-175, 10);
disp('Conic evaluation at Xu')
disp([fV]);
disp('Motion reconstruction error');
disp([err]);


%% Figure 5b-d
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

% Simulate Motion
figure(1); clf;
subplot(1,2,1);
[Xu, fV] = construct_network(Xs, Us, x0, conn, 0);
[Us, Uu, err] = construct_motion(Xs, -Xs, Xu, conn, .25, .25);
axis(1.8*[-1 1 -1 1]);
subplot(1,2,2);
[XMot, fC] = sim_motion(Xs, Xu, conn, .01, 100, -[Xs Xu]);
axis(1.8*[-1 1 -1 1]);
disp('Conic evaluation at Xu')
disp([fV]);
disp('Motion reconstruction error');
disp([err]);
disp('Simulation total error');
disp(sum(abs(fC)));


%% Figure 5e-g
% Positions and Motions
Xs = [[2 1 3 3.5 3.5];...
      [0 0 0 -1 1];...
	  [1 2 3 2.5 2.5]];
Us = [-[-1 1 -1 0 0]*.8;...
      -[0 0 0 -1 1]*.8;...
      -[1 -1 0 0 0]*.8];
% Initial guess for unspecified nodes
x0 = [3.092 0.808 2.183;...
       2.758 -0.930 1.516;...
       2.758 0.930 1.516;...
       2 -.77 2.061;...
       2 .77 2.061]';
% Connectivity
conn = [1 1; 2 1; 3 1; 4 1; 5 1;...
        1 2; 2 2; 3 2; 4 2; 5 2;...
        1 3; 2 3; 3 3; 4 3; 5 3;...
        1 4; 2 4; 3 4; 4 4;...
        1 5; 2 5; 3 5; 5 5];

% Simulate Motion
figure(1); clf;
subplot(1,2,1);
[Xu fV] = construct_network(Xs, Us, x0, conn, 0);
[Usp, Uup, err] = construct_motion(Xs, -Us, Xu, conn, .5, .5);
axis([0 4 -1.5 1.5 0.5 3.5]);
view(-15, 10);
subplot(1,2,2);
[XMot, fC] = sim_motion3D(Xs, Xu, conn, .01, 200, [Xs Xu]);
axis([0 4 -1.5 1.5 0.5 3.5]);
view(-15, 10);
disp('Conic evaluation at Xu')
disp([fV]);
disp('Motion reconstruction error');
disp([err]);
disp('Simulation total error');
disp(sum(abs(fC)));