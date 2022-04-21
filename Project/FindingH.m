clc; clear all; close all;
syms drone_x drone_y drone_z vx vy vz roll pitch yaw p q r X1 X2 X3 X4 Y1 Y2 Y3 Y4 Z1 Z2 Z3 Z4
G = [sqrt((drone_x-X1)^2 + (drone_y-Y1)^2 + (drone_z-Z1)^2);
     sqrt((drone_x-X2)^2 + (drone_y-Y2)^2 + (drone_z-Z2)^2);
     sqrt((drone_x-X3)^2 + (drone_y-Y3)^2 + (drone_z-Z3)^2);
     sqrt((drone_x-X4)^2 + (drone_y-Y4)^2 + (drone_z-Z4)^2)];
H = jacobian(G,[drone_x; drone_y; drone_z; vx; vy; vz; roll; pitch; yaw; p; q; r]);
simplify(H)