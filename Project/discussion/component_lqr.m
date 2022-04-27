%%
clc; clear; close all;

%%
% Script to simulate LQR control on an Iris quadcopter

% Physical parameters
mass = 5.27996526;
g = 9.81;
Ixx = 0.171308;
Iyy = 0.171958;
Izz = 0.327894;
I = [Ixx, Iyy, Izz]';

% Integration parameters
dt = 0.001;
t = 0:dt:20;
n = length(t);

%% Getting initial state from a previous run
data = readmatrix('C:\Users\kumasan\Downloads\2021-06-07-10-26-42.csv');

x = zeros(12, n);
x(1, 1) = data(2,1);
x(2, 1) = data(3,1);
x(3, 1) = data(4,1);
x(4, 1) = data(5,1);
x(5, 1) = data(6,1);
x(6, 1) = data(7,1);
x(7, 1) = data(8,1);
x(8, 1) = data(9,1);
x(9, 1) = data(10,1);
x(10, 1) = data(11,1);
x(11, 1) = data(12,1);
x(12, 1) = data(13,1);

%% Desired state
x_des = zeros(12, 1);
x_des(3) = 3;

% Place holder for thrust and torque inputs over time
u = zeros(4, n);
% Solve for optimal K matrix
% The results don't change much between calculating K initially and at
% every timestep for a simple hover condition. It just makes the code
% slower
K = lqr_impl(mass, g, x(:, 1), I);

%% Splitting u into components
u_rate = 50; % Control is run at this rate
latency = 0; % Latency in ms
% This is subtracted from i in the loop below to simulate latency
latency_idx = round(latency/1000/dt);

for i = 1:(n-1)
    if mod(i-1, floor((u_rate*dt)^-1)) == 0
        for j = 1:12
            current_u(:, j) = -K(:, j)*(x(j, max([1, i-latency_idx]))-x_des(j,1));
        end
    end
    u(:, i) = sum(current_u, 2) + [mass*g, 0, 0, 0]';
    x(:, i+1) = x(:, i) + dynamics(mass, g, x(:, i), u(:, i), I)*dt;
end

state_names = ["x", "y", "z", "vx", "vy", "vz", "roll", "pitch", "yaw", "p", "q", "r"];
figure(1)
for i=1:12
    subplot(4,3,i)
    hold on;
    plot(t, x(i,:))
    title(state_names(i))
    hold off;
end

u_names = ["Thrust", "Roll torque", "Pitch torque", "Yaw torque"];
figure(2)
for i=1:4
    subplot(4,1,i)
    hold on;
    plot(t, u(i, :))
    title(u_names(i))
    hold off;
end

%% Finding prop speeds
Ct = 5.84e-6;
% Found using (Ct * moment constant)/D = 5.84e-6*0.06/0.2413
Cq = 1.4521e-6;
Lrf = 0.22;
Lrb = 0.2;
Lp = 0.13;

% M from docs
% M = [Ct, Ct, Ct, Ct;
%      -Ct*Lrf, Ct*Lrb, Ct*Lrf, -Ct*Lrb;
%      Ct*Lp, -Ct*Lp, Ct*Lp, -Ct*Lp;
%      Cq, Cq, -Cq, -Cq];
 
% M after flipping pitch signs
% M = [Ct, Ct, Ct, Ct;
%      -Ct*Lrf, Ct*Lrb, Ct*Lrf, -Ct*Lrb;
%      -Ct*Lp, Ct*Lp, -Ct*Lp, Ct*Lp;
%      Cq, Cq, -Cq, -Cq];
 
% M after flipping yaw signs
M = [Ct, Ct, Ct, Ct;
     -Ct*Lrf, Ct*Lrb, Ct*Lrf, -Ct*Lrb;
     -Ct*Lp, Ct*Lp, -Ct*Lp, Ct*Lp;
     -Cq, -Cq, Cq, Cq];

inv_M = inv(M);

fprintf("inv_M << %g, %g, %g, %g,\n", inv_M(1,1), inv_M(1,2), inv_M(1,3), inv_M(1,4));
fprintf("\t\t %g, %g, %g, %g,\n", inv_M(2,1), inv_M(2,2), inv_M(2,3), inv_M(2,4));
fprintf("\t\t %g, %g, %g, %g,\n", inv_M(3,1), inv_M(3,2), inv_M(3,3), inv_M(3,4));
fprintf("\t\t %g, %g, %g, %g;\n", inv_M(4,1), inv_M(4,2), inv_M(4,3), inv_M(4,4));
 
for i = 1:n
    prop_speeds(:, i) = (inv_M * u(:, i));
end

prop_speeds = max(prop_speeds, 100^2);
prop_speeds = min(prop_speeds, 1100^2);
prop_speeds = sqrt(prop_speeds);

figure(3)
subplot(2,2,1)
plot(t, prop_speeds(3, :), 'LineWidth', 1.5)
hold on
yline(1100, '--g')
yline(100, '--r')
hold off
subplot(2,2,2)
plot(t, prop_speeds(1, :), 'LineWidth', 1.5)
hold on
yline(1100, '--g')
yline(100, '--r')
hold off
subplot(2,2,3)
plot(t, prop_speeds(2, :), 'LineWidth', 1.5)
hold on
yline(1100, '--g')
yline(100, '--r')
hold off
subplot(2,2,4)
plot(t, prop_speeds(4, :), 'LineWidth', 1.5)
hold on
yline(1100, '--g')
yline(100, '--r')
hold off

%% Remultiplying prop speeds with M to get u
u2 = M * prop_speeds.^2;

figure(4)
sgtitle('u = M * w^2')
for i=1:4
    subplot(4,1,i)
    hold on;
    plot(t, u(i, :), 'LineWidth', 1.5)
    plot(t, u2(i, :))
    title(u_names(i))
    hold off;
end

%% Functions
function x_dot = dynamics(mass, g, states, u, I)
    q4 = states(4);
    q5 = states(5);
    q6 = states(6);
    q7 = states(7);
    q8 = states(8);
    q9 = states(9);
    q10 = states(10);
    q11 = states(11);
    q12 = states(12);
    u3 = u(1);
    u4 = u(2);
    u5 = u(3);
    u6 = u(4);
    m = mass;
    Ixx = I(1);
    Iyy = I(2);
    Izz = I(3);
    x_dot =[
        q4;
        q5;
        q6;
        (u3*(sin(q7)*sin(q9) + cos(q7)*cos(q9)*sin(q8)))/m;
        -(u3*(cos(q9)*sin(q7) - cos(q7)*sin(q8)*sin(q9)))/m;
        -(g*m - u3*cos(q7)*cos(q8))/m;
        q10 + q12*cos(q7)*tan(q8) + q11*sin(q7)*tan(q8);
        q11*cos(q7) - q12*sin(q7);
        (q12*cos(q7))/cos(q8) + (q11*sin(q7))/cos(q8);
        (u4 + Iyy*q11*q12 - Izz*q11*q12)/Ixx;
        (u5 - Ixx*q10*q12 + Izz*q10*q12)/Iyy;
        (u6 + Ixx*q10*q11 - Iyy*q10*q11)/Izz;
    ];
end

function K = lqr_impl(mass, g, states, I)
    % State dependent dynamics
    u3 = mass * g;
    q7 = states(7);
    q8 = states(8);
    q9 = states(9);
    q10 = states(10);
    q11 = states(11);
    q12 = states(12);
    Ixx = I(1);
    Iyy = I(2);
    Izz = I(3);
    m = mass;
    A = [
        0, 0, 0, 1, 0, 0,                                                   0,                                                                 0,                                                  0,                        0,                       0,                        0;
        0, 0, 0, 0, 1, 0,                                                   0,                                                                 0,                                                  0,                        0,                       0,                        0;
        0, 0, 0, 0, 0, 1,                                                   0,                                                                 0,                                                  0,                        0,                       0,                        0;
        0, 0, 0, 0, 0, 0,  (u3*(cos(q7)*sin(q9) - cos(q9)*sin(q7)*sin(q8)))/m,                                    (u3*cos(q7)*cos(q8)*cos(q9))/m, (u3*(cos(q9)*sin(q7) - cos(q7)*sin(q8)*sin(q9)))/m,                        0,                       0,                        0;
        0, 0, 0, 0, 0, 0, -(u3*(cos(q7)*cos(q9) + sin(q7)*sin(q8)*sin(q9)))/m,                                    (u3*cos(q7)*cos(q8)*sin(q9))/m, (u3*(sin(q7)*sin(q9) + cos(q7)*cos(q9)*sin(q8)))/m,                        0,                       0,                        0;
        0, 0, 0, 0, 0, 0,                             -(u3*cos(q8)*sin(q7))/m,                                           -(u3*cos(q7)*sin(q8))/m,                                                  0,                        0,                       0,                        0;
        0, 0, 0, 0, 0, 0,           q11*cos(q7)*tan(q8) - q12*sin(q7)*tan(q8),         q12*cos(q7)*(tan(q8)^2 + 1) + q11*sin(q7)*(tan(q8)^2 + 1),                                                  0,                        1,         sin(q7)*tan(q8),          cos(q7)*tan(q8);
        0, 0, 0, 0, 0, 0,                         - q12*cos(q7) - q11*sin(q7),                                                                 0,                                                  0,                        0,                 cos(q7),                 -sin(q7);
        0, 0, 0, 0, 0, 0,       (q11*cos(q7))/cos(q8) - (q12*sin(q7))/cos(q8), (q12*cos(q7)*sin(q8))/cos(q8)^2 + (q11*sin(q7)*sin(q8))/cos(q8)^2,                                                  0,                        0,         sin(q7)/cos(q8),          cos(q7)/cos(q8);
        0, 0, 0, 0, 0, 0,                                                   0,                                                                 0,                                                  0,                        0, (Iyy*q12 - Izz*q12)/Ixx,  (Iyy*q11 - Izz*q11)/Ixx;
        0, 0, 0, 0, 0, 0,                                                   0,                                                                 0,                                                  0, -(Ixx*q12 - Izz*q12)/Iyy,                       0, -(Ixx*q10 - Izz*q10)/Iyy;
        0, 0, 0, 0, 0, 0,                                                   0,                                                                 0,                                                  0,  (Ixx*q11 - Iyy*q11)/Izz, (Ixx*q10 - Iyy*q10)/Izz,                        0;
    ];
    % Actuator dependent dynamics
    B = [
        0, 0, 0, 0;
        0, 0, 0, 0;
        0, 0, 0, 0;
        (sin(q7) * sin(q9) + cos(q7) * cos(q9) * sin(q8)) / m, 0, 0, 0;
        -(cos(q9) * sin(q7) - cos(q7) * sin(q8) * sin(q9)) / m, 0, 0, 0;
        (cos(q7) * cos(q8)) / m, 0, 0, 0;
        0, 0, 0, 0;
        0, 0, 0, 0;
        0, 0, 0, 0;
        0, 1 / Ixx, 0, 0;
        0, 0, 1 / Iyy, 0;
        0, 0, 0, 1 / Izz;
    ];

    % Q matrix gives weights to importance of regulating each state in the 
    % cost function
    Q = eye(12);
    for i = 1:3
        Q(i, i) = 100;
    end
    for i = 4:6
        Q(i, i) = 100;
    end
    Q(6, 6) = 200;
    for i=7:9
        Q(i, i) = 0.1;
    end
    for i=10:12
        Q(i, i) = 10;
    end

    % R matrix gives weights to importance of regulating control effort in the
    % cost function
    R = eye(4);
    R(1,1) = 1;
    for i=2:4
        R(i, i) = 1000;
    end

    % matlab solves algebraic ricatti equation for optimal K matrix
    K = lqr(A, B, Q, R);
    fprintf("K << %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g,\n", K(1,1), K(1,2), K(1,3), K(1,4), K(1,5), K(1,6), K(1,7), K(1,8), K(1,9), K(1,10), K(1,11), K(1,12));
    fprintf("\t %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g,\n", K(2,1), K(2,2), K(2,3), K(2,4), K(2,5), K(2,6), K(2,7), K(2,8), K(2,9), K(2,10), K(2,11), K(2,12));
    fprintf("\t %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g,\n", K(3,1), K(3,2), K(3,3), K(3,4), K(3,5), K(3,6), K(3,7), K(3,8), K(3,9), K(3,10), K(3,11), K(3,12));
    fprintf("\t %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g;\n", K(4,1), K(4,2), K(4,3), K(4,4), K(4,5), K(4,6), K(4,7), K(4,8), K(4,9), K(4,10), K(4,11), K(4,12));
end