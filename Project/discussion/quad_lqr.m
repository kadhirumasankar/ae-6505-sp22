clc; clear; close all;

% Script to simulate LQR control on an Iris quadcopter

% Physical parameters
mass = 1.545;
g = 9.81;
Ixx = 0.029125;
Iyy = 0.029125;
Izz = 0.055225;
I = [Ixx, Iyy, Izz]';

% Integration parameters
dt = .01;
t = 0:dt:20;
n = length(t);

% Creating a (number of states) x (number of timesteps) matrix to store
% the change in states over time
x = zeros(12, n);
% Randomly setting initial position to be regulated
x(1, 1) = -5+10*rand();
x(2, 1) = -5+10*rand();
x(3, 1) = -5+10*rand();

% Place holder for thrust and torque inputs over time
u = zeros(4, n);
% Solve for optimal K matrix
% The results don't change much between calculating K initially and at
% every timestep for a simple hover condition. It just makes the code
% slower
K = lqr_impl(mass, g, x(:, 1), I);

% Euler integration
for i = 1:(n-1)
    u(:, i) = -K*x(:, i) + [mass*g, 0, 0, 0]';
    x(:, i+1) = x(:, i) + dynamics(mass, g, x(:, i), u(:, i), I)*dt;
end

state_names = ["x", "y", "z", "vx", "vy", "vz", "roll", "pitch", "yaw", "p", "q", "r"];
for i=1:12
    subplot(4,3,i)
    hold on;
    plot(t, x(i,:))
    title(state_names(i))
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
    for i = 1:6
        Q(i, i) = .01;
    end
    for i = 10:12
        Q(i, i) = 2;
    end

    % R matrix gives weights to importance of regulating control effort in the
    % cost function
    R = eye(4);
    R(1) = 0.01;

    % matlab solves algebraic ricatti equation for optimal K matrix
    K = lqr(A, B, Q, R);
end