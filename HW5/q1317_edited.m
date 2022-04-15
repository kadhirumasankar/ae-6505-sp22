clc; clear all; close all;

Ra = 1.9; % Winding resistance
L = 0.003; % Winding inductance
lambda = 0.1; % Motor constant
J = 0.00018; % Moment of inertia
B = 0.001; % Coefficient of viscous friction

ControlNoise = 0.01; % std dev of uncertainty in control inputs
MeasNoise = 0.1; % standard deviation of measurement noise
R = [MeasNoise^2 0; 0 MeasNoise^2]; % Measurement noise covariance
xdotNoise = [0.1 0.1 0.5 0];
Q = [xdotNoise(1)^2 0 0 0; 0 xdotNoise(2)^2 0 0; 0 0 xdotNoise(3)^2 0; 0 0 0 xdotNoise(4)^2]; % Process noise covariance
P = 1*eye(4); % Initial state estimation covariance

dt = 0.001; % Integration step size
tf = 1.5; % Simulation length

x_obs = [0; 0; 0; 0]; % Initial state
x_est = x_obs; % State estimate
w = 2 * pi; % Control input frequency

dtPlot = 0.01; % How often to plot results
tPlot = -inf;

% Initialize arrays for plotting at the end of the program
xArray = [];
xhatArray = [];
trPArray = [];
tArray = [];
dt = 0.1;
xplot = [];

[t x_obs] = ode45(@obs, 0:dt:tf, x_obs);

% Begin simulation loop
dt = 0.1;
i = 1;
for t = 0 : dt : tf
    % Observed
    H = [1 0 0 0; 0 1 0 0];
    y_obs = H * x_obs(i,:)' + [MeasNoise*randn; MeasNoise*randn] / sqrt(dt);

    % Propagation of state
    ua0 = sin(w*t);
    ub0 = cos(w*t);
    xdot_est = [-Ra/L*x_est(1) + x_est(3)*lambda/L*sin(x_est(4)) + ua0/L;
                -Ra/L*x_est(2) - x_est(3)*lambda/L*cos(x_est(4)) + ub0/L;
                -3/2*lambda/J*x_est(1)*sin(x_est(4)) + 3/2*lambda/J*x_est(2)*cos(x_est(4)) - B/J*x_est(3);
                x_est(3)];
    x_est = x_est + xdot_est * dt;
    x_est(4) = mod(x_est(4), 2*pi);

    % Propagation of state covariance
    A = [-Ra/L 0 lambda/L*sin(x_est(4)) x_est(3)*lambda/L*cos(x_est(4));
        0 -Ra/L -lambda/L*cos(x_est(4)) x_est(3)*lambda/L*sin(x_est(4));
        -3/2*lambda/J*sin(x_est(4)) 3/2*lambda/J*cos(x_est(4)) -B/J -3/2*lambda/J*(x_est(1)*cos(x_est(4))+x_est(2)*sin(x_est(4)));
        0 0 1 0];
    Pdot = A * P + P * A' + Q;
    P = P + Pdot * dt;
    
    % Assembling y_computed
    y_comp = [x_est(1); x_est(2)];

    % Kalman Gain
    K = P*H'*inv(H*P*H'+R);

    % Measurement update
    x_est = x_est + K * (y_obs - y_comp);
    x_est(4) = mod(x_est(4), 2*pi);
    P = (eye(4)-K*H)*P*(eye(4)-K*H)' + K*R*K';
    figure(1)
    hold on
%     scatter3(x_obs(1), x_obs(2), x_obs(3), 20, 'red')
%     scatter3(x_est(1), x_est(2), x_est(3), 20, 'black')
%     scatter(t, x_est, 20, 'black')
    scatter(t, x_est(1), 20, 'red')
%     ylim([-5 5])
end


function xdot = obs(t, x_obs)
    w = 2 * pi; % Control input frequency
    Ra = 1.9; % Winding resistance
    L = 0.003; % Winding inductance
    lambda = 0.1; % Motor constant
    J = 0.00018; % Moment of inertia
    B = 0.001; % Coefficient of viscous friction
    ua0 = sin(w*t);
    ub0 = cos(w*t);
    xdot = [-Ra/L*x_obs(1) + x_obs(3)*lambda/L*sin(x_obs(4)) + ua0/L;
        -Ra/L*x_obs(2) - x_obs(3)*lambda/L*cos(x_obs(4)) + ub0/L;
        -3/2*lambda/J*x_obs(1)*sin(x_obs(4)) + 3/2*lambda/J*x_obs(2)*cos(x_obs(4)) - B/J*x_obs(3);
        x_obs(3)];
end