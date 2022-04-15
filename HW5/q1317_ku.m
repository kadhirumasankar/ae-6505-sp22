clc; clear all; close all;

% Optimal State Estimation, by Dan Simon
% Corrected June 18, 2012, and November 10, 2013, thanks to Jean-Michel Papy
%
% Continuous time extended Kalman filter simulation for two-phase step motor.
% Estimate the stator currents, and the rotor position and velocity, on the
% basis of noisy measurements of the stator currents.

Ra = 1.9; % Winding resistance
L = 0.003; % Winding inductance
lambda = 0.1; % Motor constant
J = 0.00018; % Moment of inertia
B = 0.001; % Coefficient of viscous friction

ControlNoise = 0.01; % std dev of uncertainty in control inputs
MeasNoise = 0.1; % standard deviation of measurement noise
R = [MeasNoise^2 0; 0 MeasNoise^2]; % Measurement noise covariance
xdotNoise = [ControlNoise/L ControlNoise/L 0.5 0];
Q = [xdotNoise(1)^2 0 0 0; 0 xdotNoise(2)^2 0 0; 0 0 xdotNoise(3)^2 0; 0 0 0 xdotNoise(4)^2]; % Process noise covariance
P = 1*eye(4); % Initial state estimation covariance

dt = 0.001; % Integration step size
tf = 1.5; % Simulation length

x = [0; 0; 0; 0]; % Initial state
xhat = x; % State estimate
w = 2 * pi; % Control input frequency

dtPlot = 0.01; % How often to plot results
tPlot = -inf;

% Initialize arrays for plotting at the end of the program
xArray = [];
xhatArray = [];
trPArray = [];
tArray = [];

% Begin simulation loop
for t = 0 : dt : tf
    % Nonlinear simulation
    ua0 = sin(w*t);
    ub0 = cos(w*t);
    xdot = [-Ra/L*x(1) + x(3)*lambda/L*sin(x(4)) + ua0/L;
        -Ra/L*x(2) - x(3)*lambda/L*cos(x(4)) + ub0/L;
        -3/2*lambda/J*x(1)*sin(x(4)) + 3/2*lambda/J*x(2)*cos(x(4)) - B/J*x(3);
        x(3)];
    x = x + xdot * dt + [xdotNoise(1)*randn; xdotNoise(2)*randn; xdotNoise(3)*randn; xdotNoise(4)*randn] * sqrt(dt);
    x(4) = mod(x(4), 2*pi);
    % Kalman filter    
    H = [1 0 0 0; 0 1 0 0];
    y_obs = H * x + [MeasNoise*randn; MeasNoise*randn] / sqrt(dt);
   
    F = [-Ra/L 0 lambda/L*sin(xhat(4)) xhat(3)*lambda/L*cos(xhat(4));
        0 -Ra/L -lambda/L*cos(xhat(4)) xhat(3)*lambda/L*sin(xhat(4));
        -3/2*lambda/J*sin(xhat(4)) 3/2*lambda/J*cos(xhat(4)) -B/J -3/2*lambda/J*(xhat(1)*cos(xhat(4))+xhat(2)*sin(xhat(4)));
        0 0 1 0];
    xhatdot = [-Ra/L*xhat(1) + xhat(3)*lambda/L*sin(xhat(4)) + ua0/L;
        -Ra/L*xhat(2) - xhat(3)*lambda/L*cos(xhat(4)) + ub0/L;
        -3/2*lambda/J*xhat(1)*sin(xhat(4)) + 3/2*lambda/J*xhat(2)*cos(xhat(4)) - B/J*xhat(3);
        xhat(3)];
    Pdot = F * P + P * F' + Q;
    xhat = xhat + xhatdot*dt;
    xhat(4) = mod(xhat(4), 2*pi);
    P = P + Pdot*dt;
    y_comp = H * xhat + [MeasNoise*randn; MeasNoise*randn] / sqrt(dt);
    
    K = P * H' * inv(H*P*H' + R);
    xhat = xhat + K*(y_obs - y_comp);
    xhat(4) = mod(xhat(4), 2*pi);
    P = (eye(4)-K*H)*P*(eye(4)-K*H)' + K*R*K';
 
    figure(1)
    hold on
    scatter(t, y_comp(1))
end
