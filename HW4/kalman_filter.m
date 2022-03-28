%%
clc; clear all; close all;

% SET CONSTANTS
g = 9.81;

% Measurement noise covariance
% 5 ft and 0.1 deg in m and rad
R = [(1.524)^2 0 0;
     0 (0.1*pi/180)^2 0;
     0 0 (0.1*pi/180)^2];

data = importdata('homerun_data_HW4.txt');
data(1,:) = [];
tspan = data(:, 1);
rho = data(:, 2)/3.281; % in m
alpha = data(:, 3)*pi/180; % in rad
beta = data(:, 4)*pi/180; % in rad

% Assembling Y_obs using read in data
Y_obs = [];
for i = 1:length(rho)
    Y_obs = vertcat(Y_obs, [rho(i) alpha(i) beta(i)]);
end
m = size(Y_obs, 1);

% Visualizing data
figure(1)
[x y z] = sph2cart(alpha, beta, rho);
scatter3(x, y, z, 10, 'm', 'filled')

t0 = 0;
X0star = [0.4921 0.4921 2.0013 -26.2467 114.3051 65.9941]'/3.281; % in SI
xbar_0 = [0 0 0 0 0 0]';
Pbar_0 = diag([4 4 4 0.1 0.1 0.1])/3.281^2;
i = 1;

while i<=m
    if i == 1
        xhat_i_1 = X0star;
        P_i_1 = Pbar_0;
        t_i_1 = t0;
    else
        xhat_i_1 = xhat_i;
        P_i_1 = P_i;
        t_i_1 = tspan(i-1);
    end

    Phi = [1 0 0 tspan(i)-t_i_1 0 0;
           0 1 0 0 tspan(i)-t_i_1 0;
           0 0 1 0 0 tspan(i)-t_i_1;
           0 0 0 1 0 0;
           0 0 0 0 1 0;
           0 0 0 0 0 1];

    B = [0; 0; -1/2*(tspan(i)-t_i_1)^2; 0; 0; -(tspan(i)-t_i_1)];
    xbar_i = Phi * xhat_i_1 + B * g;
    Pbar_i = Phi * P_i_1 * Phi';
%     + Q to previous line to add filter from divering
    
    X = xbar_i(1);
    Y = xbar_i(2);
    Z = xbar_i(3);

    rho = sqrt(X^2 + Y^2 + Z^2);
    alpha = atan2(Y, X);
    beta = atan2(Z, sqrt(X^2 + Y^2));

    Y_computed = [rho; alpha; beta];
    y_i = Y_obs(i,:)' - Y_computed;
    H_tilda_i = [X/rho Y/rho Z/rho 0 0 0;
                 (-Y/X^2)/(1 + (Y/X)^2) (1/X)/(1+(Y/X)^2) 0 0 0 0;
                 ((-X*Z)/(X^2 + Y^2)^(3/2))/(1+(Z^2)/(X^2 + Y^2)) ((-Y*Z)/(X^2 + Y^2)^(3/2))/(1 + (Z^2)/(X^2 + Y^2)) ((1)/(X^2 + Y^2)^(1/2))/(1 + (Z^2)/(X^2 + Y^2)) 0 0 0];
    K_i = Pbar_i*H_tilda_i'*inv(H_tilda_i*Pbar_i*H_tilda_i' + R);
    
    xhat_i = xbar_i + K_i*(y_i);
    P_i = (eye(6) - K_i*H_tilda_i)*Pbar_i;
    X_calculated(:, i) = xhat_i;
    
    figure(1)
    hold on
    scatter3(X_calculated(1, i), X_calculated(2, i), X_calculated(3, i), 20, 'black')
    
    i = i+1;
end

xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
title('Trajectory of Baseball')

%% Batch approach
% Hyperparameters
n = 6; % number of params
k = 3; % number of obs per epoch
max_iter = 8;
g = 9.81; % m/s^2
dt = 0.1;
t_all = [0:dt:3.5];

% Measurement noise covariance
% 5 ft and 0.1 deg in m and rad
R = [(1.524)^2 0 0;
     0 (0.1*pi/180)^2 0;
     0 0 (0.1*pi/180)^2];

% Initial state
X_0 = [0.4921 0.4921 2.0013 -26.2467 114.3051 65.9941]'; % in ft and ft/s
X_0 = X_0/3.281; % converting to SI

% IF NO APRIORI COMMENT THIS OUT
P_bar = diag([4 4 4 0.1 0.1 0.1])/3.281^2;
x_bar = [0; 0; 0; 0; 0; 0];

% Initial measurements
data = importdata('../HW3/homerun_data.txt');
tspan = data(:, 1);
rho = data(:, 2)/3.281; % in m
alpha = data(:, 3)*pi/180; % in rad
beta = data(:, 4)*pi/180; % in rad

% Assembling Y_obs using read in data
Y_obs = [];
for i = 1:length(rho)
    Y_obs = vertcat(Y_obs, [rho(i) alpha(i) beta(i)]);
end
m = size(Y_obs, 1);

% Visualizing data
figure(2)
[x y z] = sph2cart(alpha, beta, rho)
scatter3(x, y, z, 10, 'm', 'filled')

% Run iterations
for i = 1:max_iter
    HtH = zeros(n, n);
    Hty = zeros(n, 1);
    
    % Loop over all observations
    % Compute the nominal trajectory and computed observations
    % Compute Normal equations
    for j=1:m
        Phi = [1 0 0 tspan(j) 0 0;
               0 1 0 0 tspan(j) 0;
               0 0 1 0 0 tspan(j);
               0 0 0 1 0 0;
               0 0 0 0 1 0;
               0 0 0 0 0 1];
        
        % Compute nominal trajectory
        B = [0; 0; -1/2*tspan(j)^2; 0; 0; -tspan(j)];
        X_nom = Phi * X_0 + B * g;
        
        X = X_nom(1);
        Y = X_nom(2);
        Z = X_nom(3);
        
        % Compute partials
        rho = sqrt(X^2 + Y^2 + Z^2);
        
        H_tilda = [X/rho Y/rho Z/rho 0 0 0;
                   (-Y/X^2)/(1 + (Y/X)^2) (1/X)/(1+(Y/X)^2) 0 0 0 0;
                   ((-X*Z)/(X^2 + Y^2)^(3/2))/(1+(Z^2)/(X^2 + Y^2)) ((-Y*Z)/(X^2 + Y^2)^(3/2))/(1 + (Z^2)/(X^2 + Y^2)) ((1)/(X^2 + Y^2)^(1/2))/(1 + (Z^2)/(X^2 + Y^2)) 0 0 0];
        
        % Generate observation deviation (i.e. observed - computed)
        rho = sqrt(X^2 + Y^2 + Z^2);
        alpha = atan2(Y, X);
        beta = atan2(Z, sqrt(X^2 + Y^2));
        
        Y_computed = [rho; alpha; beta];
        y = Y_obs(j,:)' - Y_computed;
        
        % Use STM to map partials back to initial coords
        H = H_tilda * Phi;
        
        % Compute rank-k update of normal equations for this observation
        HtH = HtH + H'*inv(R)*H;
        Hty = Hty + H'*inv(R)*y;
    end
    
    % Compute state deviation
    P = inv(HtH);
    
    % Uncomment line 1 if no apriori, uncomment line 2 if given apriori info
%     x_hat = P * Hty;
    x_hat = inv(HtH + inv(P_bar))*(Hty + inv(P_bar)*x_bar);
    
    % Add state deviation to get final estimate
    X_new = X_0 + x_hat
    
    % Use the new estimates as starting point for next iter
    X_0 = X_new;
    P_bar = inv(HtH + inv(P_bar));
    
    %Plot estimated trajectory for current iteration
    figure(2);
    hold on;
    if i==1 || i == 3
        for k=1:size(t_all, 2)
            X0 = X_0(1);
            Y0 = X_0(2);
            Z0 = X_0(3);
            Vx0 = X_0(4);
            Vy0 = X_0(5);
            Vz0 = X_0(6);
            X = X0 + Vx0 * t_all(k);
            Y = Y0 + Vy0 * t_all(k);
            Z = Z0 + Vz0 * t_all(k) - (1/2)*g*t_all(k)^2;
            X_KF(i, k) = X;
            Y_KF(i, k) = Y;
            Z_KF(i, k) = Z;
            scatter3(X, Y, Z, 20)
        end
    end
end

xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
title('Trajectory of Baseball')

%%
close all;

figure()
scatter3(X_calculated(1,:), X_calculated(2,:), X_calculated(3,:), 20, 'black')
hold on;
scatter3(X_KF(1, :), Y_KF(1, :), Z_KF(1, :), 10, '*')
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
title('Trajectory of Baseball')
legend('Kalman Filtering', 'Batch 1st Iteration', 'Location', 'best')
hold off;

figure()
scatter3(X_calculated(1,:), X_calculated(2,:), X_calculated(3,:), 20, 'black')
hold on;
scatter3(X_KF(3, :), Y_KF(3, :), Z_KF(3, :),  10, '*')
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
title('Trajectory of Baseball')
legend('Kalman Filtering', 'Batch 3rd Iteration', 'Location', 'best')
hold off;