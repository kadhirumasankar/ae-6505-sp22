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
% data(1,:) = []; % COME BACK AND UNCOMMENT THIS LATER
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
[x y z] = sph2cart(alpha, beta, rho)
scatter3(x, y, z, 10, 'm', 'filled')

t0 = 0;
X0star = [0.4921 0.4921 2.0013 -26.2467 114.3051 65.9941]'/3.281; % in SI
xbar_0 = [0 0 0 0 0 0]';
Pbar_0 = [1.2192^2 0 0 0 0 0;
          0 1.2192^2 0 0 0 0;
          0 0 1.2192^2 0 0 0;
          0 0 0 0.03048^2 0 0;
          0 0 0 0 0.03048^2 0;
          0 0 0 0 0 0.03048^2];
i = 1;
Xstar_i_1 = X0star;
% xhat_i_1(1) = xbar_0;
% P_i_1(1) = Pbar_0;

while i<=m
    if i == 1
        xhat_i_1 = xbar_0;
        P_i_1 = Pbar_0;
        t_i_1 = t0;
    else
        xhat_i_1 = xhat_i;
        P_i_1 = P_i;
        t_i_1 = tspan(i-1);
    end
    % MIGHT NEED TO CHANGE THIS BC IT'S NO LONGER PHI(T, T0) IT'S PHI(T_I,
    % T_I-1)
    Phi = [1 0 0 tspan(i)-t_i_1 0 0;
           0 1 0 0 tspan(i)-t_i_1 0;
           0 0 1 0 0 tspan(i)-t_i_1;
           0 0 0 1 0 0;
           0 0 0 0 1 0;
           0 0 0 0 0 1];

    B = [0; 0; -1/2*tspan(i)^2; 0; 0; -tspan(i)];
    X_nom = Phi * X0star + B * g;
    
    xbar_i = Phi * xhat_i_1;
    Pbar_i = Phi * P_i_1 * Phi';
    
    X = X_nom(1);
    Y = X_nom(2);
    Z = X_nom(3);

    rho = sqrt(X^2 + Y^2 + Z^2);
    alpha = atan2(Y, X);
    beta = atan2(Z, sqrt(X^2 + Y^2));

    Y_computed = [rho; alpha; beta];
    y_i = Y_obs(i,:)' - Y_computed;
    H_tilda_i = [X/rho Y/rho Z/rho 0 0 0;
                 (-Y/X^2)/(1 + (Y/X)^2) (1/X)/(1+(Y/X)^2) 0 0 0 0;
                 ((-X*Z)/(X^2 + Y^2)^(3/2))/(1+(Z^2)/(X^2 + Y^2)) ((-Y*Z)/(X^2 + Y^2)^(3/2))/(1 + (Z^2)/(X^2 + Y^2)) ((1)/(X^2 + Y^2)^(1/2))/(1 + (Z^2)/(X^2 + Y^2)) 0 0 0];
    K_i = Pbar_i*H_tilda_i'*inv(H_tilda_i*Pbar_i*H_tilda_i' + R);
    
    xhat_i = xbar_i + K_i*(y_i - H_tilda_i*xbar_i);
    P_i = (eye(6) - K_i*H_tilda_i)*Pbar_i;
    i = i+1;
    X_calculated(:, i) = X_nom + xhat_i;
end

figure(1)
hold on
scatter3(X_calculated(1, :), X_calculated(2, :), X_calculated(3, :), 20)