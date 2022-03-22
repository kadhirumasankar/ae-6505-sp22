clc; clear all; close all;

% Given info

% Measurement noise covariance
R0 = [5^2 0 0;
      0 (0.1*pi/180)^2 0;
      0 0 (0.1*pi/180)^2];
% Initial state
X0 = [0.4921; 0.4921; 2.0013; -26.2467; 114.3051; 65.9941];
P0_bar = [inf 0 0 0 0 0;
          0 inf 0 0 0 0;
          0 0 inf 0 0 0;
          0 0 0 inf 0 0;
          0 0 0 0 inf 0;
          0 0 0 0 0 inf];
X0_bar = [0; 0; 0; 0; 0; 0];

% Initial measurements
data = importdata('homerun_data.txt');
tspan = data(:, 1);
rho = data(:, 2);
alpha = data(:, 3)*pi/180;
beta = data(:, 4)*pi/180;

Y_obs = [];
for i = 1:length(rho)
    Y_obs = vertcat(Y_obs, rho(i));
    Y_obs = vertcat(Y_obs, alpha(i));
    Y_obs = vertcat(Y_obs, beta(i));
end

x_hist = [X0(1)];
y_hist = [X0(2)];
z_hist = [X0(3)];
all_X = [];
all_Y = [];
all_Z = [];
for iter = 1:10
    [~, X_computed] = ode45(@baseball_dynamics, tspan, X0);

    for i = 1:length(X_computed)
        x_i = X_computed(i, 1);
        y_i = X_computed(i, 2);
        z_i = X_computed(i, 3);
    end

    Y_computed = [];
    X_computed = [];

    H = [];
    R = [];
    for t = 0:0.1:1.1
        Phi = [1 0 0 t 0 0;
               0 1 0 0 t 0;
               0 0 1 0 0 t;
               0 0 0 1 0 0;
               0 0 0 0 1 0;
               0 0 0 0 0 1];
        B = [0; 0; -0.5*t^2; 0; 0; -t];
        X_i = Phi * X0 + B * 32.185;
        x = X_i(1);
        y = X_i(2);
        z = X_i(3);
        Htilda = [x/sqrt(x^2+y^2+z^2), y/sqrt(x^2+y^2+z^2), z/sqrt(x^2+y^2+z^2), 0, 0, 0;
                  -y/(x^2*(y^2/x^2 + 1)), 1/(x*(y^2/x^2 + 1)), 0, 0, 0, 0;
                  -(x*z)/((z^2/(x^2 + y^2) + 1)*(x^2 + y^2)^(3/2)), -(y*z)/((z^2/(x^2 + y^2) + 1)*(x^2 + y^2)^(3/2)), 1/((z^2/(x^2 + y^2) + 1)*(x^2 + y^2)^(1/2)), 0, 0, 0];
        H = vertcat(H, Htilda * Phi);
        R = blkdiag(R, R0);
        Y_computed = vertcat(Y_computed, [sqrt(x^2 + y^2 + z^2); atan2(y,x); atan2(z,sqrt(x^2 + y^2))]);
%         X_computed = vertcat(X_computed, 
        [temp_Y temp_X] = ode45(@baseball_dynamics, linspace(0,4), X0);
        all_X = horzcat(all_X, temp_X(:, 1));
        all_Y = horzcat(all_Y, temp_X(:, 2));
        all_Z = horzcat(all_Z, temp_X(:, 3));
    end
    y_dev = Y_obs - Y_computed;
    X0_hat = inv(H'*inv(R)*H + inv(P0_bar)) * (H'*inv(R)*y_dev + inv(P0_bar)*X0_bar);
%     X0_hat = inv(H'*inv(R)*H) * (H'*inv(R)*y_dev);
%     P0_hat = inv(H'*inv(R)*H);
    P0_hat = inv(H'*inv(R)*H + inv(P0_bar));
    
    X0 = X0 + X0_hat;
%     Pk = (eye(6) - K*Htilda * Phi
    P0_bar = P0_hat;
    x_hist = vertcat(x_hist, X0(1));
    y_hist = vertcat(y_hist, X0(2));
    z_hist = vertcat(z_hist, X0(3));
end

figure()
plot(x_hist)
hold on
plot(y_hist)
plot(z_hist)
hold off

figure()
[x y z] = sph2cart(alpha, beta, rho)
scatter3(x, y, z)
hold on
% for i=1:size(all_X, 2)
%     plot3(all_X(:, i), all_Y(:, i), all_Z(:, i))
% end
plot3(all_X(:, end), all_Y(:, end), all_Z(:, end))
xlabel('x (ft)')
ylabel('y (ft)')
zlabel('z (ft)')
title('Trajectory of Baseball')
% plot3(all_X(:, 2), all_Y(:, 2), all_Z(:, 2))
% plot3(all_X(:, 3), all_Y(:, 3), all_Z(:, 3))
% plot3(all_X(:, 6), all_Y(:, 4), all_Z(:, 4))
% legend('s', '1', '2', '3', '4', '5')
hold off

[~, idx] = min(all_Z(:, end).^2);
sqrt(all_X(idx, end)^2 + all_Y(idx, end)^2 + all_Z(idx, end)^2)
norm(P0_bar, 'fro')
%% Functions
function Xdot = baseball_dynamics(t, state)
    g = 32.185;
%     Xdot = [state(4);
%             state(5);
%             state(6) - g*t;
%             0;
%             0;
%             -g];
    Xdot = [state(4);
            state(5);
            state(6);
            0;
            0;
            -g];
end