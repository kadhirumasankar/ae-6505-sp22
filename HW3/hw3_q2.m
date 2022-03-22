clc; clear all; close all;

% Given info
h = 6.4;
k1 = 2.6;
k2 = 3.5;
m = 1.8;
w = sqrt((k1+k2)/m);
% Measurement noise covariance
R0 = [0.0625 0;
     0 0.01];
% Initial state
X0 = [3.9;
      0.5];
tspan = 0:10;
P0_bar = [1000 0;
          0 100];
X0_bar = [0; 0];

% Initial measurements
Y_obs = [];
rho = [6.3118 6.0217 5.1295 6.3685 5.5422 5.5095 5.9740 5.4930 6.8653 6.6661 5.0692]';
rhoDot = [0.3035 1.3854 -1.5512 0.6064 0.8642 -1.5736 1.1588 0.4580 -1.2328 1.4348 -0.4229]';
for i = 1:length(rho)
    Y_obs = vertcat(Y_obs, rho(i));
    Y_obs = vertcat(Y_obs, rhoDot(i));
end

x_hist = [X0(1)];
v_hist = [X0(2)];
for iter = 1:50
    [~, X_computed] = ode45(@spring_mass_dynamics, tspan, X0);
    
    Y_computed = [];
    for i = 1:length(X_computed)
        c_i = X_computed(i, 1);
        d_i = X_computed(i, 2);
        rho_i = sqrt(c_i^2 + h^2);
        Y_computed = vertcat(Y_computed, [rho_i; c_i*d_i/rho_i]);
    end

    y_dev = Y_obs - Y_computed;
    
    H = [];
    R = [];
    for t = 0:10
        Phi = [cos(w*t) 1/w*sin(w*t);
               -w*sin(w*t) cos(w*t)];
        X_i = Phi * X0;
        c_i = X_i(1);
        d_i = X_i(2);
        rho_i = sqrt(c_i^2 + h^2);
        Htilda = [c_i/rho_i 0;
                  d_i/rho_i - c_i^2*d_i/rho_i^3 c_i/rho_i];
        H = vertcat(H, Htilda * Phi);
        R = blkdiag(R, R0);
    end
    X0_hat = inv(H'*inv(R)*H + inv(P0_bar)) * (H'*inv(R)*y_dev + inv(P0_bar)*X0_bar);
    P0_hat = inv(H'*inv(R)*H + inv(P0_bar));
    
    X0 = X0 + X0_hat;
    P0_bar = P0_hat
    x_hist = vertcat(x_hist, X0(1));
    v_hist = vertcat(v_hist, X0(2));
end

plot(0:length(x_hist)-1, x_hist, 'LineWidth', 1.5)
hold on
plot(0:length(v_hist)-1, v_hist, 'LineWidth', 1.5)
title('X_0 Estimate over Multiple Iterations')
legend('x_0 (in m)', 'v_0 (in m/s)')
xlabel('Iterations')
grid on
box on
hold off

%% Functions
function Xdot = spring_mass_dynamics(t, state)
    k1 = 2.6;
    k2 = 3.5;
    m = 1.8;
    w = sqrt((k1+k2)/m);
    Xdot = [state(2);
            -w^2*state(1)];
end