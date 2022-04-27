clear all; close all;

data = importdata('homerun_data_HW4.txt');
% data(1,:) = [];
tspan = data(:, 1);
rho = data(:, 2)/3.281; % in m
alpha = data(:, 3)*pi/180; % in rad
beta = data(:, 4)*pi/180; % in rad
g = 9.81;
dt = 0.1;
store_x = [];

% Visualizing data
figure(1)
[plotx ploty plotz] = sph2cart(alpha, beta, rho);
scatter3(plotx, ploty, plotz, 10, 'm', 'filled')

% Initial Conditions
x_(:,1) = [0.4921 0.4921 2.0013 -26.2467 114.3051 65.9941]'/3.281;  
% x_(:,1) = [-0.6562; 0.9843; 1.4764;  -32.8084; 119.6219; 55.7806]/3.281;
n = size(x_, 1);
P = diag([4 4 4 0.1 0.1 0.1])/3.281;
R = [(1.524)^2 0 0;
     0 (0.1*pi/180)^2 0;
     0 0 (0.1*pi/180)^2];      %the error covariance constant to be used
Q = diag([0.5^2,0.5^2,0.5^2,0.01^2,0.01^2,0.01^2]);
% Q = zeros(6);
M = eye(3); % COMBAK: is M eye(3) correct?
W = ones(2*n,1) / (2*n); % UKF weights
 
for i =2:length(rho)
    % Observed
    y_obs(:,i) = [rho(i); alpha(i); beta(i)] + sqrt(diag(R)).*randn(size(diag(R)));% + randn*sigmaw;
    
    xhat_k_1 = [];
    % Sigma points
    for j = 1:2*n
        if j<=n
            xtilda = chol(n*P);
            xtilda = xtilda(j,:)';
        else
            xtilda = -chol(n*P);
            xtilda = xtilda(j-n,:)';
        end
        xhat_k_1 = [xhat_k_1 x_(:,i-1)+xtilda];
    end
    
    % Propagation of state
    xhat_k = [];
    for j = 1:2*n
        current_x = xhat_k_1(:,j);
        current_xhatdot = [current_x(4) current_x(5) current_x(6)-g*dt 0 0 -g]';
        xhat_k = [xhat_k current_x + current_xhatdot*dt]; % COMBAK: is rectangular integration ok?
    end
    xhatk_ = sum(xhat_k, 2)/(2*n);
    temp_Pk_ = 0;
    for j = 1:2*n
        temp_Pk_ = temp_Pk_ + (xhat_k(:,j) - xhatk_)*(xhat_k(:,j) - xhatk_)';
    end
    Pk_ = temp_Pk_/(2*n) + Q;
    
    xhat_k = [];
    % Sigma points
    for j = 1:2*n
        if j<=n
            xtilda = chol(n*P);
            xtilda = xtilda(j,:)';
        else
            xtilda = -chol(n*P);
            xtilda = xtilda(j-n,:)';
        end
        xhat_k = [xhat_k xhatk_+xtilda];
    end
    
    % Assembling y_computed
    y_comp = [];
    for j = 1:2*n
        current_x = xhat_k(:,j);
        X_ = current_x(1,end);
        Y_ = current_x(2,end);
        Z_ = current_x(3,end);
        rho_ = sqrt(X_^2 + Y_^2 + Z_^2);
        alpha_ = atan2(Y_, X_);
        beta_ = atan2(Z_, sqrt(X_^2 + Y_^2));

        y_comp = [y_comp [rho_; alpha_; beta_]];
    end
    y_comp_hat = sum(y_comp, 2)/(2*n);

    temp_Py = 0;
    temp_Pxy = 0;
    for j = 1:2*n
        temp_Py = temp_Py + (y_comp(:,j) - y_comp_hat)*(y_comp(:,j) - y_comp_hat)';
        temp_Pxy = temp_Pxy + (xhat_k(:,j)-xhatk_)*(y_comp(:,j)-y_comp_hat)';
    end
    Py = temp_Py/(2*n) + R;
    Pxy = temp_Pxy/(2*n);
    
    K = Pxy*inv(Py);
    x_(:,i) = xhatk_ + K*(y_obs(:,i) - y_comp_hat);
    P = Pk_ - K*Py*K';

    store_x = [store_x x_(:, i)];

    figure(1)
    hold on
    scatter3(x_(1,i), x_(2,i), x_(3,i), 20, 'black')
end

xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
title('Trajectory of Baseball')
legend('Observed', 'Predicted', 'Location', 'best')
hold off;

writematrix(x_, 'baseball_ukf.csv')