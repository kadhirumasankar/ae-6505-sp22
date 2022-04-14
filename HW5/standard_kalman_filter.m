clear all; close all;

data = importdata('homerun_data_HW4.txt');
% data(1,:) = [];
tspan = data(:, 1);
rho = data(:, 2)/3.281; % in m
alpha = data(:, 3)*pi/180; % in rad
beta = data(:, 4)*pi/180; % in rad
g = 9.81;
store_x = [];

% Visualizing data
figure(1)
[plotx ploty plotz] = sph2cart(alpha, beta, rho);
scatter3(plotx, ploty, plotz, 10, 'm', 'filled')

% Initial Conditions
x_(:,1) = [0.4921 0.4921 2.0013 -26.2467 114.3051 65.9941]'/3.281;  
% x_(:,1) = zeros(6,1);% Given initial conditions
P = diag([4 4 4 0.1 0.1 0.1])/3.281^2;
R = [(1.524)^2 0 0;
     0 (0.1*pi/180)^2 0;
     0 0 (0.1*pi/180)^2];      %the error covariance constant to be used

for i =2:length(rho)
  % Observed
  y_obs(:,i) = [rho(i); alpha(i); beta(i)];% + randn*sigmaw;

  dt = 0.1;
  Phi = [1 0 0 dt 0 0;
         0 1 0 0 dt 0;
         0 0 1 0 0 dt;
         0 0 0 1 0 0;
         0 0 0 0 1 0;
         0 0 0 0 0 1];
        
  % Compute nominal trajectory
  B = [0; 0; -1/2*dt^2; 0; 0; -dt];
  x_(:,i) = Phi * x_(:, i-1) + B * g;
  
  % Assembling y_computed
  X_ = x_(1,end);
  Y_ = x_(2,end);
  Z_ = x_(3,end);
  rho_ = sqrt(X_^2 + Y_^2 + Z_^2);
  alpha_ = atan2(Y_, X_);
  beta_ = atan2(Z_, sqrt(X_^2 + Y_^2));
  
  y_comp(:,i) = [rho_; alpha_; beta_];

  % Propagation of state covariance
  P = Phi*P*Phi'; % Add Q here to improve tracking

  % Computing H using y_comp
  H = [X_/rho_ Y_/rho_ Z_/rho_ 0 0 0;
       (-Y_/X_^2)/(1 + (Y_/X_)^2) (1/X_)/(1+(Y_/X_)^2) 0 0 0 0;
       ((-X_*Z_)/(X_^2 + Y_^2)^(3/2))/(1+(Z_^2)/(X_^2 + Y_^2)) ((-Y_*Z_)/(X_^2 + Y_^2)^(3/2))/(1 + (Z_^2)/(X_^2 + Y_^2)) ((1)/(X_^2 + Y_^2)^(1/2))/(1 + (Z_^2)/(X_^2 + Y_^2)) 0 0 0];
  % Kalman gain
  K = P*H'*inv(H*P*H'+R);
  % Measurement update
  H*x_(:,i)
  x_(:,i) = x_(:,i) + K * (y_obs(:,i) - y_comp(:,i));
  P = (eye(6)-K*H)*P;
  
  store_x = [store_x x_(:, i)];

  figure(1)
  hold on
  scatter3(X_, Y_, Z_, 20, 'black')
end

xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
title('Trajectory of Baseball')
legend('Observed', 'Predicted', 'Location', 'best')
hold off;
