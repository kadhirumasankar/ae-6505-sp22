% Hybrid extended Kalman filter example.
% Track a body falling through the atmosphere.
% This example is taken from [Jul00], which was based on [Ath68].

clc; clear all; close all;

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

x = [0.4921 0.4921 2.0013 -26.2467 114.3051 65.9941]'/3.281;
xhat = [0.4921 0.4921 2.0013 -26.2467 114.3051 65.9941]'/3.281;
xhatukf = xhat;

P = diag([4 4 4 0.1 0.1 0.1])/3.281^2;
Pukf = P;
R = [(1.524)^2 0 0;
     0 (0.1*pi/180)^2 0;
     0 0 (0.1*pi/180)^2];      %the error covariance constant to be used
Q = eye(6) * 0; % process noise covariance
Q = diag([0.5^2,0.5^2,0.5^2,0.01^2,0.01^2,0.01^2]);

T = tspan(2)-tspan(1); % measurement time step
randn('state',sum(100*clock)); % random number generator seed

tf = tspan(end); % simulation length (seconds)
xArray = x;
xhatArray = xhat;
xhatukfArray = xhatukf;
Parray = diag(P);
Pukfarray = diag(Pukf);

W = ones(12,1) / 12; % UKF weights

counter = 2;
for t = 2:length(tspan)
   % Simulate the noisy measurement.
   z = [rho(t); alpha(t); beta(t)];
   
   [root,p] = chol(3*Pukf);
   
   for i = 1 : 6
       sigma(:,i) = xhatukf + root(i,:)';
       sigma(:,i+6) = xhatukf - root(i,:)';
   end
   for i = 1 : 12
       xbreve(:,i) = sigma(:,i);
   end
   % UKF time update
   for i = 1 : 12
%        for tau = dt : dt : T
%           xbrevedot = [xbreve(4, i) xbreve(5, i) xbreve(6, i)-g*dt 0 0 -g]';
%           xbreve(:,i) = xbreve(:,i) + xbrevedot * dt;
%       end
        Phi = [1 0 0 dt 0 0;
         0 1 0 0 dt 0;
         0 0 1 0 0 dt;
         0 0 0 1 0 0;
         0 0 0 0 1 0;
         0 0 0 0 0 1];
        B = [0; 0; -1/2*dt^2; 0; 0; -dt];
        xbreve(:,i) = Phi * xbreve(:, i) + B * g;
  end
  xhatukf = zeros(6,1);
  for i = 1 : 12
      xhatukf = xhatukf + W(i) * xbreve(:,i);
  end
  Pukf = zeros(6,6);
  for i = 1 : 12
      Pukf = Pukf + W(i) * (xbreve(:,i) - xhatukf) * (xbreve(:,i) - xhatukf)';
  end
  Pukf = Pukf + Q;
  % UKF measurement update
  for i = 1 : 12
      current_x = xbreve(:,i);
      X_ = current_x(1);
      Y_ = current_x(2);
      Z_ = current_x(3);
      rho_ = sqrt(X_^2 + Y_^2 + Z_^2);
      alpha_ = atan2(Y_, X_);
      beta_ = atan2(Z_, sqrt(X_^2 + Y_^2));
      zukf(:,i) = [rho_; alpha_; beta_];
  end
  zhat = 0;
  for i = 1 : 12
      zhat = zhat + W(i) * zukf(:,i);
  end
  zhat
  Py = 0;
  Pxy = zeros(6,1);
  for i = 1 : 12
      Py = Py + W(i) * (zukf(:,i) - zhat) * (zukf(:,i) - zhat)';
      Pxy = Pxy + W(i) * (xbreve(:,i) - xhat) * (zukf(:,i) - zhat)';
  end
  Py = Py + R;
  Kukf = Pxy * inv(Py);
  xhatukf = xhatukf + Kukf * (z - zhat);
  Pukf = Pukf - Kukf * Py * Kukf';      
   
   % Save data for plotting.
   xArray = [xArray x];
   xhatArray = [xhatArray xhat];
   xhatukfArray = [xhatukfArray xhatukf];
   Parray = [Parray diag(P)];
   Pukfarray = [Pukfarray diag(Pukf)];
   counter = counter + 1;
   
   figure(1)
   hold on
   scatter3(xhatukf(1), xhatukf(2), xhatukf(3), 20, 'black')
end

%{
close all;
t = 0 : T : tf;
figure; 
semilogy(t, abs(xArray(1,:) - xhatArray(1,:)), 'b'); hold;
%plot(t, sqrt(Parray(1,:)), 'b--');
semilogy(t, abs(xArray(1,:) - xhatukfArray(1,:)), 'r:');
%plot(t, sqrt(Pukfarray(1,:)), 'r--');
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Seconds');
ylabel('Position Estimation Error');
legend('Kalman filter', 'Unscented filter');

figure; 
semilogy(t, abs(xArray(2,:) - xhatArray(2,:)), 'b'); hold;
%plot(t, sqrt(Parray(2,:)), 'b--');
semilogy(t, abs(xArray(2,:) - xhatukfArray(2,:)), 'r:');
%plot(t, sqrt(Pukfarray(2,:)), 'r--');
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Seconds');
ylabel('Velocity Estimation Error');
legend('Kalman filter', 'Unscented filter');

figure; 
semilogy(t, abs(xArray(3,:) - xhatArray(3,:)), 'b'); hold;
%plot(t, sqrt(Parray(3,:)), 'b--');
semilogy(t, abs(xArray(3,:) - xhatukfArray(3,:)), 'r:');
%plot(t, sqrt(Pukfarray(3,:)), 'r--');
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Seconds');
ylabel('Ballistic Coefficient Estimation Error');
legend('Kalman filter', 'Unscented filter');

figure;
plot(t, xArray(1,:));
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Seconds');
ylabel('True Position');

figure;
plot(t, xArray(2,:));
title('Falling Body Simulation', 'FontSize', 12);
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Seconds');
ylabel('True Velocity');
%}