clc; clear all; close all;

data = importdata('homerun_data_HW4.txt');
% data(1,:) = [];
tspan = data(:, 1);
rho = data(:, 2)/3.281; % in m
alpha = data(:, 3)*pi/180; % in rad
beta = data(:, 4)*pi/180; % in rad
batch_data = importdata('baseball_batch.csv');
kf_data = importdata('baseball_kf.csv');
ekf_data = importdata('baseball_ekf.csv');
ukf_data = importdata('baseball_ukf.csv');
pf_data = importdata('baseball_pf.csv');

true_data = [];
t_all = tspan;
X_0 = [-0.6562; 0.9843; 1.4764;  -32.8084; 119.6219; 55.7806]/3.281;
g = 9.81;
for k=1:length(t_all)
    X0 = X_0(1);
    Y0 = X_0(2);
    Z0 = X_0(3);
    Vx0 = X_0(4);
    Vy0 = X_0(5);
    Vz0 = X_0(6);
    X = X0 + Vx0 * t_all(k);
    Y = Y0 + Vy0 * t_all(k);
    Z = Z0 + Vz0 * t_all(k) - (1/2)*g*t_all(k)^2;
    Vx = Vx0;
    Vy = Vy0;
    Vz = Vz0 - g*t_all(k);
    true_data = [true_data [X; Y; Z; Vx; Vy; Vz]];
end

% Visualizing data
figure(1)
[plotx ploty plotz] = sph2cart(alpha, beta, rho);
scatter3(plotx, ploty, plotz, 10, 'm', 'filled')
hold on
scatter3(ekf_data(1,:),ekf_data(2,:),ekf_data(3,:), 20, 'black')
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
title('Trajectory of Baseball')
legend('Observed', 'Extended KF', 'Location', 'best')
hold off;

figure(2)
scatter3(plotx, ploty, plotz, 10, 'm', 'filled')
hold on
scatter3(ukf_data(1,:),ukf_data(2,:),ukf_data(3,:), 20, 'black')
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
title('Trajectory of Baseball')
legend('Observed', 'Unscented KF', 'Location', 'best')
hold off;

figure(3)
scatter3(plotx, ploty, plotz, 10, 'm', 'filled')
hold on
scatter3(pf_data(1,:),pf_data(2,:),pf_data(3,:), 20, 'black')
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
title('Trajectory of Baseball')
legend('Observed', 'Particle Filter', 'Location', 'best')
hold off;

figure(4)
scatter3(true_data(1,:),true_data(2,:),true_data(3,:), 10, 'm', 'filled')
hold on
scatter3(batch_data(1,:),batch_data(2,:),batch_data(3,:), 10, 's', 'filled', 'red')
scatter3(kf_data(1,:),kf_data(2,:),kf_data(3,:), 10, 'd', 'filled', 'green')
scatter3(ekf_data(1,:),ekf_data(2,:),ekf_data(3,:), 20, 'p', 'filled', 'cyan')
scatter3(ukf_data(1,:),ukf_data(2,:),ukf_data(3,:), 10, '^', 'filled', 'black')
scatter3(pf_data(1,:),pf_data(2,:),pf_data(3,:), 10, 'h', 'filled', 'blue')
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
title('Trajectory of Baseball')
legend('True', 'Batch Estimation', 'Standard Kalman Filter', 'Extended Kalman Filter', 'Unscented Kalman Filter', 'Particle Filter', 'Location', 'best')
hold off;

figure(5)
% subplot(3,1,1)
sgtitle('Error from Truth')
plot(t_all, vecnorm(batch_data(1:3,:)-true_data(1:3,:)))
hold on
plot(t_all, vecnorm(kf_data(1:3,:)-true_data(1:3,:)))
plot(t_all, vecnorm(ekf_data(1:3,:)-true_data(1:3,:)))
plot(t_all, vecnorm(ukf_data(1:3,:)-true_data(1:3,:)))
plot(t_all, vecnorm(pf_data(1:3,:)-true_data(1:3,:)))
ylabel('Error magnitude (ft)')
xlabel('Time (s)')
legend('Batch Estimation', 'Standard Kalman Filter', 'Extended Kalman Filter', 'Unscented Kalman Filter', 'Particle Filter', 'Location', 'best')
hold off

fprintf("Std Dev of Estimation Error of position with batch estimation \t%g\n", std(vecnorm(batch_data(1:3,:)-true_data(1:3,:))))
fprintf("Std Dev of Estimation Error of position with standard KF \t\t%g\n", std(vecnorm(kf_data(1:3,:)-true_data(1:3,:))))
fprintf("Std Dev of Estimation Error of position with EKF \t\t\t\t%g\n", std(vecnorm(ekf_data(1:3,:)-true_data(1:3,:))))
fprintf("Std Dev of Estimation Error of position with UKF \t\t\t\t%g\n", std(vecnorm(ukf_data(1:3,:)-true_data(1:3,:))))
fprintf("Std Dev of Estimation Error of position with particle filter \t%g\n", std(vecnorm(pf_data(1:3,:)-true_data(1:3,:))))