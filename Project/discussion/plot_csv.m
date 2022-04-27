%%
clc; clear all; close all;

%%
data = readmatrix('C:\Users\kumasan\Downloads\2021-06-07-10-26-42.csv');

%%
t = data(1,:);
x_pos = data(2,:);
y_pos = data(3,:);
z_pos = data(4,:);
x_vel = data(5,:);
y_vel = data(6,:);
z_vel = data(7,:);
roll = data(8,:);
pitch = data(9,:);
yaw = data(10,:);
p = data(11,:);
q = data(12,:);
r = data(13,:);
thrust = data(14,:);
roll_torque = data(15,:);
pitch_torque = data(16,:);
yaw_torque = data(17,:);
w_setpoint_t = data(18, :);
w1_setpoint = data(19, :);
w2_setpoint = data(20, :);
w3_setpoint = data(21, :);
w4_setpoint = data(22, :);
% thrust,roll_torque,pitch_torque,yaw_torque,w_setpoint_t,w1_setpoint,w2_setpoint,w3_setpoint,w4_setpoint,w_t,w1,w2,w3,w4

hold on
% plot(t, x_pos, 'DisplayName', 'x pos')
% plot(t, y_pos, 'DisplayName', 'y pos')
% plot(t, z_pos, 'DisplayName', 'z pos')
% plot(t, x_vel, 'DisplayName', 'x vel')
% plot(t, y_vel, 'DisplayName', 'y vel')
% plot(t, z_vel, 'DisplayName', 'z vel')
% plot(t, roll, 'DisplayName', 'roll')
% plot(t, pitch, 'DisplayName', 'pitch')
% plot(t, yaw, 'DisplayName', 'yaw')
% plot(t, p, 'DisplayName', 'p')
% plot(t, q, 'DisplayName', 'q')
% plot(t, r, 'DisplayName', 'r')
plot(t, thrust, 'DisplayName', 'thrust')
% plot(t, roll_torque, 'DisplayName', 'roll torque')
% plot(t, pitch_torque, 'DisplayName', 'pitch torque')
% plot(t, yaw_torque, 'DisplayName', 'yaw torque')
legend
grid on
hold off

%% Test: finding prop speeds
u = [thrust; roll_torque; pitch_torque; yaw_torque];
zero_vec = zeros(1, size(u, 2));
% u(1, :) = zero_vec;
% u(2, :) = zero_vec;
% u(3, :) = zero_vec;
% u(4, :) = zero_vec;

Ct = 5.84e-6;
Cq = 1.75e-4;
Lrf = 0.22;
Lrb = 0.2;
Lp = 0.13;

% M from docs
% M = [Ct, Ct, Ct, Ct;
%      -Ct*Lrf, Ct*Lrb, Ct*Lrf, -Ct*Lrb;
%      Ct*Lp, -Ct*Lp, Ct*Lp, -Ct*Lp;
%      Cq, Cq, -Cq, -Cq];
 
% M after flipping pitch signs
% M = [Ct, Ct, Ct, Ct;
%      -Ct*Lrf, Ct*Lrb, Ct*Lrf, -Ct*Lrb;
%      -Ct*Lp, Ct*Lp, -Ct*Lp, Ct*Lp;
%      Cq, Cq, -Cq, -Cq];
 
% M after flipping yaw signs
M = [Ct, Ct, Ct, Ct;
     -Ct*Lrf, Ct*Lrb, Ct*Lrf, -Ct*Lrb;
     -Ct*Lp, Ct*Lp, -Ct*Lp, Ct*Lp;
     -Cq, -Cq, Cq, Cq];

inv_M = inv(M)
 
for i = 1:size(u, 2)
    prop_speeds(:, i) = (inv_M * u(:, i));
end

prop_speeds = max(prop_speeds, 100^2);
prop_speeds = min(prop_speeds, 1100^2);
prop_speeds = sqrt(prop_speeds);

figure(3)
sgtitle('Blue=MATLAB, Orange=ROS')
subplot(2,2,1)
plot(t, prop_speeds(3, :), 'LineWidth', 1)
hold on
plot(w_setpoint_t, w3_setpoint, 'LineWidth', 1)
yline(1100, '--g')
yline(100, '--r')
hold off
subplot(2,2,2)
plot(t, prop_speeds(1, :), 'LineWidth', 1)
hold on
plot(w_setpoint_t, w1_setpoint, 'LineWidth', 1)
yline(1100, '--g')
yline(100, '--r')
hold off
subplot(2,2,3)
plot(t, prop_speeds(2, :), 'LineWidth', 1)
hold on
plot(w_setpoint_t, w2_setpoint, 'LineWidth', 1)
yline(1100, '--g')
yline(100, '--r')
hold off
subplot(2,2,4)
plot(t, prop_speeds(4, :), 'LineWidth', 1)
hold on
plot(w_setpoint_t, w4_setpoint, 'LineWidth', 1)
yline(1100, '--g')
yline(100, '--r')
hold off