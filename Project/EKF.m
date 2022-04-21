clc; clear all; close all;

%% Prepping state data
% Baseball code
% data = importdata('homerun_data_HW4.txt');
% % data(1,:) = [];
% tspan = data(:, 1);
% rho = data(:, 2)/3.281; % in m
% alpha = data(:, 3)*pi/180; % in rad
% beta = data(:, 4)*pi/180; % in rad
% g = 9.81;
% dt = 0.1;
% store_x = [];

% My code
m = 1.545;
Ixx = 0.029125;
Iyy = 0.029125;
Izz = 0.055225;
g = 9.81;

bag = rosbag('2022-04-20-15-32-08.bag');

statesMsgs = readMessages(select(bag, 'Topic', '/lqr_controller/states'),'DataFormat','struct');
states = [];
target_states = [];
time = [];
for i = 1:length(statesMsgs)
    if mod(i,1)==0
        current_state = statesMsgs(i);
        current_state = current_state{1}.Data;
        states = [states current_state(1:12)];
        target_states = [target_states current_state(13:24)];
        time = [time current_state(end)];
    end
end
time = (time)./1e9;
disp("Everything in SI, angles in radians")
% plot(time, states(1,:))
% hold on
% plot(time, target_states(1,:))

% target_states = [];
% time = [];
% statesMsgs = readMessages(select(bag, 'Topic', '/lqr_controller/target_states'),'DataFormat','struct');
% for i = 1:length(statesMsgs)
%     current_state = statesMsgs(i);
%     current_state = current_state{1};
% %     states = [states current_state(1:12)];
%     target_states = [target_states [current_state.Pose.Pose.Position.X; current_state.Pose.Pose.Position.Y; current_state.Pose.Pose.Position.Z]];
%     time = [time current_state.Header.Stamp.Sec + current_state.Header.Stamp.Nsec/1e9 ];
% end
% time = (time);
% plot(time, target_states(1,:))

%% Prepping measurement data
cam1loc = [10, 10, 10]';
cam2loc = [-10, 10, 10]';
cam3loc = [-10, -10, 10]';
cam4loc = [10, -10, 10]';
y1 = vecnorm(states(1:3,:)-cam1loc);
y2 = vecnorm(states(1:3,:)-cam2loc);
y3 = vecnorm(states(1:3,:)-cam3loc);
y4 = vecnorm(states(1:3,:)-cam4loc);

%%
% Baseball code
% % Visualizing data
% figure(1)
% [plotx ploty plotz] = sph2cart(alpha, beta, rho);
% scatter3(plotx, ploty, plotz, 10, 'm', 'filled')

% My code
% Visualizing data
figure(1)
scatter3(states(1,:), states(2,:), states(3,:), 10, 'm', 'filled')

% Baseball code
% % Initial Conditions
% x_(:,1) = [0.4921 0.4921 2.0013 -26.2467 114.3051 65.9941]'/3.281;  
% % x_(:,1) = [-0.6562; 0.9843; 1.4764;  -32.8084; 119.6219; 55.7806]/3.281;
% P = diag([4 4 4 0.1 0.1 0.1])/3.281^2;
% R = [(1.524)^2 0 0;
%      0 (0.1*pi/180)^2 0;
%      0 0 (0.1*pi/180)^2];      %the error covariance constant to be used
% % Q = diag([0.5^2,0.5^2,0.5^2,0.01^2,0.01^2,0.01^2]);
% Q = zeros(6);
% M = eye(3); % COMBAK: is M eye(3) correct?

% My code
% Initial conditions
x_(:,1) = states(:,1); % Using the first column from the data
P = zeros(12); % COMBAK: perfect knowledge of initial state so zero
Q = eye(12); % COMBAK: need to change this later to fit the function
R = eye(4);

for i =2:length(y1)
    dt = time(i) - time(i-1);
%     Baseball code
%   % Observed
%   y_obs(:,i) = [rho(i); alpha(i); beta(i)];% + randn*sigmaw;
%   My code
%   Observed
    y_obs(:,i) = [y1(i); y2(i); y3(i); y4(i)];
  
%     Baseball code
%   % Propagation of state
%   xhatdot = [x_(4,i-1) x_(5,i-1) x_(6,i-1)-g*dt 0 0 -g]';
%   x_(:,i) = x_(:,i-1) + xhatdot*dt; % COMBAK: is rectangular integration ok?
  
%   My code
%   Propagation of state
    x_(:,i) = propagate_state(x_(:,i-1), target_states(:,i-1), dt, m, Ixx, Iyy, Izz, g);
    
%     My code
    % Propagation of state covariance
    A = find_A(x_(:,i-1), m, Ixx, Iyy, Izz, g);
    Pdot = A*P + P*A' + Q; % + LQL' COMBAK: I didn't use Q here bc we weren't given one. OK?
    P = P + Pdot*dt;
  
%     Baseball code
%     % Assembling y_computed
%     X_ = x_(1,end);
%     Y_ = x_(2,end);
%     Z_ = x_(3,end);
%     rho_ = sqrt(X_^2 + Y_^2 + Z_^2);
%     alpha_ = atan2(Y_, X_);
%     beta_ = atan2(Z_, sqrt(X_^2 + Y_^2));
%   
%   y_comp(:,i) = [rho_; alpha_; beta_];

%     My code
%   Assembling y_comp
    y_comp(:,i) = [vecnorm(x_(1:3,end)-cam1loc);
                   vecnorm(x_(1:3,end)-cam2loc);
                   vecnorm(x_(1:3,end)-cam3loc);
                   vecnorm(x_(1:3,end)-cam4loc)];
%     Baseball code
%   % Computing H using y_comp
%   H = [X_/rho_ Y_/rho_ Z_/rho_ 0 0 0;
%        (-Y_/X_^2)/(1 + (Y_/X_)^2) (1/X_)/(1+(Y_/X_)^2) 0 0 0 0;
%        ((-X_*Z_)/(X_^2 + Y_^2)^(3/2))/(1+(Z_^2)/(X_^2 + Y_^2)) ((-Y_*Z_)/(X_^2 + Y_^2)^(3/2))/(1 + (Z_^2)/(X_^2 + Y_^2)) ((1)/(X_^2 + Y_^2)^(1/2))/(1 + (Z_^2)/(X_^2 + Y_^2)) 0 0 0];

%   My code
%   Computing H using y_comp
    X1 = cam1loc(1);
    Y1 = cam1loc(2);
    Z1 = cam1loc(3);
    X2 = cam2loc(1);
    Y2 = cam2loc(2);
    Z2 = cam2loc(3);
    X3 = cam3loc(1);
    Y3 = cam3loc(2);
    Z3 = cam3loc(3);
    X4 = cam4loc(1);
    Y4 = cam4loc(2);
    Z4 = cam4loc(3);
    drone_x = x_(1,end);
    drone_y = x_(2,end);
    drone_z = x_(3,end);
    H = [-(X1 - drone_x)/((X1 - drone_x)^2 + (Y1 - drone_y)^2 + (Z1 - drone_z)^2)^(1/2), -(Y1 - drone_y)/((X1 - drone_x)^2 + (Y1 - drone_y)^2 + (Z1 - drone_z)^2)^(1/2), -(Z1 - drone_z)/((X1 - drone_x)^2 + (Y1 - drone_y)^2 + (Z1 - drone_z)^2)^(1/2), 0, 0, 0, 0, 0, 0, 0, 0, 0;
         -(X2 - drone_x)/((X2 - drone_x)^2 + (Y2 - drone_y)^2 + (Z2 - drone_z)^2)^(1/2), -(Y2 - drone_y)/((X2 - drone_x)^2 + (Y2 - drone_y)^2 + (Z2 - drone_z)^2)^(1/2), -(Z2 - drone_z)/((X2 - drone_x)^2 + (Y2 - drone_y)^2 + (Z2 - drone_z)^2)^(1/2), 0, 0, 0, 0, 0, 0, 0, 0, 0;
         -(X3 - drone_x)/((X3 - drone_x)^2 + (Y3 - drone_y)^2 + (Z3 - drone_z)^2)^(1/2), -(Y3 - drone_y)/((X3 - drone_x)^2 + (Y3 - drone_y)^2 + (Z3 - drone_z)^2)^(1/2), -(Z3 - drone_z)/((X3 - drone_x)^2 + (Y3 - drone_y)^2 + (Z3 - drone_z)^2)^(1/2), 0, 0, 0, 0, 0, 0, 0, 0, 0;
         -(X4 - drone_x)/((X4 - drone_x)^2 + (Y4 - drone_y)^2 + (Z4 - drone_z)^2)^(1/2), -(Y4 - drone_y)/((X4 - drone_x)^2 + (Y4 - drone_y)^2 + (Z4 - drone_z)^2)^(1/2), -(Z4 - drone_z)/((X4 - drone_x)^2 + (Y4 - drone_y)^2 + (Z4 - drone_z)^2)^(1/2), 0, 0, 0, 0, 0, 0, 0, 0, 0];
 
% Baseball code
    % Kalman gain
%     K = P*H'*inv(H*P*H'+ R);
%     % Measurement update
%     x_(:,i) = x_(:,i) + K * (y_obs(:,i) - y_comp(:,i));
%     P = (eye(6)-K*H)*P*(eye(6)-K*H)' + K*R*K';

%   My code
%   Kalman gain
    K = P*H'*inv(H*P*H'+R);
    
    % Measurement update
    x_(:,i) = x_(:,i) + K * (y_obs(:,i) - y_comp(:,i));
    P = (eye(12)-K*H)*P*(eye(12)-K*H)' + K*R*K';
  
%     store_x = [store_x x_(:, i)];

    figure(1)
    hold on
    scatter3(x_(1,end), x_(2,end), x_(3,end), 20, 'black')
end

xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
title('Trajectory of Baseball')
legend('Observed', 'Predicted', 'Location', 'best')
hold off;
writematrix(x_, 'baseball_ekf.csv')

%% Functions
function [u3 u4 u5 u6] = get_u(states, target_states, m, Ixx, Iyy, Izz, g)
    q7 = states(7);
    q8 = states(8);
    q9 = states(9);
    q10 = states(10);
    q11 = states(11);
    q12 = states(12);
    u3 = m*g;
    A = find_A(states, m, Ixx, Iyy, Izz, g);
    B = [0, 0, 0, 0;
         0, 0, 0, 0;
         0, 0, 0, 0;
         (sin(q7) * sin(q9) + cos(q7) * cos(q9) * sin(q8)) / m, 0, 0, 0;
         -(cos(q9) * sin(q7) - cos(q7) * sin(q8) * sin(q9)) / m, 0, 0, 0;
         (cos(q7) * cos(q8)) / m, 0, 0, 0;
         0, 0, 0, 0;
         0, 0, 0, 0;
         0, 0, 0, 0;
         0, 1 / Ixx, 0, 0;
         0, 0, 1 / Iyy, 0;
         0, 0, 0, 1 / Izz];
    Q = diag([100 100 100 100 100 100 0.1 0.1 0.1 10 10 10]);
    R = diag([0.1 1000 1000 1000]);
    K = lqr(A, B, Q, R);
    u = -K*(states - target_states);
    u(1) = u(1) + m*g;
    u3 = u(1);
    u4 = u(2);
    u5 = u(3);
    u6 = u(4);
end

function A = find_A(states, m, Ixx, Iyy, Izz, g)
    q7 = states(7);
    q8 = states(8);
    q9 = states(9);
    q10 = states(10);
    q11 = states(11);
    q12 = states(12);
    u3 = m*g;
    A = [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0;
         0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0;
         0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0;
         0, 0, 0, 0, 0, 0, (u3 * (cos(q7) * sin(q9) - cos(q9) * sin(q7) * sin(q8))) / m, ...
                           (u3 * cos(q7) * cos(q8) * cos(q9)) / m, ...
                           (u3 * (cos(q9) * sin(q7) - cos(q7) * sin(q8) * sin(q9))) / m, ...
                           0, 0, 0;
         0, 0, 0, 0, 0, 0, -(u3 * (cos(q7) * cos(q9) + sin(q7) * sin(q8) * sin(q9))) / m, ...
                           (u3 * cos(q7) * cos(q8) * sin(q9)) / m, ...
                           (u3 * (sin(q7) * sin(q9) + cos(q7) * cos(q9) * sin(q8))) / m, ...
                           0, 0, 0;
         0, 0, 0, 0, 0, 0, -(u3 * cos(q8) * sin(q7)) / m, ...
                           -(u3 * cos(q7) * sin(q8)) / m, ...
                           0, 0, 0, 0;
         0, 0, 0, 0, 0, 0, q11 * cos(q7) * tan(q8) - q12 * sin(q7) * tan(q8), ...
                           q12 * cos(q7) * (tan(q8)^2 + 1) + q11 * sin(q7) * (tan(q8)^2 + 1), ...
                           0, 1, sin(q7) * tan(q8), cos(q7) * tan(q8);
         0, 0, 0, 0, 0, 0, -q12 * cos(q7) - q11 * sin(q7), 0, 0, 0, cos(q7), -sin(q7);
         0, 0, 0, 0, 0, 0, (q11 * cos(q7)) / cos(q8) - (q12 * sin(q7)) / cos(q8),  ...
                           (q12 * cos(q7) * sin(q8)) / cos(q8)^2 + (q11 * sin(q7) * sin(q8)) / cos(q8)^2, ...
                           0, 0, sin(q7) / cos(q8), cos(q7) / cos(q8);
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (Iyy * q12 - Izz * q12) / Ixx, ...
                                       (Iyy * q11 - Izz * q11) / Ixx;
         0, 0, 0, 0, 0, 0, 0, 0, 0, -(Ixx * q12 - Izz * q12) / Iyy, ...
                                    0, -(Ixx * q10 - Izz * q10) / Iyy;
         0, 0, 0, 0, 0, 0, 0, 0, 0, (Ixx * q11 - Iyy * q11) / Izz, ...
                                    (Ixx * q10 - Iyy * q10) / Izz, 0];
end

function x = propagate_state(last_state, target_state, dt, m, Ixx, Iyy, Izz, g)
    x = last_state;
    last_t = 0;
    for t = linspace(1e-99,dt,20)
        x_pos = x(1);
        y_pos = x(2);
        z_pos = x(3);
        x_vel = x(4);
        y_vel = x(5);
        z_vel = x(6);
        roll = x(7);
        pitch = x(8);
        yaw = x(9);
        roll_rate = x(10);
        pitch_rate = x(11);
        yaw_rate = x(12);
        [u3, u4, u5, u6] = get_u(x, target_state, m, Ixx, Iyy, Izz, g);

        xhatdot =[x_vel;
                  y_vel;
                  z_vel;
                  (u3*(sin(roll)*sin(yaw) + cos(roll)*cos(yaw)*sin(pitch)))/m;
                  -(u3*(cos(yaw)*sin(roll) - cos(roll)*sin(pitch)*sin(yaw)))/m;
                  -(g*m - u3*cos(pitch)*cos(roll))/m;
                  roll_rate + yaw_rate*cos(roll)*tan(pitch) + pitch_rate*tan(pitch)*sin(roll);
                  pitch_rate*cos(roll) - yaw_rate*sin(roll);
                  (yaw_rate*cos(roll))/cos(pitch) + (pitch_rate*sin(roll))/cos(pitch);
                  (u4 + Iyy*pitch_rate*yaw_rate - Izz*pitch_rate*yaw_rate)/Ixx;
                  (u5 - Ixx*roll_rate*yaw_rate + Izz*roll_rate*yaw_rate)/Iyy;
                  (u6 + Ixx*pitch_rate*roll_rate - Iyy*pitch_rate*roll_rate)/Izz];

        x = x + xhatdot*(t-last_t);
        last_t = t;
    end
end