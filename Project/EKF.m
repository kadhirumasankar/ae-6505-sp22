clc; clear all; close all;

%% Prepping state data

bag_list = ["circle", "hover", "sine", "square"];
% bag_list = ["hover"];
freq_list = [1 2 5 10 20 60 120];
%freq_list = [20 60 120];

for bag_idx = 1:length(bag_list)
for freq_idx = 1:length(freq_list)
clc; close all; clearvars -except bag_idx freq_idx bag_list freq_list
m = 1.545;
Ixx = 0.029125;
Iyy = 0.029125;
Izz = 0.055225;
g = 9.81;
shape = bag_list(bag_idx);
sample_freq = freq_list(freq_idx);
bag = rosbag(strcat(shape, '.bag'));
sample_freq
shape

statesMsgs = readMessages(select(bag, 'Topic', '/lqr_controller/states'),'DataFormat','struct');
states = [];
target_states = [];
marker_locs = [];
time = [];
for i = 1:length(statesMsgs)
    current_state = statesMsgs(i);
    current_state = current_state{1}.Data;
    if current_state(3) < 0 && isempty(states)
        
    else
        if mod(i,sample_freq)==0
            states = [states current_state(1:12)];
            target_states = [target_states current_state(13:24)];
            marker_locs = [marker_locs get_marker_locs(current_state(1:12))];
            time = [time current_state(end)];
        end
    end
    if i > length(statesMsgs)/1
        break
    end
end
time = (time)./1e9;
disp("Everything in SI, angles in radians")


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
% Visualizing data
figure(1)
scatter3(states(1,:), states(2,:), states(3,:), 10, 'm', 'filled')
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
title('Trajectory of Drone')
legend('Observed', 'Predicted', 'Location', 'best')
if shape == "hover"
    xlim([-1 1])
    ylim([-1 1])
end
if sample_freq == 1
saveas(gcf, strcat(shape, num2str(sample_freq), '_just_traj.png'))
end
% hold on
% scatter3(marker_locs(1,:), marker_locs(2,:), marker_locs(3,:), 5, 'r', 'filled')
% scatter3(marker_locs(4,:), marker_locs(5,:), marker_locs(6,:), 5, 'b', 'filled')
% scatter3(marker_locs(7,:), marker_locs(8,:), marker_locs(9,:), 5, 'g', 'filled')
% scatter3(marker_locs(10,:), marker_locs(11,:), marker_locs(12,:), 5, 'k', 'filled')
% scatter3(marker_locs(13,:), marker_locs(14,:), marker_locs(15,:), 5, 'k', 'filled')

% Initial conditions
x_(:,1) = states(:,1); % Using the first column from the data
x_(4:6, 1) = [0; 0; 0];
x_(10:12, 1) = [0; 0; 0];
P = diag([(0.001)^2 (0.001)^2 (0.001)^2 (0.001)^2 (0.001)^2 (0.001)^2 (3.8785e-5)^2 (3.8785e-5)^2 (3.8785e-5)^2 (3.8785e-5)^2 (3.8785e-5)^2 (3.8785e-5)^2]);
Q = eye(12).*1e-3;
R = diag([0.0015^2 0.015^2 0.002^2 0.1^2 (3.8785e-5)^2 (3.8785e-5)^2 (3.8785e-5)^2]);

store_x = [x_(:,1)];
store_P = [norm(diag(P))];

timeElapsedArray = [];
for i=2:length(y1)
    tic;
    dt = time(i) - time(i-1);
%   Observed
    y_obs(:,i) = [y1(i); y2(i); y3(i); y4(i); states(7,i); states(8,i); states(9,i)] + sqrt(diag(R)).*randn(size(diag(R)));
  
%   Propagation of state
    if dt ~= 0
        [t_out, y_out] = ode45(@(t,y) drone_dynamics(t, y, target_states(:,i-1), m, Ixx, Iyy, Izz, g), [0 dt], x_(:, i-1), odeset('RelTol',1e-2,'AbsTol',1e-4));
    end

    x_(:,i) = y_out(end,:)';

    % Propagation of state covariance
    A = find_A(x_(:,i-1), m, Ixx, Iyy, Izz, g);
    Pdot = A*P + P*A' + Q; % + LQL' COMBAK: I didn't use Q here bc we weren't given one. OK?
    P = P + Pdot*dt;
  
%   Assembling y_comp
    y_comp(:,i) = [vecnorm(x_(1:3,end)-cam1loc);
                   vecnorm(x_(1:3,end)-cam2loc);
                   vecnorm(x_(1:3,end)-cam3loc);
                   vecnorm(x_(1:3,end)-cam4loc);
                   x_(7,end);
                   x_(8,end);
                   x_(9,end)];

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
         -(X4 - drone_x)/((X4 - drone_x)^2 + (Y4 - drone_y)^2 + (Z4 - drone_z)^2)^(1/2), -(Y4 - drone_y)/((X4 - drone_x)^2 + (Y4 - drone_y)^2 + (Z4 - drone_z)^2)^(1/2), -(Z4 - drone_z)/((X4 - drone_x)^2 + (Y4 - drone_y)^2 + (Z4 - drone_z)^2)^(1/2), 0, 0, 0, 0, 0, 0, 0, 0, 0;
         0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0;
         0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;
         0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0]; 

%   Kalman gain
    K = P*H'*inv(H*P*H'+R);
    
    % Measurement update
    x_(:,i) = x_(:,i) + K * (y_obs(:,i) - y_comp(:,i));
    
    P = (eye(12)-K*H)*P*(eye(12)-K*H)' + K*R*K';
  
    store_x = [store_x x_(:, i)];
    store_P = [store_P norm(diag(P))];

    timeElapsed = toc;
    timeElapsedArray = [timeElapsedArray timeElapsed];
    
    fprintf("%g%% done\n", round(i/length(y1)*100))
end
figure(1)
hold on
scatter3(x_(1,:), x_(2,:), x_(3,:), 20, 'black')
legend('Observed', 'Predicted', 'Location', 'best')
hold off;
if shape == "hover"
    xlim([-1 1])
    ylim([-1 1])
end
saveas(gcf, strcat(shape, num2str(sample_freq), '_traj.png'))

err = store_x - states;

figure(2)
subplot(3,1,1)
hold on
plot(time, store_x(1,:), 'LineWidth', 1, 'Color', '#da7e30')
plot(time, states(1,:), '--', 'LineWidth', 1, 'Color', '#6b4c9a')
box on
grid on
ylabel("x position (m)")
subplot(3,1,2)
hold on
plot(time, store_x(2,:), 'LineWidth', 1, 'Color', '#da7e30')
plot(time, states(2,:), '--', 'LineWidth', 1, 'Color', '#6b4c9a')
box on
grid on
ylabel("y position (m)")
subplot(3,1,3)
hold on
plot(time, store_x(3,:), 'LineWidth', 1, 'Color', '#da7e30')
plot(time, states(3,:), '--', 'LineWidth', 1, 'Color', '#6b4c9a')
box on
grid on
ylabel("z position (m)")
xlabel("Time (s)")
legend("Estimated", "Actual", "Location", "best")
sgtitle("Comparison of Actual and Estimated states")
saveas(gcf, strcat(shape, num2str(sample_freq), '_compare.png'))

figure(3)
subplot(3,1,1)
hold on
plot(time, err(1,:), 'LineWidth', 1, 'Color', '#5ca793')
box on
grid on
ylabel("x position error (m)")
subplot(3,1,2)
hold on
plot(time, err(2,:), 'LineWidth', 1, 'Color', '#5ca793')
box on
grid on
ylabel("y position error (m)")
subplot(3,1,3)
hold on
plot(time, err(3,:), 'LineWidth', 1, 'Color', '#5ca793')
box on
grid on
ylabel("z position error (m)")
xlabel("Time (s)")
sgtitle("Error between Actual and Estimated states")
saveas(gcf, strcat(shape, num2str(sample_freq), '_err.png'))

writematrix([std(err(1,:));
             std(err(2,:));
             std(err(3,:));
             strcat(num2str(1/mean(time(2:end) - time(1:end-1))*1.2), "Hz");
             strcat(num2str(mean(timeElapsedArray)),"s per timestep")], strcat(shape, num2str(sample_freq), '_std.csv'))
end
end
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

function xhatdot = drone_dynamics(t, x, target_state, m, Ixx, Iyy, Izz, g)
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

end

function out = get_marker_locs(state)
    x = state(1);
    y = state(2);
    z = state(3);
    roll = state(7);
    pitch = state(8);
    yaw = state(9);
    marker_1_loc = [0; 0; 0.05];
    marker_2_loc = [0.105; 0; 0];
    marker_3_loc = [0; 0.06; 0];
    marker_4_loc = [-0.115; 0; 0];
    marker_5_loc = [0; -0.06; 0];
    R = [cos(pitch)*cos(yaw) cos(pitch)*sin(yaw) -sin(pitch);
         sin(roll)*sin(pitch)*cos(yaw)-cos(roll)*sin(yaw) sin(roll)*sin(pitch)*sin(yaw)+cos(roll)*cos(yaw) sin(roll)*cos(pitch);
         cos(roll)*sin(pitch)*cos(yaw)+sin(roll)*sin(yaw) cos(roll)*sin(pitch)*sin(yaw)-sin(roll)*cos(yaw) cos(roll)*cos(pitch)]';
     out = [state(1:3) + R*marker_1_loc;
           state(1:3) + R*marker_2_loc;
           state(1:3) + R*marker_3_loc;
           state(1:3) + R*marker_4_loc;
           state(1:3) + R*marker_5_loc];
end
