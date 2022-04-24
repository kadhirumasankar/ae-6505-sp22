clc; clear all; close all;

%% Prepping state data

m = 1.545;
Ixx = 0.029125;
Iyy = 0.029125;
Izz = 0.055225;
g = 9.81;

bag = rosbag('circle.bag');

statesMsgs = readMessages(select(bag, 'Topic', '/lqr_controller/states'),'DataFormat','struct');
states = [];
target_states = [];
marker_locs = [];
time = [];
for i = 1:length(statesMsgs)
    if mod(i,1)==0
        current_state = statesMsgs(i);
        current_state = current_state{1}.Data;
        states = [states get_states(current_state(1:12))];
        target_states = [target_states get_states(current_state(13:24))];
        marker_locs = [marker_locs get_marker_locs(current_state(1:12))];
        time = [time current_state(end)];
%         figure(1)
%         scatter3(current_state(1), current_state(2), current_state(3), 10, 'm', 'filled')
%         hold on
%         scatter3(marker_locs(1, end), marker_locs(2, end), marker_locs(3, end), 5, 'black', 'filled')
%         scatter3(marker_locs(4, end), marker_locs(5, end), marker_locs(6, end), 5, 'black', 'filled')
%         scatter3(marker_locs(7, end), marker_locs(8, end), marker_locs(9, end), 5, 'black', 'filled')
%         scatter3(marker_locs(10, end), marker_locs(11, end), marker_locs(12, end), 5, 'black', 'filled')
%         scatter3(marker_locs(13, end), marker_locs(14, end), marker_locs(15, end), 5, 'black', 'filled')
%         xlim([current_state(1)-.2 current_state(1)+.2])
%         ylim([current_state(2)-.2 current_state(2)+.2])
%         zlim([current_state(3)-.2 current_state(3)+.2])
%         hold off
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
y1 = mean([vecnorm(states(1:3,:)-cam1loc); vecnorm(states(4:6,:)-cam1loc); vecnorm(states(7:9,:)-cam1loc); vecnorm(states(10:12,:)-cam1loc)], 1);
y2 = mean([vecnorm(states(1:3,:)-cam2loc); vecnorm(states(4:6,:)-cam2loc); vecnorm(states(7:9,:)-cam2loc); vecnorm(states(10:12,:)-cam2loc)], 1);
y3 = mean([vecnorm(states(1:3,:)-cam3loc); vecnorm(states(4:6,:)-cam3loc); vecnorm(states(7:9,:)-cam3loc); vecnorm(states(10:12,:)-cam3loc)], 1);
y4 = mean([vecnorm(states(1:3,:)-cam4loc); vecnorm(states(4:6,:)-cam4loc); vecnorm(states(7:9,:)-cam4loc); vecnorm(states(10:12,:)-cam4loc)], 1);

%%
% Visualizing data
figure(1)
scatter3(mean(states(1:3:12,:),1), mean(states(2:3:12,:),1), mean(states(3:3:12,:),1), 10, 'm', 'filled')
% hold on
% scatter3(marker_locs(1,:), marker_locs(2,:), marker_locs(3,:), 5, 'r', 'filled')
% scatter3(marker_locs(4,:), marker_locs(5,:), marker_locs(6,:), 5, 'b', 'filled')
% scatter3(marker_locs(7,:), marker_locs(8,:), marker_locs(9,:), 5, 'g', 'filled')
% scatter3(marker_locs(10,:), marker_locs(11,:), marker_locs(12,:), 5, 'k', 'filled')
% scatter3(marker_locs(13,:), marker_locs(14,:), marker_locs(15,:), 5, 'k', 'filled')

% Initial conditions
x_(:,1) = states(:,1); % Using the first column from the data
P = zeros(12); % COMBAK: perfect knowledge of initial state so zero
Q = eye(12).*1e-7;
R = eye(7).*(.001)^2; % COMBAK: need to change this later to fit the function, add radian error to this

store_x = [x_(:,1)];
store_P = [norm(diag(P))];
1/(time(2)-time(1))
for i =2:length(y1)
    dt = time(i) - time(i-1);
    1/dt
%   Observed
    y_obs(:,i) = [y1(i); y2(i); y3(i); y4(i); states(7,i); states(8,i); states(9,i)] + sqrt(diag(R)).*randn(size(diag(R)));
  
%   Propagation of state
    [t_out, y_out] = ode45(@(t,y) drone_dynamics(t, y, target_states(:,i-1), m, Ixx, Iyy, Izz, g), [0 dt], x_(:, i-1), odeset('RelTol',1e-2,'AbsTol',1e-4));

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
%     H = [-(X1 - drone_x)/((X1 - drone_x)^2 + (Y1 - drone_y)^2 + (Z1 - drone_z)^2)^(1/2), -(Y1 - drone_y)/((X1 - drone_x)^2 + (Y1 - drone_y)^2 + (Z1 - drone_z)^2)^(1/2), -(Z1 - drone_z)/((X1 - drone_x)^2 + (Y1 - drone_y)^2 + (Z1 - drone_z)^2)^(1/2), 0, 0, 0, 0, 0, 0, 0, 0, 0;
%          -(X2 - drone_x)/((X2 - drone_x)^2 + (Y2 - drone_y)^2 + (Z2 - drone_z)^2)^(1/2), -(Y2 - drone_y)/((X2 - drone_x)^2 + (Y2 - drone_y)^2 + (Z2 - drone_z)^2)^(1/2), -(Z2 - drone_z)/((X2 - drone_x)^2 + (Y2 - drone_y)^2 + (Z2 - drone_z)^2)^(1/2), 0, 0, 0, 0, 0, 0, 0, 0, 0;
%          -(X3 - drone_x)/((X3 - drone_x)^2 + (Y3 - drone_y)^2 + (Z3 - drone_z)^2)^(1/2), -(Y3 - drone_y)/((X3 - drone_x)^2 + (Y3 - drone_y)^2 + (Z3 - drone_z)^2)^(1/2), -(Z3 - drone_z)/((X3 - drone_x)^2 + (Y3 - drone_y)^2 + (Z3 - drone_z)^2)^(1/2), 0, 0, 0, 0, 0, 0, 0, 0, 0;
%          -(X4 - drone_x)/((X4 - drone_x)^2 + (Y4 - drone_y)^2 + (Z4 - drone_z)^2)^(1/2), -(Y4 - drone_y)/((X4 - drone_x)^2 + (Y4 - drone_y)^2 + (Z4 - drone_z)^2)^(1/2), -(Z4 - drone_z)/((X4 - drone_x)^2 + (Y4 - drone_y)^2 + (Z4 - drone_z)^2)^(1/2), 0, 0, 0, 0, 0, 0, 0, 0, 0];
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

    figure(1)
    hold on
    scatter3(x_(1,end), x_(2,end), x_(3,end), 20, 'black')
    fprintf("%g%% done\n", (i/length(y1)*100))
    
    err = store_x(:,i) - states(:,i);
    
end
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
title('Trajectory of Baseball')
legend('Observed', 'Predicted', 'Location', 'best')
hold off;

err = store_x - states;

figure(2)
subplot(3,1,1)
hold on
plot(time, store_x(1,:))
plot(time, states(1,:))
subplot(3,1,2)
hold on
plot(time, store_x(2,:))
plot(time, states(2,:))
subplot(3,1,3)
hold on
plot(time, store_x(3,:))
plot(time, states(3,:))

writematrix(x_, 'baseball_ekf.csv')
%% Functions
function [u3 u4 u5 u6] = get_u(states, target_states, m, Ixx, Iyy, Izz, g)
    roll = states(7);
    pitch = states(8);
    yaw = states(9);
    roll_rate = states(10);
    pitch_rate = states(11);
    yaw_rate = states(12);
    u3 = m*g;
    A = find_A(states, m, Ixx, Iyy, Izz, g);
%     B = [0, 0, 0, 0;
%          0, 0, 0, 0;
%          0, 0, 0, 0;
%          (sin(q7) * sin(q9) + cos(q7) * cos(q9) * sin(q8)) / m, 0, 0, 0;
%          -(cos(q9) * sin(q7) - cos(q7) * sin(q8) * sin(q9)) / m, 0, 0, 0;
%          (cos(q7) * cos(q8)) / m, 0, 0, 0;
%          0, 0, 0, 0;
%          0, 0, 0, 0;
%          0, 0, 0, 0;
%          0, 1 / Ixx, 0, 0;
%          0, 0, 1 / Iyy, 0;
%          0, 0, 0, 1 / Izz];
    B = [                                             0,     0,     0,     0;
                                                      0,     0,     0,     0;
                                                      0,     0,     0,     0;
                                                      0,     0,     0,     0;
                                                      0,     0,     0,     0;
                                                      0,     0,     0,     0;
                                                      0,     0,     0,     0;
                                                      0,     0,     0,     0;
                                                      0,     0,     0,     0;
                                                      0,     0,     0,     0;
                                                      0,     0,     0,     0;
                                                      0,     0,     0,     0;
 (sin(roll)*sin(yaw) + cos(roll)*cos(yaw)*sin(pitch))/m,     0,     0,     0;
-(cos(yaw)*sin(roll) - cos(roll)*sin(pitch)*sin(yaw))/m,     0,     0,     0;
                               (cos(pitch)*cos(roll))/m,     0,     0,     0;
                                                      0,     0,     0,     0;
                                                      0,     0,     0,     0;
                                                      0,     0,     0,     0;
                                                      0, 1/Ixx,     0,     0;
                                                      0,     0, 1/Iyy,     0;
                                                      0,     0,     0, 1/Izz];
 
    Q = diag([100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 0.1 0.1 0.1 10 10 10]);
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
    roll = states(7);
    pitch = states(8);
    yaw = states(9);
    roll_rate = states(10);
    pitch_rate = states(11);
    yaw_rate = states(12);
    u3 = m*g;
%     A = [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0;
%          0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0;
%          0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0;
%          0, 0, 0, 0, 0, 0, (u3 * (cos(q7) * sin(q9) - cos(q9) * sin(q7) * sin(q8))) / m, ...
%                            (u3 * cos(q7) * cos(q8) * cos(q9)) / m, ...
%                            (u3 * (cos(q9) * sin(q7) - cos(q7) * sin(q8) * sin(q9))) / m, ...
%                            0, 0, 0;
%          0, 0, 0, 0, 0, 0, -(u3 * (cos(q7) * cos(q9) + sin(q7) * sin(q8) * sin(q9))) / m, ...
%                            (u3 * cos(q7) * cos(q8) * sin(q9)) / m, ...
%                            (u3 * (sin(q7) * sin(q9) + cos(q7) * cos(q9) * sin(q8))) / m, ...
%                            0, 0, 0;
%          0, 0, 0, 0, 0, 0, -(u3 * cos(q8) * sin(q7)) / m, ...
%                            -(u3 * cos(q7) * sin(q8)) / m, ...
%                            0, 0, 0, 0;
%          0, 0, 0, 0, 0, 0, q11 * cos(q7) * tan(q8) - q12 * sin(q7) * tan(q8), ...
%                            q12 * cos(q7) * (tan(q8)^2 + 1) + q11 * sin(q7) * (tan(q8)^2 + 1), ...
%                            0, 1, sin(q7) * tan(q8), cos(q7) * tan(q8);
%          0, 0, 0, 0, 0, 0, -q12 * cos(q7) - q11 * sin(q7), 0, 0, 0, cos(q7), -sin(q7);
%          0, 0, 0, 0, 0, 0, (q11 * cos(q7)) / cos(q8) - (q12 * sin(q7)) / cos(q8),  ...
%                            (q12 * cos(q7) * sin(q8)) / cos(q8)^2 + (q11 * sin(q7) * sin(q8)) / cos(q8)^2, ...
%                            0, 0, sin(q7) / cos(q8), cos(q7) / cos(q8);
%          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (Iyy * q12 - Izz * q12) / Ixx, ...
%                                        (Iyy * q11 - Izz * q11) / Ixx;
%          0, 0, 0, 0, 0, 0, 0, 0, 0, -(Ixx * q12 - Izz * q12) / Iyy, ...
%                                     0, -(Ixx * q10 - Izz * q10) / Iyy;
%          0, 0, 0, 0, 0, 0, 0, 0, 0, (Ixx * q11 - Iyy * q11) / Izz, ...
%                                     (Ixx * q10 - Iyy * q10) / Izz, 0];
    A = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,                                                                   0,                                                                                             0,                                                           0,                                     0,                                   0,                                     0;
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,                                                                   0,                                                                                             0,                                                           0,                                     0,                                   0,                                     0;
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,                                                                   0,                                                                                             0,                                                           0,                                     0,                                   0,                                     0;
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,                                                                   0,                                                                                             0,                                                           0,                                     0,                                   0,                                     0;
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,                                                                   0,                                                                                             0,                                                           0,                                     0,                                   0,                                     0;
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,                                                                   0,                                                                                             0,                                                           0,                                     0,                                   0,                                     0;
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,                                                                   0,                                                                                             0,                                                           0,                                     0,                                   0,                                     0;
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,                                                                   0,                                                                                             0,                                                           0,                                     0,                                   0,                                     0;
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,                                                                   0,                                                                                             0,                                                           0,                                     0,                                   0,                                     0;
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,                                                                   0,                                                                                             0,                                                           0,                                     0,                                   0,                                     0;
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,                                                                   0,                                                                                             0,                                                           0,                                     0,                                   0,                                     0;
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,                                                                   0,                                                                                             0,                                                           0,                                     0,                                   0,                                     0;
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         (u3*(cos(roll)*sin(yaw) - cos(yaw)*sin(pitch)*sin(roll)))/m,                                                          (u3*cos(pitch)*cos(roll)*cos(yaw))/m, (u3*(cos(yaw)*sin(roll) - cos(roll)*sin(pitch)*sin(yaw)))/m,                                     0,                                   0,                                     0;
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,        -(u3*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw)))/m,                                                          (u3*cos(pitch)*cos(roll)*sin(yaw))/m, (u3*(sin(roll)*sin(yaw) + cos(roll)*cos(yaw)*sin(pitch)))/m,                                     0,                                   0,                                     0;
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                                        -(u3*cos(pitch)*sin(roll))/m,                                                                  -(u3*cos(roll)*sin(pitch))/m,                                                           0,                                     0,                                   0,                                     0;
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,     pitch_rate*cos(roll)*tan(pitch) - yaw_rate*tan(pitch)*sin(roll),               yaw_rate*cos(roll)*(tan(pitch)^2 + 1) + pitch_rate*sin(roll)*(tan(pitch)^2 + 1),                                                           0,                                     1,                tan(pitch)*sin(roll),                  cos(roll)*tan(pitch);
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                         - yaw_rate*cos(roll) - pitch_rate*sin(roll),                                                                                             0,                                                           0,                                     0,                           cos(roll),                            -sin(roll);
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (pitch_rate*cos(roll))/cos(pitch) - (yaw_rate*sin(roll))/cos(pitch), (yaw_rate*cos(roll)*sin(pitch))/cos(pitch)^2 + (pitch_rate*sin(pitch)*sin(roll))/cos(pitch)^2,                                                           0,                                     0,                sin(roll)/cos(pitch),                  cos(roll)/cos(pitch);
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                                                                   0,                                                                                             0,                                                           0,                                     0,   (Iyy*yaw_rate - Izz*yaw_rate)/Ixx, (Iyy*pitch_rate - Izz*pitch_rate)/Ixx;
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                                                                   0,                                                                                             0,                                                           0,    -(Ixx*yaw_rate - Izz*yaw_rate)/Iyy,                                   0,  -(Ixx*roll_rate - Izz*roll_rate)/Iyy;
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                                                                   0,                                                                                             0,                                                           0, (Ixx*pitch_rate - Iyy*pitch_rate)/Izz, (Ixx*roll_rate - Iyy*roll_rate)/Izz,                                     0];
 
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
    x_pos = mean(x(1:3:12));
    y_pos = mean(x(2:3:12));
    z_pos = mean(x(3:3:12));
    x_vel = x(13);
    y_vel = x(14);
    z_vel = x(15);
    roll = x(16);
    pitch = x(17);
    yaw = x(18);
    roll_rate = x(19);
    pitch_rate = x(20);
    yaw_rate = x(21);
    [u3, u4, u5, u6] = get_u(x, target_state, m, Ixx, Iyy, Izz, g);

    xhatdot =[x_vel;
              y_vel;
              z_vel;
              x_vel;
              y_vel;
              z_vel;
              x_vel;
              y_vel;
              z_vel;
              x_vel;
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
    marker_1_loc = [0.105; 0; 0];
    marker_2_loc = [0; 0.06; 0];
    marker_3_loc = [-0.115; 0; 0];
    marker_4_loc = [0; -0.06; 0];
    R = [cos(pitch)*cos(yaw) cos(pitch)*sin(yaw) -sin(pitch);
         sin(roll)*sin(pitch)*cos(yaw)-cos(roll)*sin(yaw) sin(roll)*sin(pitch)*sin(yaw)+cos(roll)*cos(yaw) sin(roll)*cos(pitch);
         cos(roll)*sin(pitch)*cos(yaw)+sin(roll)*sin(yaw) cos(roll)*sin(pitch)*sin(yaw)-sin(roll)*cos(yaw) cos(roll)*cos(pitch)]';
     out = [state(1:3) + R*marker_1_loc;
            state(1:3) + R*marker_2_loc;
            state(1:3) + R*marker_3_loc;
            state(1:3) + R*marker_4_loc];
end

function out = get_states(state)
    marker_locs = get_marker_locs(state(1:12));
    out(1:12, :) = marker_locs;
    out = vertcat(out, state(4:end));
end