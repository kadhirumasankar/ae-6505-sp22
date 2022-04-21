clc; clear all; close all;
rng(12)
store_x = [];
dt = 0.1;
tf = 60;
MeasNoise = 1;
xdotNoise = [0,0,4,4];
Q = diag(xdotNoise);
P = diag([1 1 1 1]); % IF PERFECTLY KNOWN IS IT ZERO COV OR 1 COV
R = diag([1 1]);

x_obs(:, 1) = [0; 0; 50; 50];
N1 = 20;
E1 = 0;
N2 = 0;
E2 = 20;
y_obs(:, 1) = [sqrt((x_obs(1,1)-N1)^2 + (x_obs(2,1)-E1)^2);
               sqrt((x_obs(1,1)-N2)^2 + (x_obs(2,1)-E2)^2)];
A = [1 0 dt 0;
     0 1 0 dt;
     0 0 1 0;
     0 0 0 1];

x_obs = [0; 0; 50; 50]; % initial state
xhat = [0; 0; 50; 50]; % initial state estimate

T = 0.1; % measurement time step
tf = 60; % simulation length
dt = tf / 40000; % time step for integration
xArray = x_obs;
xhatArray = xhat;
for t = T : T : tf
   % Simulate the system.
   A = [1 0 T 0;
        0 1 0 T;
        0 0 1 0;
        0 0 0 1];
   x_obs = A*x_obs;
   
   % Simulate the measurement.
   y_obs = [sqrt((x_obs(1)-N1)^2 + (x_obs(2)-E1)^2);
            sqrt((x_obs(1)-N2)^2 + (x_obs(2)-E2)^2)] + sqrt(diag(R)).*randn(2,1);
   % Simulate the continuous-time part of the filter.
   A = [1 0 T 0;
        0 1 0 T;
        0 0 1 0;
        0 0 0 1];
   xhat = A*xhat;
   P = A*P*A' + Q;
   
   x1 = xhat(1);
   x2 = xhat(2);
   H = [-(N1 - x1)/((E1 - x2)^2 + (N1 - x1)^2)^(1/2) -(E1 - x2)/((E1 - x2)^2 + (N1 - x1)^2)^(1/2) 0 0;
        -(N2 - x1)/((E2 - x2)^2 + (N2 - x1)^2)^(1/2) -(E2 - x2)/((E2 - x2)^2 + (N2 - x1)^2)^(1/2) 0 0];
   y_comp = [sqrt((x1-N1)^2 + (x2-E1)^2);
             sqrt((x1-N2)^2 + (x2-E2)^2)];
   K = P * H' * inv(H * P * H' + R);
   xhat = xhat + K * (y_obs - y_comp);
   P = (eye(4) - K * H) * P;
   % Save data for plotting.
   xArray = [xArray x_obs];
   xhatArray = [xhatArray xhat];
%    Parray = [Parray P(3,3)];
end

% Plot data
%     close all;
t = 0 : T : tf;

err = xArray - xhatArray;

figure(1)
subplot(2,2,1)
hold on
plot(t, err(1,:))
title('North Position Estimation Error')
subplot(2,2,2)
hold on
plot(t, err(2,:))
title('East Position Estimation Error')
subplot(2,2,3)
hold on
plot(t, err(3,:))
title('North Velocity Estimation Error')
subplot(2,2,4)
hold on
plot(t, err(4,:))
title('East Velocity Estimation Error')

figure(2)
hold on
plot(t, xArray(1,:))
plot(t, xhatArray(1,:))
legend('Obs', 'Comp')
