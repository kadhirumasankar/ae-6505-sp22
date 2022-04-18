clc; clear all; close all;

store_x = [];
dt = 0.1;
tf = 60;
MeasNoise = 1;
xdotNoise = [0,0,4,4];
Q = diag(xdotNoise);
P_k_1_plus = diag([1 1 1 1]); % IF PERFECTLY KNOWN IS IT ZERO COV OR 1 COV
R = diag([1 1]);

x_obs = [0; 0; 50; 50]; % initial state
xhat_k_1_plus = [0; 0; 50; 50]; % initial state estimate
N1 = 20;
E1 = 0;
N2 = 0;
E2 = 20;
y_obs(:, 1) = [sqrt((x_obs(1,1)-N1)^2 + (x_obs(2,1)-E1)^2);
               sqrt((x_obs(1,1)-N2)^2 + (x_obs(2,1)-E2)^2)];
F = [1 0 dt 0;
     0 1 0 dt;
     0 0 1 0;
     0 0 0 1];



T = 0.1; % measurement time step
tf = 60; % simulation length
dt = tf / 40000; % time step for integration
xArray = x_obs;
xhatArray = xhat_k_1_plus;
for t = T : T : tf
   % Simulate the system.
   F = [1 0 T 0;
        0 1 0 T;
        0 0 1 0;
        0 0 0 1];
   x_obs = F*x_obs;
   
   % Simulate the measurement.
   y_obs = [sqrt((x_obs(1)-N1)^2 + (x_obs(2)-E1)^2);
            sqrt((x_obs(1)-N2)^2 + (x_obs(2)-E2)^2)] + [MeasNoise*randn; MeasNoise*randn] / sqrt(dt);
   % Simulate the continuous-time part of the filter.
   F = [1 0 T 0;
        0 1 0 T;
        0 0 1 0;
        0 0 0 1];
   xhat_k_minus = F*xhat_k_1_plus;
   P_k_minus = F*P_k_1_plus*F' + Q;
   
   xhat_k_i_plus = xhat_k_minus;
   P_k_i_plus = P_k_minus;
   
   N = 100;
   for i=1:N
       x1 = xhat_k_i_plus(1);
       x2 = xhat_k_i_plus(2);
       H = [-(N1 - x1)/((E1 - x2)^2 + (N1 - x1)^2)^(1/2) -(E1 - x2)/((E1 - x2)^2 + (N1 - x1)^2)^(1/2) 0 0;
            -(N2 - x1)/((E2 - x2)^2 + (N2 - x1)^2)^(1/2) -(E2 - x2)/((E2 - x2)^2 + (N2 - x1)^2)^(1/2) 0 0];
       y_comp = [sqrt((x1-N1)^2 + (x2-E1)^2);
                 sqrt((x1-N2)^2 + (x2-E2)^2)];
       K = P_k_minus * H' * inv(H * P_k_minus * H' + R);
       P_k_ip1_plus = (eye(4) - K * H) * P_k_minus;
       xhat_k_ip1_plus = xhat_k_minus + K * (y_obs - y_comp);
       xhat_k_i_plus = xhat_k_ip1_plus;
       P_k_i_plus = P_k_ip1_plus;
   end
   
   % Save data for plotting.
   xArray = [xArray x_obs];
   xhatArray = [xhatArray xhat_k_i_plus];
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
