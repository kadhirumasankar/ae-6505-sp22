clc; clear all; close all;

for q = 0:0.01:0.1
    rng('default')
    rng(12);
    Ra = 1.9; % Winding resistance
    L = 0.003; % Winding inductance
    lambda = 0.1; % Motor constant
    J = 0.00018; % Moment of inertia
    B = 0.001; % Coefficient of viscous friction

    R = diag([0.1^2 0.1^2]);
    ControlNoise = q; % std dev of uncertainty in control inputs
    xdotNoise = [ControlNoise/L ControlNoise/L 0.5 0]; % THIS IS FROM THE EXAMPLE
    Q = diag([xdotNoise(1)^2 xdotNoise(2)^2 xdotNoise(3)^2 xdotNoise(4)^2]);
    P = 1*eye(4); % Initial state estimation covariance
    x = [0; 0; 0; 0]; % initial state
    xhat = [0; 0; 0; 0]; % initial state estimate
    w = 2 * pi; % Control input frequency
    H = [1 0 0 0; 0 1 0 0]; % measurement matrix

    T = 0.1; % measurement time step
    tf = 1.5; % simulation length
    dt = tf / 40000; % time step for integration
    xArray = x;
    xhatArray = xhat;
    Parray = P(3,3);
    for t = T : T : tf
       % Simulate the system.
       for tau = dt : dt : T
          ua0 = sin(w*tau);
          ub0 = cos(w*tau);
          xdot(1,1) = -Ra/L*x(1) + x(3)*lambda/L*sin(x(4)) + ua0/L;
          xdot(2,1) = -Ra/L*x(2) - x(3)*lambda/L*cos(x(4)) + ub0/L;
          xdot(3,1) = -3/2*lambda/J*x(1)*sin(x(4)) + 3/2*lambda/J*x(2)*cos(x(4)) - B/J*x(3);
          xdot(4,1) = x(3);
          x = x + xdot * dt;
       end
       % Simulate the measurement.
       z = H * x + [sqrt(R(1,1)); sqrt(R(2,2))] * randn;
       % Simulate the continuous-time part of the filter.
       for tau = dt : dt : T
          ua0 = sin(w*tau);
          ub0 = cos(w*tau);
          xhatdot(1,1) = -Ra/L*xhat(1) + xhat(3)*lambda/L*sin(xhat(4)) + ua0/L;
          xhatdot(2,1) = -Ra/L*xhat(2) - xhat(3)*lambda/L*cos(xhat(4)) + ub0/L;
          xhatdot(3,1) = -3/2*lambda/J*xhat(1)*sin(xhat(4)) + 3/2*lambda/J*xhat(2)*cos(xhat(4)) - B/J*xhat(3);
          xhatdot(4,1) = xhat(3);
          xhat = xhat + xhatdot * dt;
          F = [-Ra/L 0 lambda/L*sin(xhat(4)) xhat(3)*lambda/L*cos(xhat(4));
               0 -Ra/L -lambda/L*cos(xhat(4)) xhat(3)*lambda/L*sin(xhat(4));
               -3/2*lambda/J*sin(xhat(4)) 3/2*lambda/J*cos(xhat(4)) -B/J -3/2*lambda/J*(xhat(1)*cos(xhat(4))+xhat(2)*sin(xhat(4)));
               0 0 1 0];
          Pdot = F * P + P * F' + Q;
          P = P + Pdot * dt;
       end
       % Simulate the discrete-time part of the filter.
       K = P * H' * inv(H * P * H' + R);
       xhat = xhat + K * (z - H * xhat);
       P = (eye(4) - K * H) * P * (eye(4) - K * H)' + K * R * K';
       % Save data for plotting.
       xArray = [xArray x];
       xhatArray = [xhatArray xhat];
       Parray = [Parray P(3,3)];
    end

    % Plot data
%     close all;
    t = 0 : T : tf;

    figure(1)
    hold on
    plot(t, Parray)

    N = size(xArray, 2);
    N2 = round(N / 2);
    N2 = 1;
    xArray = xArray(:,N2:N);
    xhatArray = xhatArray(:,N2:N);
%     wEstErr = sqrt(norm(xArray(3,:)-xhatArray(3,:))^2 / size(xArray,2))
    fprintf("Std Dev of Estimation Error of motor velocity for sigma_p = %G = %G\n", q, std(xArray(3,:)-xhatArray(3,:)))
%     std(xArray(3,:)-xhatArray(3,:))
end
title('Motor Velocity Estimation Error vs Time for different values of \sigma_p')
legend('\sigma_p = 0', '\sigma_p = 0.01', '\sigma_p = 0.02', '\sigma_p = 0.03', '\sigma_p = 0.04', '\sigma_p = 0.05', '\sigma_p = 0.06', '\sigma_p = 0.07', '\sigma_p = 0.08', '\sigma_p = 0.09', '\sigma_p = 0.1', 'Location', 'best')
% figure;
% plot(t, xArray(1,:) - xhatArray(1,:));
% set(gca,'FontSize',12); set(gcf,'Color','White');
% xlabel('Time (seconds)');
% ylabel('Altitude Estimation Error (feet)');
% 
% AltErr = std(xArray(1,:) - xhatArray(1,:));
% disp(['Hybrid EKF RMS altitude estimation error = ', num2str(AltErr)]);
% 
% figure;
% plot(t, xArray(2,:) - xhatArray(2,:));
% set(gca,'FontSize',12); set(gcf,'Color','White');
% xlabel('Time (seconds)');
% ylabel('Velocity Estimation Error (feet/sec)');
% 
% VelErr = std(xArray(2,:) - xhatArray(2,:));
% disp(['Hybrid EKF RMS velocity estimation error = ', num2str(VelErr)]);
% 
% figure;
% plot(t, xArray(3,:) - xhatArray(3,:));
% set(gca,'FontSize',12); set(gcf,'Color','White');
% xlabel('Time (seconds)');
% ylabel('Ballistic Coefficient Estimation Error');
% 
% BallErr = std(xArray(3,:) - xhatArray(3,:));
% disp(['Hybrid EKF RMS ballistic coefficient estimation error = ', num2str(BallErr)]);
% 
% figure;
% plot(t, xArray(1,:)/1000);
% set(gca,'FontSize',12); set(gcf,'Color','White');
% xlabel('Time (seconds)');
% ylabel('True Altitude (thousands of feet)');
% 
% figure;
% plot(t, xArray(2,:));
% set(gca,'FontSize',12); set(gcf,'Color','White');
% xlabel('Time (seconds)');
% ylabel('True Velocity (feet/sec)');