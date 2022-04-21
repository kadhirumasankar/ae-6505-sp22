clc; clear all; close all;
rng(12)
tf = 60; % Simulation length
dt = 0.1; % Simulation step size
x_obs = [0 0 50 50]'; % Initial state
n = length(x_obs);
xhat = x_obs;
xhat_ukf = xhat;

% Tracking station location
N1 = 20;
E1 = 0;
N2 = 0; 
E2 = 20;

P = eye(4);
Pukf = P;
R = diag([1,1]);
Q = diag([0, 0, 4, 4]);
W = ones((2*n),1) / (2*n);

A = [1 0 dt 0;
     0 1 0 dt;
     0 0 1 0;
     0 0 0 1];

x_obs_list = x_obs;
xhatukf_list= xhat_ukf;
Pukf_list = diag(Pukf);

for t = dt : dt : tf
    % Simulate observations
    x_obs = A*x_obs;
    x_obs_list = [x_obs_list x_obs];
    x1 = x_obs(1);
    x2 = x_obs(2);
    H = [-(N1 - x1)/((E1 - x2)^2 + (N1 - x1)^2)^(1/2) -(E1 - x2)/((E1 - x2)^2 + (N1 - x1)^2)^(1/2) 0 0;
         -(N2 - x1)/((E2 - x2)^2 + (N2 - x1)^2)^(1/2) -(E2 - x2)/((E2 - x2)^2 + (N2 - x1)^2)^(1/2) 0 0];
    z_obs = H * x_obs + sqrt(diag(R)).*randn(2,1);

    % Generate the UKF sigma points.
    [root,p] = chol(n*Pukf);
    for i = 1 : n
       xhat_k_1(:,i) = xhat_ukf + root(i,:)';
       xhat_k_1(:,i+n) = xhat_ukf - root(i,:)';
    end
    % Time update
    for i = 1 : (2*n)
          xhat_k_1(:,i) = A*xhat_k_1(:,i);
    end
    xhat_ukf = zeros(n,1);
    for i = 1 : (2*n)
       xhat_ukf = xhat_ukf + W(i) * xhat_k_1(:,i);
    end
    Pukf = zeros(n,n);
    for i = 1 : (2*n)
       Pukf = Pukf + W(i) * (xhat_k_1(:,i) - xhat_ukf) * (xhat_k_1(:,i) - xhat_ukf)';
    end
    Pukf = Pukf + Q;

    % Measurement update
    for i = 1 : (2*n)
      zukf(:,i) = H*xhat_k_1(:,i);
    end
    zhat = 0;
    for i = 1 : (2*n)
      zhat = zhat + W(i) * zukf(:,i);
    end
    Py = 0;
    Pxy = zeros(n,1);
    for i = 1 : (2*n)
      Py = Py + W(i) * (zukf(:,i) - zhat) * (zukf(:,i) - zhat)';
      Pxy = Pxy + W(i) * (xhat_k_1(:,i) - xhat_ukf) * (zukf(:,i) - zhat)';
    end
    Py = Py + R;
    Kukf = Pxy * inv(Py);
    xhat_ukf = xhat_ukf + Kukf * (z_obs - zhat);
    Pukf = Pukf - Kukf * Py * Kukf';   
    xhatukf_list = [xhatukf_list xhat_ukf];
    Pukf_list = [Pukf_list diag(Pukf)];
    
end

% visualize
x_err = x_obs_list - xhatukf_list;
t = 0 : dt : tf;
figure()
subplot(2,2,1)
plot(t,x_err(1,:))
xlabel('Time (s)'); 
title('North Position Estimation Error');
subplot(2,2,2)
plot(t,x_err(2,:))
xlabel('Time (s)'); 
title('East Position Estimation Error');
subplot(2,2,3)
plot(t,x_err(3,:))
xlabel('Time (s)'); 
title('North Velocity Estimation Error');
subplot(2,2,4)
plot(t,x_err(4,:))
xlabel('Time (s)'); 
title('East Velocity Estimation Error');
