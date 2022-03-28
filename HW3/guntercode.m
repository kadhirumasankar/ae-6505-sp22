clear all; close all; 

n = 6; % number of params
k = 3; % number of obs per epoch
max_iter = 8;
g = 9.81; % m/s^2

dt = 0.1;
t_all = [0:dt:3.5];

% Define initial conditions
X_0 = [0.15 0.15 0.61 -8 40.23*cosd(30) 40.23*sind(30)]';
Y_obs = dlmread('homerun_data.txt', '\t');

t = Y_obs(:,1); 
Y_obs = [Y_obs(:,2)/3.28084 Y_obs(:,3)*pi/180 Y_obs(:,4)*pi/180];

m = size(Y_obs, 1);

sig1 = 1.524;
sig2 = pi/180/10;

% Visualize data
for i=1:12
    a = [Y_obs(i, 1) 0 0]';
    alpha = -Y_obs(i, 2);
    beta = Y_obs(i, 3);
    
    Ry = [cos(beta) 0 -sin(beta);
          0 1 0;
          sin(beta) 0 cos(beta)];
    Rz = [cos(alpha) sin(alpha) 0;
          -sin(alpha) cos(alpha) 0;
          0 0 1];
    b = Rz * Ry * a;
    scatter3(b(1), b(2), b(3), 10, 'm', 'filled');
    hold on;
end

% Run iterations
for i = 1:max_iter
    HtH = zeros(n, n);
    Hty = zeros(n, 1);
    
    % Loop over all observations
    % Compute the nominal trajectory and computed observations
    % Compute Normal equations
    for j=1:m
        Phi = [1 0 0 t(j) 0 0;
               0 1 0 0 t(j) 0;
               0 0 1 0 0 t(j);
               0 0 0 1 0 0;
               0 0 0 0 1 0;
               0 0 0 0 0 1];
        
        % Compute nominal trajectory
        X_nom = Phi * X_0 + [0; 0; -1/2*g*t(j)^2; 0; 0; -g*t(j)];
        
        X = X_nom(1);
        Y = X_nom(2);
        Z = X_nom(3);
        
        % Compute partials
        rho = sqrt(X^2 + Y^2 + Z^2);
        
        H11 = X/rho;
        H21 = (-Y/X^2)/(1 + (Y/X)^2);
        H31 = ((-X*Z)/(X^2 + Y^2)^(3/2))/(1+(Z^2)/(X^2 + Y^2));
        H12 = Y/rho;
        H22 = (1/X)/(1+(Y/X)^2);
        H32 = ((-Y*Z)/(X^2 + Y^2)^(3/2))/(1 + (Z^2)/(X^2 + Y^2));
        H13 = Z/rho;
        H23 = 0;
        H33 = ((1)/(X^2 + Y^2)^(1/2))/(1 + (Z^2)/(X^2 + Y^2));
        
        H_tilda = [H11 H12 H13 0 0 0;
                   H21 H22 H23 0 0 0;
                   H31 H32 H33 0 0 0];
               
        % Generate observation deviation (i.e. observed minus computed)
        rho = sqrt(X^2 + Y^2 + Z^2);
        alpha = atan2(Y, X);
        beta = atan2(Z, sqrt(X^2 + Y^2));
        
        Y_computed = [rho; alpha; beta];
        y = Y_obs(j,:)' - Y_computed;
        
        % Use state transition matrix to map partials back to initial
        % coords
        H = H_tilda * Phi;
        
        W = [1/sig1^2 0 0; 0 1/sig2^2 0; 0 0 1/sig2^2];
        
        % Compute rank-k update of normal equations for this observation
        HtH = HtH + H'*W*H;
        Hty = Hty + H'*W*y;
        
    end
    
    % Compute state deviation
    P = inv(HtH);
    x_hat = P * Hty;
    
    % Add state deviation to get final estimate
    X_new = X_0 + x_hat
    
    % Use the new estimates as the starting point for the next iteration
    X_0 = X_new;
    
    %Plot estimated trajectory for current iteration
    figure(1);
    hold on;
    for k=1:size(t_all, 2)
        X0 = X_0(1);
        Y0 = X_0(2);
        Z0 = X_0(3);
        Vx0 = X_0(4);
        Vy0 = X_0(5);
        Vz0 = X_0(6);
        X = X0 + Vx0 * t_all(k);
        Y = Y0 + Vy0 * t_all(k);
        Z = Z0 + Vz0 * t_all(k) - (1/2)*g*t_all(k)^2;
        scatter3(X, Y, Z, 20, [.9-.8*(i/max_iter) .9-.8*(i/max_iter) .9-.8*(i/max_iter)])
    end
    
end