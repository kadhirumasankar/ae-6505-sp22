clc; clear all; close all;

data = importdata('homerun_data_HW4.txt');
% data(1,:) = [];
tspan = data(:, 1);
rho = data(:, 2)/3.281; % in m
alpha = data(:, 3)*pi/180; % in rad
beta = data(:, 4)*pi/180; % in rad
g = 9.81;
dt = 0.1;
store_x = [];

% Visualizing data
figure(1)
[plotx ploty plotz] = sph2cart(alpha, beta, rho);
scatter3(plotx, ploty, plotz, 10, 'm', 'filled')

% Initial Conditions
% x_(:,1) = [0.4921 0.4921 2.0013 -26.2467 114.3051 65.9941]'/3.281;  
x_(:,1) = [-0.6562; 0.9843; 1.4764;  -32.8084; 119.6219; 55.7806]/3.281;  
% x_(:,1) = zeros(6,1);% Given initial conditions
P = diag([4 4 4 0.1 0.1 0.1])/3.281^2;
R = [(1.524)^2 0 0;
     0 (0.1*pi/180)^2 0;
     0 0 (0.1*pi/180)^2];      %the error covariance constant to be used
% Q = diag([0.5^2,0.5^2,0.5^2,0.01^2,0.01^2,0.01^2]);
Q = zeros(6);

N = 1000;

xpartArray = [];
PpartArray = [];

Phi = [1 0 0 dt 0 0;
       0 1 0 0 dt 0;
       0 0 1 0 0 dt;
       0 0 0 1 0 0;
       0 0 0 0 1 0;
       0 0 0 0 0 1];
     
for i = 2:length(rho)

    for j = 1 : N
        xpart(:,j) = x_(:,i-1) + sqrt(diag(P)) .* randn(size(diag(P)));
        xpartminus(:,j) = Phi*xpart(:,j) + [0; 0; -dt^2/2; 0; 0; -dt]*g + sqrt(diag(Q))* randn;
        
        X_ = xpartminus(1,j);
        Y_ = xpartminus(2,j);
        Z_ = xpartminus(3,j);
        rho_    = sqrt(X_^2 + Y_^2 + Z_^2);
        alpha_  = atan2(Y_,X_);
        beta_   = atan2( Z_, sqrt(X_^2+Y_^2) );
        y_comp = [rho_; alpha_; beta_];
        y_comp = y_comp + sqrt(diag(R)).*randn(size(diag(R)));
        
        vhat = [rho(i); alpha(i); beta(i)] - y_comp;
        q1(j) = (1 / sqrt(R(1,1)) / sqrt(2*pi)) * exp(-vhat(1)^2*R(1,1)^-1/2);
        q2(j) = (1 / sqrt(R(2,2)) / sqrt(2*pi)) * exp(-vhat(2)^2*R(2,2)^-1/2);
        q3(j) = (1 / sqrt(R(3,3)) / sqrt(2*pi)) * exp(-vhat(3)^2*R(3,3)^-1/2);
    end

    qsum1 = sum(q1);
    qsum2 = sum(q2);
    qsum3 = sum(q3);
    tolerance = 1e-99;
    if qsum1 > tolerance
        q1 = q1./qsum1;
    end
    if qsum2 > tolerance
        q2 = q2./qsum2;
    end
    if qsum3 > tolerance
        q3 = q3./qsum3;
    end
    q = mean([q1; q2; q3], 1);
%     for j=1:N
%         if qsum1 > tolerance
%             q1(j) = q1(j) / qsum1;
%         end
%         if qsum2 > tolerance
%             q2(j) = q2(j) / qsum2;
%         end
%         if qsum3 > tolerance
%             q3(j) = q3(j) / qsum3;
%         end
%         q(j) = mean([q1(j) q2(j) q3(j)]);
%     end
       
    for j = 1 : N
        u = rand;
        qtempsum = 0;
        for k = 1 : N
            qtempsum = qtempsum + q(k);
            if qtempsum >= u
                i;
                resampled_xpart(:,j) = xpartminus(:,j);
                break;
            end
        end
    end
    
    x_(:,i) = mean(resampled_xpart,2);
    P = diag(var(resampled_xpart,0,2));
    
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
writematrix(x_, 'baseball_pf.csv')