%%
clc; clear all; close all;
maxIter = 10;
foodW = normrnd(0,sqrt(10),[1,maxIter]);
obsW = normrnd(0,sqrt(10),[1,maxIter]); 
Xactual = [650; 250];
F = [1/2 2; 0 1];
Pplus(:,:,1) = [500 0; 0 200];
Q = [0 0; 0 10];
H = [1 0];
R = 10;
Xplus(:,1)= [600; 200];
n = 2;

y(:,1) = Xactual(1,1) + obsW(1);
std1 = [];
std2 = [];
K1 = [];
K2 = [];

% discrete kalman filter + simulation
for k = 2:maxIter
     Xactual(:,k) = F*Xactual(:,k-1) + [0 foodW(k)]';
     y(:,k) = Xactual(1,k) + obsW(k);
     
     Pminus(:,:,k) = F*Pplus(:,:,k-1)*F' + Q;
     K(:,:,k) = Pminus(:,:,k)*H'/(H*Pminus(:,:,k)*H' + R);
     Xminus(:,k) = F*Xplus(:,k-1);
     Xplus(:,k) = Xminus(:,k) + K(:,:,k)*(y(:,k)-H*Xminus(:,k)); 
     Pplus(:,:,k) = (eye(n) - K(:,:,k)*H)*Pminus(:,:,k);
     std1 = vertcat(std1, sqrt(Pplus(1,1,k)));
     std2 = vertcat(std2, sqrt(Pplus(2,2,k)));
     K1 = vertcat(K1, K(1,1,k));
     K2 = vertcat(K2, K(2,1,k));
end
figure()
subplot(2,2,1)
plot(Xactual(1,:))
hold on
plot(Xplus(1,:))
xlabel('Time')
ylabel('Wombats')
legend('Actual', 'Estimated', 'Location', 'best')
title('Population')

subplot(2,2,2)
plot(Xactual(2,:))
hold on
plot(Xplus(2,:))
xlabel('Time')
ylabel('Food')
legend('Actual', 'Estimated', 'Location', 'best')
title('Food Supply')

subplot(2,2,3)
plot(std1)
hold on
plot(std2)
xlabel('Time')
legend('Population', 'Food', 'Location', 'best')
title('Std Dev')

subplot(2,2,4)
plot(K1)
hold on
plot(K2)
xlabel('Time')
legend('K1', 'K2', 'Location', 'best')
title('Gain Matrix Elements')
