%%
clc; clear all; close all;

%%
data = readtable('Log_2021-07-06_163046.csv');

x = data.MotorOpticalSpeed_RPM_;
y = data.Thrust_N_;
x = x(y > 0.3044);
y = y(y > 0.3044);

p = polyfit(x.^2, y, 1);

x2 = linspace(min(x), max(x));
y2 = polyval(p, x2.^2);
scatter(x.^2, y, '.')
hold on;
plot(x2.^2, y2)
hold off;

%%
% Ct = 