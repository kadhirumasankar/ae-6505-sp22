%% 3.5
clc; clear all; close all;
syms y1 y2
H = [1; 1];
W = [1 0; 0 1/4];
y = [y1; y2];
inv(H'*W*H)*H'*W*y

%% 3.5.b
clc; clear all; close all
P0 = [1];
H1 = [1];
R1 = [1];
K1 = P0*H1'*inv(H1*P0*H1'*R1)
P1 = (1 - K1*H1)*P0*(1-K1*H1)' + K1*R1*K1'

H2 = [1 1]';
R2 = [1 0;
      0 4]';
K2 = P1*H2'*inv(H2*P1*H2'*R2)

%%
clc; clear all; close all;
P(1) = [1];
H = [1; 1; 1];
R = [1; 1; 4];

for k=2:3
    K(k) = P(k-1)*H(k)'*inv((H(k)*P(k-1)*H(k)' + R(k)));
    P(k) = ((eye(length(K(k))) - K(k)*H(k)) * P(k-1) * (eye(length(K(k))) - K(k)*H(k))') + K(k)*R(k)*K(k)'
end

%% 3.7
clc; clear all; close all;
actual = [1 2 3];
estimate_1 = [3 4 1];
estimate_2 = [1 2 6];
estimate_3 = [5 6 7];
RMS_1 = sqrt(sum(estimate_1.^2)/3)
RMS_2 = sqrt(sum(estimate_2.^2)/3)
RMS_3 = sqrt(sum(estimate_3.^2)/3)
MAE_1 = sum(abs(estimate_1 - actual))/3;
MAE_2 = sum(abs(estimate_2 - actual))/3;
MAE_3 = sum(abs(estimate_3 - actual))/3;
std_1 = std(estimate_1)
std_2 = std(estimate_2)
std_3 = std(estimate_3)

%% 3.13
clc; clear all; close all;

years = 0:10;
data = [66.6 84.9 88.6 78.0 96.8 105.2 93.2 111.6 88.3 117.0 115.2];
Y = data';
X = [ones(length(data), 1) years'];
A = inv(X'*X)*X'*Y;

figure()
scatter(years, data, 'LineWidth', 1)
hold on
fplot(@(x) A(1) + A(2)*x, [years(1) years(end)], 'LineWidth', 1)
xlabel('Years since 1946')
ylabel('Steel produced (in millions of tons)')
legend('Original data', 'Linear curve fit', 'Location', 'NorthWest')
title('Linear curve fit')
hold off

RMS_1 = norm((A(1) + A(2)*years) - data)/sqrt(length(years))
pred = A(1) + A(2)*12

X = [ones(length(data), 1) years' (years.^2)'];
A = inv(X'*X)*X'*Y;

figure()
scatter(years, data, 'LineWidth', 1)
hold on
fplot(@(x) A(1) + A(2)*x + A(3)*x.^2, [years(1) years(end)], 'LineWidth', 1)
xlabel('Years since 1946')
ylabel('Steel produced (in millions of tons)')
legend('Original data', 'Quadratic curve fit', 'Location', 'NorthWest')
title('Quadratic curve fit')
hold off

RMS_2 = norm((A(1) + A(2)*years + A(3)*years.^2) - data)/sqrt(length(years))
pred = A(1) + A(2)*12 + A(3)*12.^2

X = [ones(length(data), 1) years' (years.^2)' (years.^3)'];
A = inv(X'*X)*X'*Y;

figure()
scatter(years, data, 'LineWidth', 1)
hold on
fplot(@(x) A(1) + A(2)*x + A(3)*x.^2 + A(4)*x.^3, [years(1) years(end)], 'LineWidth', 1)
xlabel('Years since 1946')
ylabel('Steel produced (in millions of tons)')
legend('Original data', 'Cubic curve fit', 'Location', 'NorthWest')
title('Cubic curve fit')
hold off

RMS_3 = norm((A(1) + A(2)*years + A(3)*years.^2 + A(4)*years.^3) - data)/sqrt(length(years))
pred = A(1) + A(2)*12 + A(3)*12.^2 + A(4)*12.^3

X = [ones(length(data), 1) years' (years.^2)' (years.^3)' (years.^4)'];
A = inv(X'*X)*X'*Y;

figure()
scatter(years, data, 'LineWidth', 1)
hold on
fplot(@(x) A(1) + A(2)*x + A(3)*x.^2 + A(4)*x.^3 + A(5)*x.^4, [years(1) years(end)], 'LineWidth', 1)
xlabel('Years since 1946')
ylabel('Steel produced (in millions of tons)')
legend('Original data', 'Quartic curve fit', 'Location', 'NorthWest')
title('Quartic curve fit')
hold off

RMS_4 = norm((A(1) + A(2)*years + A(3)*years.^2 + A(4)*years.^3 + A(5)*years.^4) - data)/sqrt(length(years))
pred = A(1) + A(2)*12 + A(3)*12.^2 + A(4)*12.^3 + A(5)*12.^4

%% 2
clc; clear all; close all;
data = load("emp_acc_data.txt");
scatter(data(:, 1), data(:, 2), '.')
Y = data(:,2);
X = [ones(length(data(:, 1)), 1) data(:, 1) cos(2*pi*data(:, 1)/5677) sin(2*pi*data(:, 1)/5677) cos(2*pi*data(:, 1)/2838.5) sin(2*pi*data(:, 1)/2838.5)];
A = inv(X'*X)*X'*Y;
hold on
fplot(@(x) A(1) + A(2)*x + A(3)*cos(2*pi*x/5677) + A(4)*sin(2*pi*x/5677) + A(5)*cos(2*pi*x/2838.5) + A(6)*sin(2*pi*x/2838.5), [data(1, 1) data(end, 1)], 'LineWidth', 1)
legend("Original data", "Final estimated curve")
hold off

perturbed = A(1) + A(2)*data(:, 1) + A(3)*cos(2*pi*data(:, 1)/5677) + A(4)*sin(2*pi*data(:, 1)/5677) + A(5)*cos(2*pi*data(:, 1)/2838.5) + A(6)*sin(2*pi*data(:, 1)/2838.5);

r = rms(data(:, 2))
s = rms(data(:, 2) - perturbed)
(r^2 - s^2)/r^2