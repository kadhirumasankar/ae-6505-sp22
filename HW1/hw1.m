%%
clc; clear all; close all;
A = [2 3; 3 4];
B = [4 2; 2 -4];
A*B

%% 1.8
clc; clear all; close all;
syms a b c L
A = [a b; b c];
det([a-L b; b c-L])
det([L-a b; b L-c])
aa = 1;
bb = -a-c;
cc = a*c - b^2;
(-bb + sqrt(bb^2 - 4*aa*cc))/(2*aa)

%% 1.17
clc; clear all; close all;

[T,X]=ode45(@fcneom,0:0.05:5,[0 0]);
figure()
plot(T, X(:, 1), '-o', 'LineWidth', 1)
hold on
plot(T, X(:, 2), '-o', 'LineWidth', 1)
legend("Angular Pos", "Angular Vel")
title("Δt = 0.05s")
hold off

[T,X]=ode45(@fcneom,0:0.2:5,[0 0]);
figure()
plot(T, X(:, 1), '-o', 'LineWidth', 1)
hold on
plot(T, X(:, 2), '-o', 'LineWidth', 1)
legend("Angular Pos", "Angular Vel")
title("Δt = 0.2s")
hold off

[T,X]=ode45(@fcneom,0:0.5:5,[0 0]);
figure()
plot(T, X(:, 1), '-o', 'LineWidth', 1)
hold on
plot(T, X(:, 2), '-o', 'LineWidth', 1)
legend("Angular Pos", "Angular Vel")
title("Δt = 0.5s")
hold off

%% 1.17 v2
clc; clear all; close all;

F = 100;
J = 10;
T = 10;
x1(1) = 0;
x2(1) = 0;
tspan = 0:0.05:5;
for i = 2:length(tspan)
    xdot = [0 1; 0 -F/J]*[x1(i-1); x2(i-1)] + [0 1/J]' * T;
    x1(i) = x1(i-1) + xdot(1)*0.05;
    x2(i) = x2(i-1) + xdot(2)*0.05;
end
figure()
plot(tspan, x1, '-o', 'LineWidth', 1)
hold on
plot(tspan, x2, '-o', 'LineWidth', 1)
legend("Angular Pos", "Angular Vel")
title("Δt = 0.05s")
hold off

x1 = [];
x2 = [];
x1(1) = 0;
x2(1) = 0;
tspan = 0:0.2:5;
for i = 2:length(tspan)
    xdot = [0 1; 0 -F/J]*[x1(i-1); x2(i-1)] + [0 1/J]' * T;
    x1(i) = x1(i-1) + xdot(1)*0.2;
    x2(i) = x2(i-1) + xdot(2)*0.2;
end
figure()
plot(tspan, x1, '-o', 'LineWidth', 1)
hold on
plot(tspan, x2, '-o', 'LineWidth', 1)
legend("Angular Pos", "Angular Vel")
title("Δt = 0.2s")
hold off

deltaT = 0.5;
x1 = [];
x2 = [];
x1(1) = 0;
x2(1) = 0;
tspan = 0:deltaT:5;
for i = 2:length(tspan)
    xdot = [0 1; 0 -F/J]*[x1(i-1); x2(i-1)] + [0 1/J]' * T;
    x1(i) = x1(i-1) + xdot(1)*deltaT;
    x2(i) = x2(i-1) + xdot(2)*deltaT;
end
figure()
plot(tspan, x1, '-o', 'LineWidth', 1)
hold on
plot(tspan, x2, '-o', 'LineWidth', 1)
legend("Angular Pos", "Angular Vel")
title("Δt = 0.5s")
hold off

%% 2.9
clc; clear all; close all;
x = linspace(-2, 2);
plot(x, x, 'o')
hold on
plot(x, x.^2, 'o')
grid on
axis on
legend('y = x', 'y = x^2')
hold off

%% q2.15
clc; clear all; close all;

x = rand(50, 1);
figure()
histogram(x, 10)
fprintf("N = 50\nµ = %f\nσ = %f\n\n", mean(x), std(x))

x = rand(500, 1);
figure()
histogram(x, 10)
fprintf("N = 500\nµ = %f\nσ = %f\n\n", mean(x), std(x))

x = rand(5000, 1);
figure()
histogram(x, 10)
fprintf("N = 5000\nµ = %f\nσ = %f\n\n", mean(x), std(x))

%% functions
function eom = fcneom(t, state)
    F = 100;
    J = 10;
    T = 10;
    xdot = [0 1; 0 -F/J]*state + [0 1/J]' * T;
    eom = xdot;
end