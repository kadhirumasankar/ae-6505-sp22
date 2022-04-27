close all; clear;

% simple script demonstating state space form and its ability to easily deal
% with complex differential equations, such as F=ma with complicated force
% input

% integration paramters
dt = .01
t = 0:dt:10
n = length(t)
x = zeros(2, n)

% force input, simple or complex
% u = ones(1, n)
u = rand(1, n) .* (sin(t) + t.^5 + 3) .* exp(-t)

% euler integration
% professional solvers usually use runge-kutta 4 or 5
for i = 1:(n-1)
    x(:, i+1) = x(:, i) + dynamics(x(:, i), u(i))*dt;
end

plot(t, [x; u])
xlabel("time")
legend("pos", "vel", "force")

% state space representation of F = ma
% not yet in x'=Ax+Bu form
function x_dot = dynamics(x, u)
mass = 1
x_dot = zeros(2, 1)

x_dot(1) = x(2);
x_dot(2) = u/mass;
end
