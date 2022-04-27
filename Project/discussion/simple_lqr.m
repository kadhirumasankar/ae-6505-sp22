close all; clear;

% simple script showing lqr formulation to control F = ma

% integration parameters
dt = .01
t = 0:dt:10
n = length(t)
x = zeros(2, n)
% initial conditions to be regulated
x(:, 1) = [12; -24]; 
% place holder for thrust and torque inputs
u = zeros(1, n)
% solve for optimal K matrix
K = lqr_impl()

% euler integration
% professional solvers usually use runge-kutta 4 or 5
for i = 1:(n-1)
    u(i) = -K*x(:, i);
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

function K = lqr_impl()
% state dependent dynamics
A = [0, 1;
     0, 0];
% actuator dependent dynamics
B = [0;
     1];

% Q matrix gives weights to importance of regulating each state in the cost
% function
Q = [1, 0;
     0, 1];

% R matrix gives weights to importance of regulating control effort in the
% cost function
% primarily to get control effort within actuator limits
R = [1];

% matlab solves algebraic ricatti equation for optimal K matrix
K = lqr(A, B, Q, R)
end