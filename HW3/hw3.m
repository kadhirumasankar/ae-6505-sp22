%% 4.5
clc; clear all; close all;
syms tk tau tk1
% tk = 1e-5;
% tk1 = 0;
A = [-1 0; 0 -2]*(tk-tau)
eye(2)*(tk-tk1)
% integral(@(tau)(expm(A)*eye(2)*expm(A)), tk1, tk)
B = [1 1; 1 2];
simplify(expm(A)*B*eye(2)*B')
simplify(expm(A)*B*eye(2)*B'*expm(A))
%% 4.8
clc; clear all; close all;
F = eye(2)*1/2;
Q = [1 1; 0 1];
P = 0;
for i=0:100
    i
    F^i*Q*F'^i
    P = P + F^i*Q*F'^i
end

%% 4.11
clc; clear all; close all;

m0 = 1;
P0 = 2;
f = -0.5;
qc = 1;
tf = 5;
deltat = 0.01;
Q = [];
xbar = [m0];
P = [P0];

for k = deltat:deltat:tf
    Q(int16((k-deltat)/deltat + 1)) = qc/(2*f)*(exp(2*f*deltat) - 1);
    xbar(int16(k/deltat)+1) = exp(floor(k/deltat)*f*deltat)*xbar(1);
    P(int16(k/deltat)+1) = (2*f*P(int16(k/deltat)) + qc)*deltat + P(int16(k/deltat));
end

figure()
subplot(2,1,1)
plot(0:deltat:tf, xbar, 'LineWidth', 1.5)
hold on
sgtitle('$\bar{x}$ and $P$ when $P_0 = 2$', 'Interpreter', 'latex')
ylabel('$\bar{x}$', 'Interpreter', 'latex')
xlabel('Time (s)')
subplot(2,1,2)
plot([0 deltat:deltat:tf], P, 'LineWidth', 1.5)
ylabel('$P$', 'Interpreter', 'latex')
xlabel('Time (s)')
hold off

P0 = 0;
Q = [];
xbar = [m0];
P = [P0]

for k = deltat:deltat:tf
    Q(int16((k-deltat)/deltat + 1)) = qc/(2*f)*(exp(2*f*deltat) - 1);
    xbar(int16(k/deltat)+1) = exp(floor(k/deltat)*f*deltat)*xbar(1);
    P(int16(k/deltat)+1) = (2*f*P(int16(k/deltat)) + qc)*deltat + P(int16(k/deltat));
end

figure()
subplot(2,1,1)
plot(0:deltat:tf, xbar, 'LineWidth', 1.5)
hold on
sgtitle('$\bar{x}$ and $P$ when $P_0 = 0$', 'Interpreter', 'latex')
ylabel('$\bar{x}$', 'Interpreter', 'latex')
xlabel('Time (s)')
subplot(2,1,2)
plot([0 deltat:deltat:tf], P, 'LineWidth', 1.5)
ylabel('$P$', 'Interpreter', 'latex')
xlabel('Time (s)')
hold off

%% 4.13
clc; clear all; close all;
syms T
AHat = [-1 0; 0 -2];
Q = [-0.8944 -0.7071; 0.4472 0.7071];
A = [0 2; -1 -3];
Q*expm(AHat*T)*inv(Q)

%%
clc; clear all; close all;
A = [0 2; -1 -3];
syms p11 p12 p21 p22
P = [p11 p12; p21 p22];
A*P + P*A'
solve(A*P + P*A' + [0 0; 0 1] == 0)

%%
clc; clear all; close all;
syms P1(t) P2(t)
A = [0 2; -1 -3];
Qc = eye(2);
P = [P1 0;
     0 P2];
dsolve(diff(P, t) == A*P + P*A' + Qc);
%%
clc; clear all; close all;

% Given info
h = 6.4;
k1 = 2.6;
k2 = 3.5;
m = 1.8;
w = sqrt((k1+k2)/m);
% Measurement noise covariance
R0 = [0.0625 0;
     0 0.01];
% Initial state
X0 = [3.9;
      0.5];
tspan = 0:10;
P0_bar = [1000 0;
          0 100];
X0_bar = [0; 0];

% Initial measurements
Y_obs = [];
rho = [6.3118 6.0217 5.1295 6.3685 5.5422 5.5095 5.9740 5.4930 6.8653 6.6661 5.0692]';
rhoDot = [0.3035 1.3854 -1.5512 0.6064 0.8642 -1.5736 1.1588 0.4580 -1.2328 1.4348 -0.4229]';
for i = 1:length(rho)
    Y_obs = vertcat(Y_obs, rho(i));
    Y_obs = vertcat(Y_obs, rhoDot(i));
end

x_hist = [X0(1)];
v_hist = [X0(2)];
for iter = 1:50
    [~, X_computed] = ode45(@spring_mass_dynamics, tspan, X0);
    
    Y_computed = [];
    for i = 1:length(X_computed)
        c_i = X_computed(i, 1);
        d_i = X_computed(i, 2);
        rho_i = sqrt(c_i^2 + h^2);
        Y_computed = vertcat(Y_computed, [rho_i; c_i*d_i/rho_i]);
    end

    y_dev = Y_obs - Y_computed;
    
    H = [];
    R = [];
    for t = 0:10
        Phi = [cos(w*t) 1/w*sin(w*t);
               -w*sin(w*t) cos(w*t)];
        X_i = Phi * X0;
        c_i = X_i(1);
        d_i = X_i(2);
        rho_i = sqrt(c_i^2 + h^2);
        Htilda = [c_i/rho_i 0;
                  d_i/rho_i - c_i^2*d_i/rho_i^3 c_i/rho_i];
        H = vertcat(H, Htilda * Phi);
        R = blkdiag(R, R0);
    end
    X0_hat = inv(H'*inv(R)*H + inv(P0_bar)) * (H'*inv(R)*y_dev + inv(P0_bar)*X0_bar);
    P0_hat = inv(H'*inv(R)*H + inv(P0_bar));
    
    X0 = X0 + X0_hat;
    P0_bar = P0_hat
    x_hist = vertcat(x_hist, X0(1));
    v_hist = vertcat(v_hist, X0(2));
end

plot(0:length(x_hist)-1, x_hist, 'LineWidth', 1.5)
hold on
plot(0:length(v_hist)-1, v_hist, 'LineWidth', 1.5)
title('X_0 Estimate over Multiple Iterations')
legend('x_0 (in m)', 'v_0 (in m/s)')
xlabel('Iterations')
grid on
box on
hold off

%% 2 math
clc; clear all; close all;
syms phi11(t) phi12(t) phi21(t) phi22(t) w real
phi_mat = [phi11 phi12; phi21 phi22];
conds = [phi11(0) == 1; phi12(0) == 0; phi21(0) == 0; phi22(0) == 1];
A = [0 1; -w^2 0];
sol = dsolve(diff(phi_mat, t) == A*phi_mat, conds)

%% 3 math
clc; clear all; close all;
syms x y z vx vy vz g rho alpha beta real
X = [x; y; z; vx; vy; vz];
F = [vx; vy; vz; 0; 0; -g];
A = jacobian(F, X)
rho = sqrt(x^2 + y^2 + z^2);
alpha = atan(y/x);
beta = atan(z/sqrt(x^2 + y^2));
G = [rho; alpha; beta];
Htilda = jacobian(G, X)

syms phi11(t) phi12(t) phi13(t) phi14(t) phi15(t) phi16(t)
syms phi21(t) phi22(t) phi23(t) phi24(t) phi25(t) phi26(t)
syms phi31(t) phi32(t) phi33(t) phi34(t) phi35(t) phi36(t)
syms phi41(t) phi42(t) phi43(t) phi44(t) phi45(t) phi46(t)
syms phi51(t) phi52(t) phi53(t) phi54(t) phi55(t) phi56(t)
syms phi61(t) phi62(t) phi63(t) phi64(t) phi65(t) phi66(t)

phi_mat = [phi11, phi12, phi13, phi14, phi15, phi16;
           phi21, phi22, phi23, phi24, phi25, phi26;
           phi31, phi32, phi33, phi34, phi35, phi36;
           phi41, phi42, phi43, phi44, phi45, phi46;
           phi51, phi52, phi53, phi54, phi55, phi56;
           phi61, phi62, phi63, phi64, phi65, phi66];
conds = [phi11(0) == 1;
         phi12(0) == 0;
         phi13(0) == 0;
         phi14(0) == 0;
         phi15(0) == 0;
         phi16(0) == 0;
         phi21(0) == 0;
         phi22(0) == 1;
         phi23(0) == 0;
         phi24(0) == 0;
         phi25(0) == 0;
         phi26(0) == 0;
         phi31(0) == 0;
         phi32(0) == 0;
         phi33(0) == 1;
         phi34(0) == 0;
         phi35(0) == 0;
         phi36(0) == 0;
         phi41(0) == 0;
         phi42(0) == 0;
         phi43(0) == 0;
         phi44(0) == 1;
         phi45(0) == 0;
         phi46(0) == 0;
         phi51(0) == 0;
         phi52(0) == 0;
         phi53(0) == 0;
         phi54(0) == 0;
         phi55(0) == 1;
         phi56(0) == 0;
         phi61(0) == 0;
         phi62(0) == 0;
         phi63(0) == 0;
         phi64(0) == 0;
         phi65(0) == 0;
         phi66(0) == 1];
sol = dsolve(diff(phi_mat, t) == A*phi_mat, conds);

sol_phi_mat = [sol.phi11, sol.phi12, sol.phi13, sol.phi14, sol.phi15, sol.phi16;
               sol.phi21, sol.phi22, sol.phi23, sol.phi24, sol.phi25, sol.phi26;
               sol.phi31, sol.phi32, sol.phi33, sol.phi34, sol.phi35, sol.phi36;
               sol.phi41, sol.phi42, sol.phi43, sol.phi44, sol.phi45, sol.phi46;
               sol.phi51, sol.phi52, sol.phi53, sol.phi54, sol.phi55, sol.phi56;
               sol.phi61, sol.phi62, sol.phi63, sol.phi64, sol.phi65, sol.phi66]
           
%% 3
clc; clear all; close all;

% Given info

% Measurement noise covariance
R0 = [5^2 0 0;
      0 (0.1*pi/180)^2 0;
      0 0 (0.1*pi/180)^2];
% Initial state
X0 = [0.4921; 0.4921; 2.0013; -26.2467; 114.3051; 65.9941];
P0_bar = [inf 0 0 0 0 0;
          0 inf 0 0 0 0;
          0 0 inf 0 0 0;
          0 0 0 inf 0 0;
          0 0 0 0 inf 0;
          0 0 0 0 0 inf];
X0_bar = [0; 0; 0; 0; 0; 0];

% Initial measurements
data = importdata('homerun_data.txt');
tspan = data(:, 1);
rho = data(:, 2);
alpha = data(:, 3)*pi/180;
beta = data(:, 4)*pi/180;

Y_obs = [];
for i = 1:length(rho)
    Y_obs = vertcat(Y_obs, rho(i));
    Y_obs = vertcat(Y_obs, alpha(i));
    Y_obs = vertcat(Y_obs, beta(i));
end

x_hist = [X0(1)];
y_hist = [X0(2)];
z_hist = [X0(3)];
all_X = [];
all_Y = [];
all_Z = [];
for iter = 1:10
    [~, X_computed] = ode45(@baseball_dynamics, tspan, X0);
    
    
    for i = 1:length(X_computed)
        x_i = X_computed(i, 1);
        y_i = X_computed(i, 2);
        z_i = X_computed(i, 3);
    end

    
    Y_computed = [];
    X_computed = [];

    
    H = [];
    R = [];
    for t = 0:0.1:1.1
        Phi = [1 0 0 t 0 0;
               0 1 0 0 t 0;
               0 0 1 0 0 t;
               0 0 0 1 0 0;
               0 0 0 0 1 0;
               0 0 0 0 0 1];
        B = [0; 0; -0.5*t^2; 0; 0; -t];
        X_i = Phi * X0 + B * 32.185;
        x = X_i(1);
        y = X_i(2);
        z = X_i(3);
        Htilda = [x/sqrt(x^2+y^2+z^2), y/sqrt(x^2+y^2+z^2), z/sqrt(x^2+y^2+z^2), 0, 0, 0;
                  -y/(x^2*(y^2/x^2 + 1)), 1/(x*(y^2/x^2 + 1)), 0, 0, 0, 0;
                  -(x*z)/((z^2/(x^2 + y^2) + 1)*(x^2 + y^2)^(3/2)), -(y*z)/((z^2/(x^2 + y^2) + 1)*(x^2 + y^2)^(3/2)), 1/((z^2/(x^2 + y^2) + 1)*(x^2 + y^2)^(1/2)), 0, 0, 0];
        H = vertcat(H, Htilda * Phi);
        R = blkdiag(R, R0);
        Y_computed = vertcat(Y_computed, [sqrt(x^2 + y^2 + z^2); atan2(y,x); atan2(z,sqrt(x^2 + y^2))]);
%         X_computed = vertcat(X_computed, 
        [temp_Y temp_X] = ode45(@baseball_dynamics, linspace(0,2), X0);
        all_X = horzcat(all_X, temp_X(:, 1));
        all_Y = horzcat(all_Y, temp_X(:, 2));
        all_Z = horzcat(all_Z, temp_X(:, 3));
    end
    y_dev = Y_obs - Y_computed;
    X0_hat = inv(H'*inv(R)*H + inv(P0_bar)) * (H'*inv(R)*y_dev + inv(P0_bar)*X0_bar);
%     X0_hat = inv(H'*inv(R)*H) * (H'*inv(R)*y_dev);
%     P0_hat = inv(H'*inv(R)*H);
    P0_hat = inv(H'*inv(R)*H + inv(P0_bar));
    
    X0 = X0 + X0_hat;
    P0_bar = P0_hat;
    x_hist = vertcat(x_hist, X0(1));
    y_hist = vertcat(y_hist, X0(2));
    z_hist = vertcat(z_hist, X0(3));
end

figure()
plot(x_hist)
hold on
plot(y_hist)
plot(z_hist)
hold off

figure()
[x y z] = sph2cart(alpha, beta, rho)
scatter3(x, y, z)
hold on
% for i=1:size(all_X, 2)
%     plot3(all_X(:, i), all_Y(:, i), all_Z(:, i))
% end
plot3(all_X(:, end), all_Y(:, end), all_Z(:, end))
xlabel('x (ft)')
ylabel('y (ft)')
zlabel('z (ft)')
title('Trajectory of Baseball')
% plot3(all_X(:, 2), all_Y(:, 2), all_Z(:, 2))
% plot3(all_X(:, 3), all_Y(:, 3), all_Z(:, 3))
% plot3(all_X(:, 6), all_Y(:, 4), all_Z(:, 4))
% legend('s', '1', '2', '3', '4', '5')
hold off

[~, idx] = min(all_Z(:, end).^2);
sqrt(all_X(idx, end)^2 + all_Y(idx, end)^2 + all_Z(idx, end)^2)
norm(P0_bar, 'fro')

%% Functions
function Xdot = spring_mass_dynamics(t, state)
    k1 = 2.6;
    k2 = 3.5;
    m = 1.8;
    w = sqrt((k1+k2)/m);
    Xdot = [state(2);
            -w^2*state(1)];
end

function Xdot = baseball_dynamics(t, state)
    g = 32.185;
    Xdot = [state(4);
            state(5);
            state(6) - g*t;
            0;
            0;
            -g];
end