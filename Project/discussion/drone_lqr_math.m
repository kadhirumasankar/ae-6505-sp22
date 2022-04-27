% symbolic state spacevariables for 12 state drone dynamics
syms x_pos y_pos z_pos x_vel y_vel z_vel roll pitch yaw roll_rate pitch_rate yaw_rate
% symbolic lqr outputs, thrust, roll, pitch, and yaw torque
syms u3 u4 u5 u6
% mass and gravity
syms m g
% I3 = sym('I', [3 3])
% simplifying assumption for symbolic math
% use full nondiagnoal numerical matrix when applying to real drone
syms Ixx Iyy Izz 
I3 = [Ixx, 0,   0;
      0,   Iyy, 0;
      0,   0,   Izz]
% deriv of position variables
X1_dot = [x_vel;
          y_vel;
          z_vel]
% deriv of velocity variables in inertial frame
X2_dot = (R1(roll, pitch, yaw)*[0; 0; u3] - [0;0;m*g])/m
% deriv of roll pitch and yaw in inertial frame
X3_dot = R2(roll, pitch)*[roll_rate; pitch_rate;yaw_rate]
% deriv of angular rates in body frame
X4_dot = inv(I3)*([u4; u5; u6] - cross([roll_rate; pitch_rate;yaw_rate], I3*[roll_rate; pitch_rate;yaw_rate]))
% combined 12 by 1 nonlinear state space equations
X_dot = [X1_dot;
         X1_dot;
         X1_dot;
         X1_dot;
         X2_dot;
         X3_dot;
         X4_dot];
% jacobian takes dericative of 12 by state space eq wrt 12 by state vector
% yielding 12 by 12 nonlinear version of A matrix
A = jacobian(X_dot, [x_pos y_pos z_pos x_pos y_pos z_pos x_pos y_pos z_pos x_pos y_pos z_pos x_vel y_vel z_vel roll pitch yaw roll_rate pitch_rate yaw_rate])
% linear the A matrix about hover thrust (weight), and stable (0) other
% states
for i = 4:size(A,1)
    for j = 1:size(A,2)
        i
        j
        A(i, j)
    end
end
A = taylor(A , [x_pos y_pos z_pos x_vel y_vel z_vel roll pitch yaw roll_rate pitch_rate yaw_rate u3  u4 u5 u6], ...
               [0  0  0  0  0  0  0  0  0  0   0   0   m*g 0  0  0], 'order', 1)
% same steps wrt to controller output to yield B matrix
B = jacobian(X_dot, [u3  u4 u5 u6])
B = taylor(B , [x_pos y_pos z_pos x_vel y_vel z_vel roll pitch yaw roll_rate pitch_rate yaw_rate u3  u4 u5 u6], ...
               [0  0  0  0  0  0  0  0  0  0   0   0   m*g 0  0  0], 'order', 1)

% numerical substitution to yeild usable result
% ofcourse mass and inertial terms need to be replaced with real values
A = double(subs(A, g, 9.81))
B = double(subs(B, [m Ixx Iyy Izz], [1 1 1 1]))
% Q matrix gives weights to importance of regulating each state in the cost
% function
Q = eye(12)
% options for changing subsystem priority
%{
for i = 1:6
    Q(i, i) = .01
end
for i = 10:12
    Q(i, i) = 2
end
%}
% R matrix gives weights to importance of regulating control effort in the
% cost function
% primarily to get control effort within actuator limits
R = eye(4)
% matlab solves algebraic ricatti equation for optimal K matrix
K = lqr(A, B, Q, R)

function R = R1(roll, pitch, yaw)
%x
Rx = [1 0        0;
      0 c(roll) -s(roll);
      0 s(roll)  c(roll)];
%y
Ry = [ c(pitch) 0 s(pitch);
       0        1 0;
      -s(pitch) 0 c(pitch)];
%z
Rz = [c(yaw) -s(yaw) 0;
      s(yaw)  c(yaw) 0;
      0       0      1];
  
 R = simplify(Rz*Ry*Rx)
end

function R = R2(roll, pitch)

R = [1 s(roll)*t(pitch)  c(roll)*t(pitch);
     0 c(roll)          -s(roll);
     0 s(roll)/c(pitch)  c(roll)/c(pitch)]
end


function out = s(th)
out = sin(th);
end

function out = c(th)
out = cos(th);
end

function out = t(th)
out = tan(th);
end




