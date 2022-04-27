#!/usr/bin/env python2

import math
import rosbag
import matplotlib.pyplot as plt
from tf.transformations import euler_from_quaternion
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("path")
args = parser.parse_args()

# Change this input to the path to the bag you want to plot
bag = rosbag.Bag(args.path)

# Initialize variables to store /mavros/local_position/pose data
x_pos = []
y_pos = []
z_pos = []
roll = []
pitch = []
yaw = []
# pos_t = []
# Initialize variables to store /mavros/local_position/velocity_body data
x_vel = []
y_vel = []
z_vel = []
p = []
q = []
r = []
# vel_t = []
# Initialize variables to store /lqr_controller/u_star data
thrust = []
roll_torque = []
pitch_torque = []
yaw_torque = []
states_t = []
# Initialize variables to store /mavros_control/motor_speed_setpoint data
w1_setpoint = []
w2_setpoint = []
w3_setpoint = []
w4_setpoint = []
w_setpoint_t = []
w1 = []
w2 = []
w3 = []
w4 = []
w_t = []

x_target = []
y_target = []
z_target = []
xvel_target = []
yvel_target = []
zvel_target = []
t_target = []

# Walks through the bag and stores data in respective variables
for topic, msg, t in bag.read_messages():
    if topic == '/lqr_controller/states':
        if len(msg.data) == 17:
            x_pos.append(msg.data[0])
            y_pos.append(msg.data[1])
            z_pos.append(msg.data[2])
            x_vel.append(msg.data[3])
            y_vel.append(msg.data[4])
            z_vel.append(msg.data[5])
            roll.append(msg.data[6])
            pitch.append(msg.data[7])
            yaw.append(msg.data[8])
            p.append(msg.data[9])
            q.append(msg.data[10])
            r.append(msg.data[11])
            thrust.append(msg.data[12])
            roll_torque.append(msg.data[13])
            pitch_torque.append(msg.data[14])
            yaw_torque.append(msg.data[15])
            states_t.append(msg.data[16]/10**9)
        elif len(msg.data) == 7:
            x_target.append(msg.data[0])
            y_target.append(msg.data[1])
            z_target.append(msg.data[2])
            xvel_target.append(msg.data[3])
            yvel_target.append(msg.data[4])
            zvel_target.append(msg.data[5])
            t_target.append(msg.data[6]/10**9)
        else:
            w1_setpoint.append(msg.data[0])
            w2_setpoint.append(msg.data[1])
            w3_setpoint.append(msg.data[2])
            w4_setpoint.append(msg.data[3])
            w_setpoint_t.append(msg.data[4]/10**9)

for topic, msg, t in bag.read_messages():
    if topic == '/gazebo/rotor_speeds':
        if t.to_sec() >= states_t[0] and t.to_sec() <= states_t[-1]:
            w1.append(msg.rpm.rotor_0)
            w2.append(msg.rpm.rotor_1)
            w3.append(-msg.rpm.rotor_2)
            w4.append(-msg.rpm.rotor_3)
            w_t.append(t.to_sec())
bag.close()

# Plots the states on one window
fig1, axs = plt.subplots(4, 3)
# Removing markers so that I can see the lines
axs[0, 0].plot(states_t, x_pos)
axs[0, 0].set_title("x")
axs[0, 0].grid()
axs[0, 1].plot(states_t, y_pos)
axs[0, 1].set_title("y")
axs[0, 1].grid()
axs[0, 2].plot(states_t, z_pos)
axs[0, 2].set_title("z")
axs[0, 2].grid()
axs[1, 0].plot(states_t, x_vel)
axs[1, 0].set_title("vx")
axs[1, 0].grid()
axs[1, 1].plot(states_t, y_vel)
axs[1, 1].set_title("vy")
axs[1, 1].grid()
axs[1, 2].plot(states_t, z_vel)
axs[1, 2].set_title("vz")
axs[1, 2].grid()
axs[2, 0].plot(states_t, roll)
axs[2, 0].set_title("roll")
axs[2, 0].grid()
axs[2, 1].plot(states_t, pitch)
axs[2, 1].set_title("pitch")
axs[2, 1].grid()
axs[2, 2].plot(states_t, yaw)
axs[2, 2].set_title("yaw")
axs[2, 2].grid()
axs[3, 0].plot(states_t, p)
axs[3, 0].set_title("p")
axs[3, 0].grid()
axs[3, 1].plot(states_t, q)
axs[3, 1].set_title("q")
axs[3, 1].grid()
axs[3, 2].plot(states_t, r)
axs[3, 2].set_title("r")
axs[3, 2].grid()
# Modifies margins and maximizes the plot window
fig1.subplots_adjust(left=0.06, bottom=0.03, right=0.97, top=0.97)
mng = plt.get_current_fig_manager()
mng.resize(*mng.window.maxsize())

# Plots the inputs and motor speeds on another window
fig2, axs2 = plt.subplots(4, 2)
axs2[0, 0].plot(states_t, thrust)
axs2[0, 0].set_title("Thrust")
axs2[1, 0].plot(states_t, roll_torque)
axs2[1, 0].set_title("Roll torque")
axs2[2, 0].plot(states_t, pitch_torque)
axs2[2, 0].set_title("Pitch torque")
axs2[3, 0].plot(states_t, yaw_torque)
axs2[3, 0].set_title("Yaw torque")
axs2[0, 1].plot(w_setpoint_t, w1_setpoint)
axs2[0, 1].set_title("w1_setpoint")
axs2[1, 1].plot(w_setpoint_t, w2_setpoint)
axs2[1, 1].set_title("w2_setpoint")
axs2[2, 1].plot(w_setpoint_t, w3_setpoint)
axs2[2, 1].set_title("w3_setpoint")
axs2[3, 1].plot(w_setpoint_t, w4_setpoint)
axs2[3, 1].set_title("w4_setpoint")
# Modifies margins and maximizes the plot window
fig2.subplots_adjust(left=0.06, bottom=0.03, right=0.97, top=0.97)
mng = plt.get_current_fig_manager()
mng.resize(*mng.window.maxsize())

# Plots the target states
fig3, axs3 = plt.subplots(3, 2)
axs3[0, 0].plot(t_target, x_target)
axs3[0, 0].plot(states_t, x_pos)
axs3[0, 0].set_title("x_target")
axs3[1, 0].plot(t_target, y_target)
axs3[1, 0].plot(states_t, y_pos)
axs3[1, 0].set_title("y_target")
axs3[2, 0].plot(t_target, z_target)
# axs3[2, 0].plot(states_t, z_pos)
axs3[2, 0].set_title("z_target")
axs3[0, 1].plot(t_target, xvel_target)
axs3[0, 1].plot(states_t, x_vel)
axs3[0, 1].set_title("xvel_target")
axs3[1, 1].plot(t_target, yvel_target)
axs3[1, 1].plot(states_t, y_vel)
axs3[1, 1].set_title("yvel_target")
axs3[2, 1].plot(t_target, zvel_target)
axs3[2, 1].plot(states_t, z_vel)
axs3[2, 1].set_title("zvel_target")
# Modifies margins and maximizes the plot window
fig2.subplots_adjust(left=0.06, bottom=0.03, right=0.97, top=0.97)
mng = plt.get_current_fig_manager()
mng.resize(*mng.window.maxsize())

plt.show()
