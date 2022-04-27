"""
Code to read our ROSbags and produce CSV files.
This will make it easier to produce and analyze plots. Note that this script will
probably not work with ROSbags that are laid out differently.
"""

#!/usr/bin/env python2

# import math
import csv
import rosbag
import matplotlib.pyplot as plt

# from tf.transformations import euler_from_quaternion
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

# Walks through the bag and stores data in respective variables
for topic, msg, t in bag.read_messages():
    if topic == "/lqr_controller/states":
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
            states_t.append(msg.data[16] / 10 ** 9)
        else:
            w1_setpoint.append(msg.data[0])
            w2_setpoint.append(msg.data[1])
            w3_setpoint.append(msg.data[2])
            w4_setpoint.append(msg.data[3])
            w_setpoint_t.append(msg.data[4] / 10 ** 9)

for topic, msg, t in bag.read_messages():
    if topic == "/gazebo/rotor_speeds":
        if t.to_sec() >= states_t[0] and t.to_sec() <= states_t[-1]:
            w1.append(msg.rpm.rotor_0)
            w2.append(msg.rpm.rotor_1)
            w3.append(-msg.rpm.rotor_2)
            w4.append(-msg.rpm.rotor_3)
            w_t.append(t.to_sec())
bag.close()

with open(args.path[13:-4] + ".csv", "w") as csvfile:
    writer = csv.writer(csvfile, delimiter=",", quotechar="|", quoting=csv.QUOTE_NONE)
    writer.writerow(
        [
            "states_t",
            "x_pos",
            "y_pos",
            "z_pos",
            "x_vel",
            "y_vel",
            "z_vel",
            "roll",
            "pitch",
            "yaw",
            "p",
            "q",
            "r",
            "thrust",
            "roll_torque",
            "pitch_torque",
            "yaw_torque",
            "w_setpoint_t",
            "w1_setpoint",
            "w2_setpoint",
            "w3_setpoint",
            "w4_setpoint",
            "w_t",
            "w1",
            "w2",
            "w3",
            "w4",
        ]
    )
    writer.writerow(states_t)
    writer.writerow(x_pos)
    writer.writerow(y_pos)
    writer.writerow(z_pos)
    writer.writerow(x_vel)
    writer.writerow(y_vel)
    writer.writerow(z_vel)
    writer.writerow(roll)
    writer.writerow(pitch)
    writer.writerow(yaw)
    writer.writerow(p)
    writer.writerow(q)
    writer.writerow(r)
    writer.writerow(thrust)
    writer.writerow(roll_torque)
    writer.writerow(pitch_torque)
    writer.writerow(yaw_torque)
    writer.writerow(w_setpoint_t)
    writer.writerow(w1_setpoint)
    writer.writerow(w2_setpoint)
    writer.writerow(w3_setpoint)
    writer.writerow(w4_setpoint)
    writer.writerow(w_t)
    writer.writerow(w1)
    writer.writerow(w2)
    writer.writerow(w3)
    writer.writerow(w4)
