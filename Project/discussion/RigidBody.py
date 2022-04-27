import numpy as np
from numpy import sin as s
from numpy import cos as c
from numpy import tan as t
import matplotlib.pyplot as plt
from gym.spaces import Box
# from scipy.integrate import odeint
from RigidBodyRender import RigidBodyRender
from QuadActuation import QuadActuation


class RigidBody:
    def __init__(self):
        # physical properties
        self.g = 9.81
        self.m = 1
        self.weight = self.m * self.g
        self.I3 = np.identity(3)  # inertia matrix
        self.invI3 = np.linalg.inv(self.I3)  # inverted inertia matrix
        self.properties = (self.m, self.I3, self.invI3)
        # steps
        self.n = 500
        self.tf = 5
        self.time = np.linspace(0, self.tf, self.n)
        self.delta_t = self.tf / self.n
        # state solution place holder w/ initial conditions
        self.x = np.zeros((12, self.n))
        x0 = 2*np.random.rand(12,)-1
        self.x[:, 0] = x0
        self.i = 0
        # render flag
        self.render_is_ready = False
        self.renderer = None
        # Quad
        self.actuation = QuadActuation()
        # Open ai gym params
        self.observation_space = Box(low=-15, high=15, shape=(12,), dtype=np.float32)
        self.action_space = Box(low=0.0, high=4.0, shape=(4,), dtype=np.float32)

    def reset(self):
        self.x = np.zeros((12, self.n))
        x0 = (2*np.random.rand(12,)-1)
        x0[0:3] *= 0
        x0[3:6] *= 4
        x0[6:9] *= .75
        x0[9:12] *= 0
        self.x[:, 0] = x0
        self.i = 0
        return x0

    @staticmethod
    def dynamics(x, t, u, properties):
        # physical properties
        m, I3, invI3 = properties
        # inputs
        F = u[0:3]  # force in inertial frame
        tau = u[3:6]  # torque in body frame
        # sub states
        x = np.reshape(x, (12, 1))
        x1 = x[0:3]  # position in inertial frame
        x2 = x[3:6]  # velocity in inertial frame
        x3 = x[6:9]  # roll pitch yaw in world frame
        x4 = x[9:12]  # angular velocity in body frame
        # sub state derivatives
        x1_dot = x2  # d(pos)/dt = vel
        x2_dot = F / m  # F = ma -> a = F/m
        # x3_dot
        R = rot_w_to_inertial(x3[0], x3[1])
        x3_dot = np.matmul(R, x4)
        # x4_dot
        w_cross_Iw = np.cross(x4, np.matmul(I3, x4), axisa=0, axisb=0).reshape(3, 1)
        x4_dot = np.matmul(invI3, tau - w_cross_Iw)
        # full nonlinear state derivatives
        return np.squeeze(
            np.concatenate((x1_dot,
                            x2_dot,
                            x3_dot,
                            x4_dot), axis=0))
    # step dynamics forward, w is propeller angular speeds
    def step(self, w):

        x0 = self.x[:, self.i] # current 12 states

        # dynamics
        u = self.actuation.inertial_force_body_torque(w, x0)  # compute forces and torques from prop speeds and state
        u[2] -= self. weight  # add gravity
        # ode solver params
        t_span = self.time[self.i:(self.i + 2)]  # time step interval
        """ option for using built in ode solver """
        # params = (u, self.properties)
        # full_span = odeint(self.dynamics, x0, t_span, args=params)
        # self.x[:, i + 1] = full_span[1, :]
        x_dot = self.dynamics(x0, t_span, u, self.properties)  # get state derivatives
        x_ = x0 + x_dot * self.delta_t  # euler integration

        self.x[:, self.i + 1] = x_  # record the new state
        self.i += 1

        # obs, reward, done, info
        # if the solver has not hit a singularity and and no states are out of control
        if (not np.any(np.isnan(x_))) and (not np.any(np.abs(x_) > 50)):
            # normal return
            # -np.sum(np.square(x_)) for regulating all states
            return x_, -np.sum(np.square(x_[0:3])), self.time[self.i] >= self.tf, None
        else:
            # normal return with huge negative reward
            return self.x[:, self.i-1], -np.sum(np.square(self.x[0:3, self.i-1]))-100000, True, None

    def render(self):
        if not self.render_is_ready:
            self.render_is_ready = True
            self.renderer = RigidBodyRender()
        self.renderer.render(self.x[:, self.i])

    def plot(self):
        fig, axs = plt.subplots(2, 2)
        axs[0, 0].plot(self.time, self.x[0:3, :].T)
        axs[1, 0].plot(self.time, self.x[3:6, :].T)
        axs[0, 1].plot(self.time, self.x[6:9, :].T)
        axs[1, 1].plot(self.time, self.x[9:12, :].T)
        plt.show(block=False)


def rot_w_to_inertial(roll, pitch):
    # rotate body angular rates to the inertial frame
    R = np.identity(3)
    R[0, :] = [1, s(roll)*t(pitch),  c(roll)*t(pitch)]
    R[1, :] = [0, c(roll)         , -s(roll)         ]
    R[2, :] = [0, s(roll)/c(pitch),  c(roll)/c(pitch)]
    return R


if __name__ == '__main__':

    drone = RigidBody()
    for j in range(100):
        drone.reset()

        for i in range(drone.n - 1):
            x0 = drone.x[:, drone.i]

            # control input
            u = drone.actuation.control(x0, drone.weight)  # desired force and torque
            w = drone.actuation.u_to_w(u)  # resultant prop speeds

            drone.step(w)
            if (i % 20) == 0: # mod to increase render speed
                drone.render()

        # optional plot
        # drone.plot()
        # plt.pause(5)
        # plt.close()
        # plt.close('all')




