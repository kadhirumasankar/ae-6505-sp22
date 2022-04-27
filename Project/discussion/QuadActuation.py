import numpy as np
from scipy.spatial.transform import Rotation as Rotate


class QuadActuation:
    def __init__(self):
        # actuator properties
        Ct = 1
        Cq = 1
        d = 1
        # matrix related prop speeds to thrust and torque
        self.kinetics = np.zeros((6, 4))
        self.kinetics[0, :] = [0,     0,     0,     0]
        self.kinetics[1, :] = [0,     0,     0,     0]
        self.kinetics[2, :] = [Ct,    Ct,    Ct,    Ct]
        self.kinetics[3, :] = [0,     d*Ct,  0,    -d*Ct]
        self.kinetics[4, :] = [-d*Ct, 0,     d*Ct,  0]
        self.kinetics[5, :] = [-Cq,   Cq,   -Cq,    Cq]
        # lqr feedback matrix from matlab script
        self.K = np.zeros((4, 12))
        self.K[0, :] = [0, 0, 0.1, 0, 0, 0.4583, 0, 0, 0, 0, 0, 0]
        self.K[1, :] = [0, -0.1000, 0, 0, -0.3160, 0, 4.4078, 0, 0, 3.2887, 0, 0]
        self.K[2, :] = [0.1, 0, 0, 0.316, 0, 0, 0, 4.4078, 0, 0, 3.2887, 0]
        self.K[3, :] = [0, 0, 0, 0, 0, 0, 0, 0, 1.0, 0, 0, 2.0]

    # computes inertial force and body torque from prop speeds and orientation
    def inertial_force_body_torque(self, w, x):
        w2 = np.square(w)
        u = np.matmul(self.kinetics, w2).reshape((6, 1))
        Rbi = Rotate.from_euler('zyx', x[8:5:-1]).as_matrix()
        u[0:3, 0] = np.matmul(Rbi, u[0:3, 0])
        return u

    # applies lqr to current state, adds nonlinear gravity offset to thrust, returns desired thrust and torques
    def control(self, x, weight):
        x_des = np.zeros((12,))
        u = np.matmul(-self.K, np.subtract(x, x_des))
        u[0] += weight / (np.cos(x[6])*np.cos(x[7]))
        return u

    # converts desired thrust and torques to prop speeds
    def u_to_w(self, u):
        w2 = np.matmul(np.linalg.inv(self.kinetics[2:6, :]), u)
        w = np.sqrt(w2)
        return w
