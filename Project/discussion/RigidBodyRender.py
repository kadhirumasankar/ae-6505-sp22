from matplotlib import pyplot as plt
import numpy as np
import mpl_toolkits.mplot3d.axes3d as p3
from scipy.spatial.transform import Rotation as Rotate


class RigidBodyRender:
    def __init__(self):
        # axis
        fig = plt.figure()
        self.ax = p3.Axes3D(fig)
        # baseline shape for homogeneous transformation
        n_rotor = 4
        n_rotor += 1
        th = np.linspace(0, 2 * np.pi, n_rotor)
        x = np.cos(th).reshape(1, n_rotor)
        y = np.sin(th).reshape(1, n_rotor)
        z = np.zeros(len(th)).reshape(1, n_rotor)
        ones = np.ones(len(th)).reshape(1, n_rotor)
        self.baseline_shape = np.concatenate((x,
                                              y,
                                              z,
                                              ones), axis=0)

    def render(self, x):
        # baseline rigid body shape
        S = self.baseline_shape
        # homogeneous rotation matrix
        H = np.identity(4)
        H[0:3, 0:3] = Rotate.from_euler('zyx', x[8:5:-1]).as_matrix()  # insert rotation
        H[0:3, 3] = x[0:3]  # insert translation
        # rotate and translate the baseline rigid body shape
        S = np.matmul(H, S)
        self.render_transformed(S)

    def render_transformed(self, s):
        self.ax.cla()
        self.ax.plot3D(s[0, :],
                       s[1, :],
                       s[2, :], linewidth=3)

        # axes properties
        lim = 10
        self.ax.set_xlim3d([-lim, lim])
        self.ax.set_ylim3d([-lim, lim])
        self.ax.set_zlim3d([-lim, lim])
        self.ax.set_xlabel('X')
        self.ax.set_ylabel('Y')
        self.ax.set_zlabel('Z')

        # plot and pause
        plt.show(block=False)
        plt.pause(.000001)

