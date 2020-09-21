# This import registers the 3D projection, but is otherwise unused.
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
import math
import seaborn as sns

sns.set()


# -------------------- Helper Functions are First  -------------------- #

def M(theta):
    """
        Return the angle phi = theta mod (2 pi) such that -pi <= theta < pi.
    """
    theta = theta % (2 * math.pi)
    if theta < -math.pi:
        return theta + 2 * math.pi

    if theta >= math.pi:
        return theta - 2 * math.pi

    return theta


def polar(x, y):
    """
    Return the polar coordinates (r, theta) of the point (x, y).
    """
    r = np.hypot(x, y)  # math.sqrt(x * x + y * y)
    theta = np.arctan2(y, x)  # math.atan2(y, x)
    return r, theta


def change_of_basis(p1, p2):
    """
    Given p1 = (x1, y1, theta1) and p2 = (x2, y2, theta2) represented in a
    coordinate system with origin (0, 0) and rotation 0 (in rad), return
    the position and rotation of p2 in the coordinate system which origin
    (x1, y1) and rotation theta1.
    """
    theta1 = p1[2]
    dx = p2[0] - p1[0]
    dy = p2[1] - p1[1]

    # new_x = dx * math.cos(theta1) + dy * math.sin(theta1)
    # new_y = -dx * math.sin(theta1) + dy * math.cos(theta1)
    xb, yb = rotate_wb(dx, dy, theta1)

    dtheta = p2[2] - p1[2]

    return xb, yb, dtheta


def sign(x):
    return 1 if x >= 0 else -1


def rotate_bw(xb, yb, psi):
    s = np.sin(psi)
    c = np.cos(psi)

    xw = xb * c - yb * s
    yw = xb * s + yb * c

    return xw, yw


def rotate_wb(xw, yw, psi):
    s = np.sin(psi)
    c = np.cos(psi)

    xb = xw * c + yw * s
    yb = -xw * s + yw * c

    return xb, yb


def plot_arrows(q0, q1):
    x0, y0, yaw0 = q0
    x1, y1, yaw1 = q1

    length = 0.6
    width = 0.4
    plt.arrow(x0, y0, length * np.cos(yaw0), length * np.sin(yaw0), head_width=width, head_length=width)
    plt.plot(x0, y0, marker='s', label='start')

    plt.arrow(x1, y1, length * np.cos(yaw1), length * np.sin(yaw1), head_width=width, head_length=width)
    plt.plot(x0, y0, marker='s', label='end')


def plot_arrow(x, y, yaw, length=0.4, width=0.25, fc="r", ec="k", label=''):
    """
        Plot and arrow
    """

    plt.arrow(x, y, length * np.cos(yaw), length * np.sin(yaw),
              fc=fc, ec=ec, head_width=width, head_length=width)

    plt.plot(x, y, marker='s', label=label)


def wrapToPi(angle):
    '''
        wraps to angle -pi< angle < pi
        :param angle: in radians
        :return: wrapped angle
    '''
    temp = np.exp(1j * angle)
    wrapped_angle = np.angle(temp)

    return wrapped_angle


def heat_map3D(table3D):
    '''
        Heat Map
    '''

    # Make data.
    shapes = table3D.shape
    X = np.arange(0, shapes[0])
    Y = np.arange(0, shapes[1])
    X, Y = np.meshgrid(X, Y)

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    for i in range(shapes[2]):
        Z = table3D[:, :, i]

        # Plot the surface.
        surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                               linewidth=0, antialiased=False)

    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=5)
    # Customize the z axis.
    # ax.set_zlim(-1.01, 1.01)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

    plt.show()


def heat_map2D(table2D):
    '''
        Heat Map
    '''

    # Make data.

    # plt.imshow(table2D, cmap='hot', interpolation='nearest')
    plt.figure()
    ax = sns.heatmap(table2D)
    plt.show()


def plot_anim(start, goal, ox, oy):
    # plot obstacles, start and end points
    plt.plot(ox, oy, ".k")
    plt.grid()
    plt.title('Obstacle Coordinates')

    plot_arrow(start[0], start[1], start[2], fc='g')
    plot_arrow(goal[0], goal[1], goal[2])


def plot_xy(x, y):
    # plt.plot(x[0], y[0], 's')
    plt.plot(x, y)
    plt.show()


def check_yaw(x, y, yaw):
    dx = np.diff(x)
    dy = np.diff(y)

    yawxy = np.rad2deg(wrapToPi(np.arctan2(dy, dx)))
    yaw = np.rad2deg(wrapToPi(yaw))

    plt.plot(yaw, label='yaw from rs')
    plt.plot(yawxy, label='xy yaw')
    plt.grid()
    plt.legend()
    plt.show()


def plot_rs(q0, q1, path):
    plt.cla()
    px = path.x
    py = path.y

    plt.plot(px, py, label="final course " + str(path.ctypes))

    # plotting
    start_x, start_y, start_yaw = q0
    end_x, end_y, end_yaw = q1

    plot_arrow(start_x, start_y, start_yaw, label='start')
    plot_arrow(end_x, end_y, end_yaw, label='end')

    plt.legend()
    plt.grid(True)
    plt.axis("equal")

    minx = min(px)
    maxx = max(px)

    miny = min(py)
    maxy = max(py)

    plt.xlim(minx - 3, maxx + 3)
    plt.ylim(miny - 3, maxy + 3)

    # plt.show()


def plot_rs_controls(path):
    # Get the time
    final_time = path.final_time
    t_sim = np.linspace(0, final_time, 100)
    dt = t_sim[1] - t_sim[0]

    # Get the functions 
    acc_func = path.controls['acceleration_func']
    str_fun = path.controls['steering_func']

    # Compute values using the piecewise linear functions
    acc_vals = acc_func(t_sim)
    str_vals = str_fun(t_sim)
    vel_vals = np.cumsum(dt * acc_vals)

    # put the controls in a list to plot
    controls = [acc_vals, str_vals, vel_vals]

    labels = ['acceleration', 'steering', 'velocity']
    fig, axs = plt.subplots(3, 1, sharex='all')

    for i in range(3):
        axs[i].plot(t_sim, controls[i], label=labels[i])
        axs[i].set_title(labels[i])
        plt.legend()
    plt.show()


class PieceWise_Func(object):
    def __init__(self, upper, lower, steering_val):
        self.upper = upper
        self.lower = lower
        self.steering_val = steering_val

    def func(self):
        fa = lambda x: np.piecewise(x, [(x >= self.lower) & (x < self.upper)], [self.steering_val, 0.0])
        return fa


def format_line(name, value, unit=''):
    """
    Formats a line e.g.
    {Name:}           {value}{unit}
    """
    name += ':'
    if isinstance(value, (float, np.ndarray)):
        value = f'{value:{0}.{4}}'

    return f'{name.ljust(40)}{value}{unit}'
