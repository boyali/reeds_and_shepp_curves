import numpy as np

# -------------------- Global Parameters ----------------------#
show_animation = False
XY_GRID_RESOLUTION = 1  # [m]
MOTION_RESOLUTION = 0.1  # [m] 0.1
OB_MAP_RESOLUTION = 1  # [m]; obstacle resolution 0.1

YAW_GRID_RESOLUTION = np.deg2rad(15.0)  # [rad]
VR = 1.0  # vehicle radius
Lw = 2.9  # [m]; 7.0 Wheel base - Distance between the axes

# Steering  Parameters
MAX_STEER = np.deg2rad(34)  # [rad]
N_STEER = 3  # number of steer command;

# REEDS AND SHEPP CURVE PARAMETERS
MAX_CURVATURE = np.tan(MAX_STEER) / Lw

# RADCURV = 1
RADCURV = 1 / MAX_CURVATURE  # radius of curvature

RADCURVMUL2 = 2 * RADCURV
RADCURVMUL4 = 4 * RADCURV
SQRADCURV = RADCURV * RADCURV
SQRADCURVMUL2 = 4 * RADCURV * RADCURV

MPI = np.pi
MPIMUL2 = 2 * np.pi
MPIDIV2 = np.pi / 2

EPS = 1e-12
INFINITY = 10000
# HEURISTICS

'''

     USE_NONHOLONOMIC_WITHOUT_OBSTACLE_HEURISTIC
     High cost for the goal if the heading is wrong
     For sparse obstacle density, it works well, for the dence environment, it is sensitive 
'''
USE_NONHOLONOMIC_WITHOUT_OBSTACLE_HEURISTIC = True

'''
    USE_HOLONOMIC_WITH_OBSTACLE_HEURISTIC
        Dynamic Programming in 2D ignoring the non-holonomic nature of the motion
        fx(x, y) = argmin over theta (x, y, theta)
        Detect U-shaped obstacle regions better than the other, 

        max(of two can be used as the heuristic map)

'''

USE_HOLONOMIC_WITH_OBSTACLE_HEURISTIC = True

# Penalty Parameters

H_COST = 5.0  # Heuristic cost for tuning increase h_cost
SB_COST = 100.0  # switch back penalty cost
BACK_COST = 5.0  # backward penalty cost
STEER_CHANGE_COST = 1.0  # steer angle change penalty cost
STEER_COST = 1.0  # steer angle change penalty cost

# Vehicle Parameters
VEHICLE_RADIUS = 1.0  # [m]; radius of rear ball;

'''
    The following global variables are for computing the acceleration profiles

'''
GRAVITY = 9.81  # M/S^2
VMAX = 2  # m/s Maximum Parking Velocity
ACCMAX = 0.05 * GRAVITY  # according to the paper it is 0.12 g for lateral acc
JERKMAX = 0.1 * GRAVITY  # according to the paper it is 0.24 g for lateral acc
