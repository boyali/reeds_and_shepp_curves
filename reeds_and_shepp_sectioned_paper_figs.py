import math

import matplotlib.pyplot as plt
import numpy as np
from rs_global_params import *
from utils import *
from path_dictionary import path_dict

import pickle

## Save History to this file
logs_pickle = dict()

"""
Reeds Shepp path planner sample code
author Ali Boyali

from the paper:
Optimal Paths for a Car That Goes Both Forwards and Backwards

and the code:
Check the followings
http://msl.cs.uiuc.edu/~lavalle/cs326a/rs.c
Steven Lavalle Planning Algorithms


Paper Summary:

u is the angular turning length : pi/2 or similar
t is the turning direction

Manifold parameters: t, u, v , w  --> sequence indices (lr, ru, sv, lw) or Ct, Cu, Cv,Sw

L(t, u, v, w) --> Length function

"""

show_animation = True


class Path:
    steering_map = {'L': MAX_STEER, 'R': -MAX_STEER, 'S': 0}

    def __init__(self):
        self.curve_num = []  # curve number
        self.lengths = []  # maneuver lengths of RS motion primitives, in radian and meters n
        self.lengths_expanded = []  # homogenous lengths all in meters assigned in append control method
        self.ctypes = []  # maneuver code names Left Right Straight
        self.total_length = 0.0  # total length of the maneuver
        self.x = []  # x coordinates list
        self.y = []
        self.yaw = []
        self.directions = []  # whether forward or backward \in [-1, 1]

        # Controls, final time and final distance
        '''
            Controls are paremetrized by the curve lengths 
        '''
        self.controls = {}  # with the keys acceleration, steering and acceleration-steering
        self.final_time = None
        self.final_distance = None

        '''
            We store steering and accelaration as piecewise functions, acc function is a function of time, 
            but the steering function as a function of distance. Later to transform the steering conditions from 
            distance to time, we need to store the steering conditions list [s0, sfinal, steering_val] --> [t0, tf, 
            steering_val]
        '''
        self.steering_conditions = None

        # we also store each curves coordinates separately
        self.man_patches = []  # [[px, py, pyaw], [px, py, pyaw], ....]

    def set_coords(self, x, y, yaw):
        self.x = x
        self.y = y
        self.yaw = yaw

    def set_directions(self, dir):
        self.directions = dir


def generate_path(q0, q1):
    x, y, phi = change_of_basis(q0, q1)

    paths = []
    paths = CSC(x, y, phi, paths)
    paths = CCC(x, y, phi, paths)
    paths = CCCC(x, y, phi, paths)
    paths = CCSC(x, y, phi, paths)
    paths = CSCC(x, y, phi, paths)
    paths = CCSCC(x, y, phi, paths)

    return paths


def set_path(paths, lengths, ctypes, length_rs, num):
    path = Path()
    path.curve_num = num
    path.ctypes = ctypes

    path.lengths = lengths
    path.total_length = length_rs

    # check same path exist
    for tpath in paths:
        typeissame = (tpath.ctypes == path.ctypes)

        if typeissame:
            if sum(tpath.lengths) - sum(path.lengths) <= 0.01:
                return paths  # not insert path

    # Base.Test.@test path.L >= 0.01
    if path.total_length >= 0.01:
        paths.append(path)

    return paths


def curve_params(phi):
    sphi = np.sin(phi)
    cphi = np.cos(phi)

    ap = RADCURV * sphi
    am = -RADCURV * sphi
    b1 = RADCURV * (cphi - 1)
    b2 = RADCURV * (cphi + 1)

    return [sphi, cphi, ap, am, b1, b2]


def time_fli(args):
    x, y, phi = args

    return -x, y, -phi


def reflect(args):
    x, y, phi = args

    return x, -y, -phi


def reflect_flip(args):
    x, y, phi = args

    return -x, -y, phi


# -------------------------------------- CURVES STARTS HERE ----------------------------------------------------
# ------------------------------------------ CCC CURVES --------------------------------------------------------
# -------------------------------------- SAME TURN CURVES ------------------------------------------------------


def c_c_c(x, y, phi, rs, rc):
    a = x - rs
    b = y + rc

    if np.abs(a) < EPS and np.abs(b) < EPS:
        return False, 0.0, 0.0, 0.0, INFINITY

    u1, theta = polar(a, b)

    if u1 < RADCURVMUL4:
        alpha = np.arccos(u1 / RADCURVMUL4)
        t = wrapToPi(MPIDIV2 + alpha + theta)
        u = wrapToPi(MPI - 2 * alpha)
        v = wrapToPi(phi - t - u)
        length_rs = RADCURV * (t + v + u)

        if t >= 0 and u >= 0 and v >= 0:
            return True, t, u, v, length_rs

    return False, 0.0, 0.0, 0.0, INFINITY


def c_cc(x, y, phi, rs, rc):
    a = x - rs
    b = y + rc

    if np.abs(a) < EPS and np.abs(b) < EPS:
        return False, 0.0, 0.0, 0.0, INFINITY

    u1, theta = polar(a, b)

    if (u1 <= RADCURVMUL4):
        alpha = np.arccos(u1 / RADCURVMUL4)
        t = wrapToPi(MPIDIV2 + alpha + theta)
        u = wrapToPi(MPI - 2 * alpha)
        v = wrapToPi(t + u - phi)
        length_rs = RADCURV * (t + v + u)

        if t >= 0 and u >= 0 and v >= 0:
            return True, t, u, v, length_rs

    return False, 0.0, 0.0, 0.0, INFINITY


def cc_c(x, y, phi, rs, rc):
    a = x - rs
    b = y + rc

    if np.abs(a) < EPS and np.abs(b) < EPS:
        return False, 0.0, 0.0, 0.0, INFINITY

    u1, theta = polar(a, b)

    if (u1 <= RADCURVMUL4):
        u = np.arccos((8 * SQRADCURV - u1 * u1) / (8 * SQRADCURV))
        va = np.sin(u)

        if np.abs(va) < 0.001:
            va = 0.0

        if np.abs(va) < 0.001 and np.abs(u1) < 0.001:
            return False, 0.0, 0.0, 0.0, INFINITY

        alpha = np.arcsin(RADCURVMUL2 * va / u1)
        t = wrapToPi(MPIDIV2 - alpha + theta)
        v = wrapToPi(t - u - phi)
        length_rs = RADCURV * (t + u + v)

        if t >= 0 and u >= 0 and v >= 0:
            return True, t, u, v, length_rs

    return False, 0.0, 0.0, 0.0, INFINITY


def CCC(x, y, phi, paths):
    """
        Formula 8.3: C|C|C
    """

    sphi, cphi, ap, am, b1, b2 = curve_params(phi)
    # num 1     : ["L+", "R-", "L+"], arguments -- c_c_c(x,y,phi,ap,b1)
    flag, t, u, v, length_rs = c_c_c(x, y, phi, ap, b1)
    if flag:
        num = 1
        cytpe = path_dict[num]
        paths = set_path(paths, [t, u, v], cytpe, length_rs, num)

    #  num 2     : ["L-", "R+", "L-"], arguments -- c_c_c(-x,y,-phi,am,b1)   -> time flip of num1
    flag, t, u, v, length_rs = c_c_c(-x, y, -phi, am, b1)
    if flag:
        num = 2
        cytpe = path_dict[num]
        paths = set_path(paths, [-t, -u, -v], cytpe, length_rs, num)

    # num 3     : ["R+", "L-", "R+"], arguments -- c_c_c(x,-y,-phi,am,b1)   -> reflection of num1
    flag, t, u, v, length_rs = c_c_c(x, -y, -phi, am, b1)
    if flag:
        num = 3
        cytpe = path_dict[num]
        paths = set_path(paths, [t, u, v], cytpe, length_rs, num)

    # num 4     : ["R-", "L+", "R-"], arguments -- c_c_c(-x,-y,phi,ap,b1)   -> time flip + reflection of num 1
    flag, t, u, v, length_rs = c_c_c(-x, -y, phi, ap, b1)
    if flag:
        num = 4
        cytpe = path_dict[num]
        paths = set_path(paths, [-t, -u, -v], cytpe, length_rs, num)

    # ------------ C | CC ----_
    """
        Formula 8.4 (1): C|CC
    """

    # num 5     : ["L+", "R-", "L-"], arguments -- c_cc(x,y,phi,ap,b1)        -> is given
    flag, t, u, v, length_rs = c_cc(x, y, phi, ap, b1)
    if flag:
        num = 5
        cytpe = path_dict[num]
        paths = set_path(paths, [t, u, v], cytpe, length_rs, num)

    # num 6     : ["L-", "R+", "L+"], arguments -- c_cc(-x,y,-phi,am,b1)      -> time flip of num 5
    flag, t, u, v, length_rs = c_cc(-x, y, -phi, am, b1)
    if flag:
        num = 6
        cytpe = path_dict[num]
        paths = set_path(paths, [-t, -u, -v], cytpe, length_rs, num)

    # # num 7     : ["R+", "L-", "R-"], arguments -- c_cc(x,-y,-phi,am,b1)   -> reflection of num 5
    flag, t, u, v, length_rs = c_cc(x, -y, -phi, am, b1)
    if flag:
        num = 7
        cytpe = path_dict[num]
        paths = set_path(paths, [t, u, v], cytpe, length_rs, num)

    # # num 8     : ["R-", "L+", "R+"], arguments -- c_cc(-x,-y,phi,ap,b1)
    flag, t, u, v, length_rs = c_cc(-x, -y, phi, ap, b1)
    if flag:
        num = 8
        cytpe = path_dict[num]
        paths = set_path(paths, [t, u, v], cytpe, length_rs, num)

    # ------------ CC|C ----_
    """
        Formula 8.4 (2): CC|C
    """
    #  num 37    : ["L+", "R+", "L-"], arguments -- cc_c(x,y,phi,ap,b1)        -> given
    flag, t, u, v, length_rs = cc_c(x, y, phi, ap, b1)
    if flag:
        num = 37
        cytpe = path_dict[num]
        paths = set_path(paths, [t, u, v], cytpe, length_rs, num)

    #  num 38    : ["R+", "L+", "R-"], arguments -- cc_c(x,-y,-phi,am,b1)  -> reflection of num 37
    flag, t, u, v, length_rs = cc_c(x, -y, -phi, am, b1)
    if flag:
        num = 38
        cytpe = path_dict[num]
        paths = set_path(paths, [t, u, v], cytpe, length_rs, num)

    # num 39    : ["L-", "R-", "L+"], arguments -- cc_c(-x,y,-phi,am,b1)      -> time flip of num 37
    flag, t, u, v, length_rs = cc_c(-x, y, -phi, am, b1)
    if flag:
        num = 39
        cytpe = path_dict[num]
        paths = set_path(paths, [-t, -u, -v], cytpe, length_rs, num)

    # num 40    : ["R-", "L-", "R+"], arguments -- cc_c(-x,-y,phi,ap,b1)      -> time flip + reflection
    flag, t, u, v, length_rs = cc_c(-x, -y, phi, ap, b1)
    if flag:
        num = 40
        cytpe = path_dict[num]
        paths = set_path(paths, [-t, -u, -v], cytpe, length_rs, num)

    return paths


# ------------------------------------------ CSC CURVES  ------------------------------------------------------
def csca(x, y, phi, rs, rc):
    """
        Formula 8.1: CSC (same turns)
    """
    a = x - rs
    b = y + rc

    u, t = polar(a, b)  # polar coordinates returns r, theta
    v = wrapToPi(phi - t)

    if t >= 0 and u >= 0 and v >= 0:
        length_rs = RADCURV * (t + v) + u
        return True, t, u, v, length_rs

    return False, 0.0, 0.0, 0.0, INFINITY


def cscb(x, y, phi, rs, rc):
    """
        Formula 8.2: CSC (opposite turns)
    """
    a = x + rs
    b = y - rc

    rho, t1 = polar(a, b)
    if rho >= RADCURVMUL2:
        u = np.sqrt(rho * rho - SQRADCURVMUL2)
        alpha = np.arctan2(RADCURVMUL2, u)
        t = wrapToPi(t1 + alpha)
        v = wrapToPi(t - phi)

        if t >= 0 and u >= 0 and v >= 0:
            length_rs = RADCURV * (t + v) + u
            return True, t, u, v, length_rs

    return False, 0.0, 0.0, 0.0, INFINITY


def CSC(x, y, phi, paths):
    sphi, cphi, ap, am, b1, b2 = curve_params(phi)

    # ----------------------- SAME TURN CURVES ------------------------------------------------------
    # num 9     : ["L+", "S+", "L+"], arguments -- csca(x,y,phi,ap,b1)
    flag, t, u, v, length_rs = csca(x, y, phi, ap, b1)

    if flag:
        num = 9
        cytpe = path_dict[num]
        paths = set_path(paths, [t, u, v], cytpe, length_rs, num)

    # num 10    : ["R+", "S+", "R+"], arguments -- csca(x,-y,-phi,am,b1)
    flag, t, u, v, length_rs = csca(x, -y, -phi, am, b1)
    if flag:
        num = 10
        cytpe = path_dict[num]
        paths = set_path(paths, [t, u, v], cytpe, length_rs, num)  # reflection of num 9

    # num 11    : ["L-", "S-", "L-"], arguments -- csca(-x,y,-phi, am, b1)
    flag, t, u, v, length_rs = csca(-x, y, -phi, am, b1)
    if flag:
        num = 11
        cytpe = path_dict[num]
        paths = set_path(paths, [-t, -u, -v], cytpe, length_rs, num)  # time flip of num 9

    # num 12    : ["R-", "S-", "R-"], arguments -- csca(-x,-y,phi,ap,b1) -> time flip and reflection of num 9
    flag, t, u, v, length_rs = csca(-x, -y, phi, ap, b1)
    if flag:
        num = 12
        cytpe = path_dict[num]
        paths = set_path(paths, [-t, -u, -v], cytpe, length_rs, num)  # both time flip + reflection

    # ----------------------- OPPOSITE TURN CURVES ------------------------------------------------------
    # num 13    : ["L+", "S+", "R+"], arguments -- cscb(x,y,phi,ap,b2)
    flag, t, u, v, length_rs = cscb(x, y, phi, ap, b2)
    if flag:
        num = 13
        cytpe = path_dict[num]
        paths = set_path(paths, [t, u, v], cytpe, length_rs, num=13)

    # num 14    : ["R+", "S+", "L+"], arguments -- cscb(x,-y,-phi,am,b2) reflection of num 13
    flag, t, u, v, length_rs = cscb(x, -y, -phi, am, b2)
    if flag:
        num = 14
        cytpe = path_dict[num]
        paths = set_path(paths, [t, u, v], cytpe, length_rs, num)  # reflection of num 13

    # num 15    : ["L-", "S-", "R-"], arguments -- cscb(-x,y,-phi,am,b2) -> time flip of num 13
    flag, t, u, v, length_rs = cscb(-x, y, -phi, am, b2)
    if flag:
        num = 15
        cytpe = path_dict[num]
        paths = set_path(paths, [-t, -u, -v], cytpe, length_rs, num)  # time flip of curve 13

    # num 16    : ["R-", "S-", "L-"], arguments -- cscb(-x,-y,phi,ap,b2)
    flag, t, u, v, length_rs = cscb(-x, -y, phi, ap, b2)  # time flip and reflection of 13
    if flag:
        num = 16
        cytpe = path_dict[num]
        paths = set_path(paths, [-t, -u, -v], cytpe, length_rs, num)  # time flip of curve 14

    return paths


# ------------------------------------------ CCCC CURVES  ------------------------------------------------------

def ccu_cuc(x, y, phi, rs, rc):
    """
        Formula 8.7: CCu|CuC
    """
    a = x + rs
    b = y - rc

    if np.abs(a) < EPS and np.abs(b) < EPS:
        return False, 0.0, 0.0, 0.0, INFINITY

    u1, theta = polar(a, b)

    if u1 <= RADCURVMUL4:
        if u1 > RADCURVMUL2:
            alpha = np.arccos((u1 / 2 - RADCURV) / RADCURVMUL2)
            t = wrapToPi(MPIDIV2 + theta - alpha)
            u = wrapToPi(MPI - alpha)
            v = wrapToPi(phi - t + 2 * u)

        else:
            alpha = np.arccos((u1 / 2 + RADCURV) / RADCURVMUL2)
            t = wrapToPi(MPIDIV2 + theta + alpha)
            u = wrapToPi(alpha)
            v = wrapToPi(phi - t + 2 * u)

        length_rs = RADCURV * (2 * u + t + v)

        if t >= 0 and u >= 0 and v >= 0:
            return True, t, u, v, length_rs

    return False, 0.0, 0.0, 0.0, INFINITY


def c_cucu_c(x, y, phi, rs, rc):
    """
        Formula 8.7: CCu|CuC
    """
    a = x + rs
    b = y - rc

    if np.abs(a) < EPS and np.abs(b) < EPS:
        return False, 0.0, 0.0, 0.0, INFINITY

    u1, theta = polar(a, b)
    if u1 > 6 * RADCURV:
        return False, 0.0, 0.0, 0.0, INFINITY

    va1 = (5 * SQRADCURV - u1 * u1 / 4) / SQRADCURVMUL2

    if va1 < 0 or va1 > 1:
        return False, 0.0, 0.0, 0.0, INFINITY

    u = np.arccos(va1)
    va2 = np.sin(u)

    alpha = np.arcsin(RADCURVMUL2 * va2 / u1)
    t = wrapToPi(MPIDIV2 + theta + alpha)
    v = wrapToPi(t - phi)
    length_rs = RADCURV * (2 * u + t + v)

    if t >= 0 and u >= 0 and v >= 0:
        return True, t, u, v, length_rs

    return False, 0.0, 0.0, 0.0, INFINITY


def CCCC(x, y, phi, paths):
    sphi, cphi, ap, am, b1, b2 = curve_params(phi)

    ########## -------------------------------- C Cu | Cu C  Curves# -------------------------------------------------
    """
    Formula 8.7: CCu|CuC
    """

    # num 17     : ["L+", "R+", "L-", "R-"], arguments -- cu_cuc(x,y,phi,ap,b2)         -> is given
    flag, t, u, v, length_rs = ccu_cuc(x, y, phi, ap, b2)
    if flag:
        num = 17
        cytpe = path_dict[num]
        paths = set_path(paths, [t, u, u, v], cytpe, length_rs, num)

    # num 18     : ["R+", "L+", "R-", "L-"], arguments -- ccu_cuc(x,-y,-phi,am,b2)
    flag, t, u, v, length_rs = ccu_cuc(x, -y, -phi, am, b2)
    if flag:
        num = 18
        cytpe = path_dict[num]
        paths = set_path(paths, [t, u, u, v], cytpe, length_rs, num)  # -> reflection of num 17

    #  num 19     : ["L-", "R-", "L+", "R+"], arguments -- ccu_cuc(-x,y,-phi,am,b2)
    flag, t, u, v, length_rs = ccu_cuc(-x, y, -phi, am, b2)
    if flag:
        num = 19
        cytpe = path_dict[num]
        paths = set_path(paths, [-t, -u, -u, -v], cytpe, length_rs, num)  # -> time filip of num 17

    # num 20: ["R-", "L-", "R+", "L+"], arguments -- ccu_cuc(-x,-y,phi,ap,b2) -> time flip and reflect of num 17
    flag, t, u, v, length_rs = ccu_cuc(-x, -y, phi, ap, b2)
    if flag:
        num = 20
        cytpe = path_dict[num]
        paths = set_path(paths, [-t, -u, -u, -v], cytpe, length_rs, num)  # -> time filip of num 17

    ########## ---------------------------------- C | Cu Cu | C     Curves --------------------------------------------
    """
         Formula 8.8: C|CuCu|C
    """

    # num 21     : ["L+", "R-", "L-", "R+"], arguments -- c_cucu_c(x,y,phi,ap,b2)
    flag, t, u, v, length_rs = c_cucu_c(x, y, phi, ap, b2)
    if flag:
        num = 21
        cytpe = path_dict[num]
        paths = set_path(paths, [t, u, u, v], cytpe, length_rs, num)  # -> given

    #  num 22     : ["R+", "L-", "R-", "L+"], arguments -- c_cucu_c(x,-y,-phi,am,b2)
    flag, t, u, v, length_rs = c_cucu_c(x, -y, -phi, am, b2)
    if flag:
        num = 22
        cytpe = path_dict[num]
        paths = set_path(paths, [t, u, u, v], cytpe, length_rs, num)  # ---> reflection of num 21

    # num 23     : ["L-", "R+", "L+", "R-"], arguments -- c_cucu_c(-x,y,-phi,am,b2) --> time flip of num 21
    flag, t, u, v, length_rs = c_cucu_c(-x, y, -phi, am, b2)
    if flag:
        num = 23
        cytpe = path_dict[num]
        paths = set_path(paths, [-t, -u, -u, -v], cytpe, length_rs, num)  # -> time flip of num 21

    # num 24     : ["R-", "L+", "R+", "L-"], arguments -- c_cucu_c(-x,-y,phi,ap,b2)
    flag, t, u, v, length_rs = c_cucu_c(-x, -y, phi, ap, b2)
    if flag:
        num = 24
        cytpe = path_dict[num]
        paths = set_path(paths, [-t, -u, -u, -v], cytpe, length_rs, num)  # -> time flip and reflection

    return paths


# ------------------------------------------ CCSC CURVES  ------------------------------------------------------
def c_c2sca(x, y, phi, rs, rc):
    """
     Formula 8.9 (1): C|C[pi/2]SC
     """
    a = x - rs
    b = y + rc
    u1, theta = polar(a, b)

    if u1 >= RADCURVMUL2:
        u = np.sqrt(u1 ** 2 - SQRADCURVMUL2) - RADCURVMUL2
        if u >= 0:
            alpha = np.arctan2(RADCURVMUL2, (u + RADCURVMUL2))
            t = wrapToPi(MPIDIV2 + theta + alpha)
            v = wrapToPi(t + MPIDIV2 - phi)
            length_rs = RADCURV * (t + MPIDIV2 + v) + u

            if t >= 0 and u >= 0 and v >= 0:
                return True, t, u, v, length_rs

    return False, 0.0, 0.0, 0.0, INFINITY


def c_c2scb(x, y, phi, rs, rc):
    """
     Formula 8.9 (1): C|C[pi/2]SC
     """
    a = x + rs
    b = y - rc
    u1, theta = polar(a, b)

    if u1 >= RADCURVMUL2:

        t = wrapToPi(MPIDIV2 + theta)
        u = u1 - RADCURVMUL2
        v = wrapToPi(phi - t - MPIDIV2)

        length_rs = RADCURV * (t + MPIDIV2 + v) + u

        if t >= 0 and u >= 0 and v >= 0:
            return True, t, u, v, length_rs

    return False, 0.0, 0.0, 0.0, INFINITY


def CCSC(x, y, phi, paths):
    sphi, cphi, ap, am, b1, b2 = curve_params(phi)

    ######## ------------------------------- C | C[pi/2] S C    Curves -----------------------------------------
    """
     Formula 8.9 (1): C|C[pi/2]SC
     """
    # ----------------------------------- SAME TURN MANEUVERS -----------------------------------------------------

    #  num 25     : ["L+", "R-", "S-", "L-"], arguments --  c_c2sca(x,y,phi,ap,b1)        -> given
    flag, t, u, v, length_rs = c_c2sca(x, y, phi, ap, b1)
    if flag:
        num = 25
        cytpe = path_dict[num]
        paths = set_path(paths, [t, MPIDIV2, u, v], cytpe, length_rs, num)

    # num 26     : ["R+", "L-", "S-", "R-"], arguments --  c_c2sca(x,-y,-phi,am,b1)
    flag, t, u, v, length_rs = c_c2sca(x, -y, -phi, am, b1)
    if flag:
        num = 26
        cytpe = path_dict[num]
        paths = set_path(paths, [t, MPIDIV2, u, v], cytpe, length_rs, num)  # -> reflection of num 25

    # num 27     : ["L-", "R+", "S+", "L+"], arguments --  c_c2sca(-x,y,-phi,am,b1)      -> time flip of num 25
    flag, t, u, v, length_rs = c_c2sca(-x, y, -phi, am, b1)
    if flag:
        num = 27
        cytpe = path_dict[num]
        paths = set_path(paths, [-t, -MPIDIV2, -u, -v], cytpe, length_rs, num)

    # num 28     : ["R-", "L+", "S+", "R+"], arguments --  c_c2sca(-x,-y,phi,ap,b1) -> time flip of num 26
    flag, t, u, v, length_rs = c_c2sca(-x, -y, phi, ap, b1)
    if flag:
        num = 28
        cytpe = path_dict[num]
        paths = set_path(paths, [-t, -MPIDIV2, -u, -v], cytpe, length_rs, num)

    # ----------------------------------------OPPOSITE TURN MANEUVERS -------------------------------------------

    # num 29     : ["L+", "R-", "S-", "R-"], arguments --  c_c2scb(x,y,phi,ap,b2)        -> given
    flag, t, u, v, length_rs = c_c2scb(x, y, phi, ap, b2)
    if flag:
        num = 29
        cytpe = path_dict[num]
        paths = set_path(paths, [t, MPIDIV2, u, v], cytpe, length_rs, num)

    # num 30     : ["R+", "L-", "S-", "L-"], arguments --  c_c2scb(x,-y,-phi,am,b2)      -> reflection of num 29
    flag, t, u, v, length_rs = c_c2scb(x, -y, -phi, am, b2)
    if flag:
        num = 30
        cytpe = path_dict[num]
        paths = set_path(paths, [t, MPIDIV2, u, v], cytpe, length_rs, num)

    # num 31     : ["L-", "R+", "S+", "R+"], arguments --  c_c2scb(-x,y,-phi,am,b2)      -> time flip
    flag, t, u, v, length_rs = c_c2scb(-x, y, -phi, am, b2)
    if flag:
        num = 31
        cytpe = path_dict[num]
        paths = set_path(paths, [-t, -MPIDIV2, -u, -v], cytpe, length_rs, num)

    # num 32     : ["R-", "L+", "S+", "L+"], arguments --  c_c2scb(-x,-y,phi,ap,b2)      -> time flip + reflection
    flag, t, u, v, length_rs = c_c2scb(-x, -y, phi, ap, b2)
    if flag:
        num = 32
        cytpe = path_dict[num]
        paths = set_path(paths, [-t, -MPIDIV2, -u, -v], cytpe, length_rs, num)

    return paths


# ------------------------------------------ CSCC CURVES  ------------------------------------------------------
def csc2_ca(x, y, phi, rs, rc):
    """
        Formula 8.9 (2): CSC[pi/2]|C
    """
    a = x - rs
    b = y + rc
    u1, theta = polar(a, b)

    if u1 >= RADCURVMUL2:
        u = np.sqrt(u1 ** 2 - SQRADCURVMUL2) - RADCURVMUL2
        if u >= 0:
            alpha = np.arctan2(u + RADCURVMUL2, RADCURVMUL2)
            t = wrapToPi(MPIDIV2 + theta - alpha)
            v = wrapToPi(t - MPIDIV2 - phi)
            length_rs = RADCURV * (t + MPIDIV2 + v) + u

            if t >= 0 and u >= 0 and v >= 0:
                return True, t, u, v, length_rs

    return False, 0.0, 0.0, 0.0, INFINITY


def csc2_cb(x, y, phi, rs, rc):
    """
        Formula 8.10 (2): CSC[pi/2]|C
    """
    a = x + rs
    b = y - rc
    u1, theta = polar(a, b)

    if u1 >= RADCURVMUL2:
        t = wrapToPi(theta)
        u = u1 - RADCURVMUL2
        v = wrapToPi(-t - MPIDIV2 + phi)

        length_rs = RADCURV * (t + MPIDIV2 + v) + u

        if t >= 0 and u >= 0 and v >= 0:
            return True, t, u, v, length_rs

    return False, 0.0, 0.0, 0.0, INFINITY


def CSCC(x, y, phi, paths):
    """
        Formula 8.9 (2): CSC[pi/2]|C
    """
    sphi, cphi, ap, am, b1, b2 = curve_params(phi)

    # ----------------------------------- SAME TURN MANEUVERS -----------------------------------------------------

    #  num 41     : ["L+", "S+", "R+", "L-"], arguments -- csc2_ca(x,y,phi,ap,b1)         -> given
    flag, t, u, v, length_rs = csc2_ca(x, y, phi, ap, b1)
    if flag:
        num = 41
        cytpe = path_dict[num]
        paths = set_path(paths, [t, u, MPIDIV2, v], cytpe, length_rs, num)

    #  num 42     : ["R+", "S+", "L+", "R-"], arguments -- csc2_ca(x,-y,-phi,am,b1)       -> reflection
    flag, t, u, v, length_rs = csc2_ca(x, -y, -phi, am, b1)
    if flag:
        num = 42
        cytpe = path_dict[num]
        paths = set_path(paths, [t, u, MPIDIV2, v], cytpe, length_rs, num)

    #  num 43     : ["L-", "S-", "R-", "L+"], arguments -- csc2_ca(-x,y,-phi,am,b1)       -> time flip
    flag, t, u, v, length_rs = csc2_ca(-x, y, -phi, am, b1)
    if flag:
        num = 43
        cytpe = path_dict[num]
        paths = set_path(paths, [-t, -u, -MPIDIV2, -v], cytpe, length_rs, num)

    #  num 44     : ["R-", "S-", "L-", "R+"], arguments -- csc2_ca(-x,-y,phi,ap,b1)       -> time flip + reflection
    flag, t, u, v, length_rs = csc2_ca(-x, -y, phi, ap, b1)
    if flag:
        num = 44
        cytpe = path_dict[num]
        paths = set_path(paths, [-t, -u, -MPIDIV2, -v], cytpe, length_rs, num)

    # ----------------------------------- SAME TURN MANEUVERS -----------------------------------------------------
    #  num 45     : ["L+", "S+", "L+", "R-"], arguments -- csc2_cb(x,y,phi,ap,b2)         -> given
    flag, t, u, v, length_rs = csc2_cb(x, y, phi, ap, b2)
    if flag:
        num = 45
        cytpe = path_dict[num]
        paths = set_path(paths, [t, u, MPIDIV2, v], cytpe, length_rs, num)

    # num 46     : ["R+", "S+", "R+", "L-"], arguments -- csc2_cb(x,-y,-phi,am,b2)       -> reflection
    flag, t, u, v, length_rs = csc2_cb(x, -y, -phi, am, b2)
    if flag:
        num = 46
        cytpe = path_dict[num]
        paths = set_path(paths, [t, u, MPIDIV2, v], cytpe, length_rs, num)

    # num 47     : ["L-", "S-", "L-", "R+"], arguments -- csc2_cb(-x,y,-phi,am,b2)       -> time flip
    flag, t, u, v, length_rs = csc2_cb(-x, y, -phi, am, b2)
    if flag:
        num = 47
        cytpe = path_dict[num]
        paths = set_path(paths, [-t, -u, -MPIDIV2, -v], cytpe, length_rs, num)

    # num 48     : ["R-", "S-", "R-", "L+"], arguments -- csc2_cb(-x,-y,phi,ap,b2)       -> time flip and reflection
    flag, t, u, v, length_rs = csc2_cb(-x, -y, phi, ap, b2)
    if flag:
        num = 48
        cytpe = path_dict[num]
        paths = set_path(paths, [-t, -u, -MPIDIV2, -v], cytpe, length_rs, num)

    return paths


# ------------------------------------------ CCSCC CURVES  ------------------------------------------------------

def c_c2sc2_c(x, y, phi, rs, rc):
    """
        Formula 8.11: C|C[pi/2]SC[pi/2]|C
    """
    a = x + rs
    b = y - rc
    u1, theta = polar(a, b)

    if u1 >= RADCURVMUL4:
        u = np.sqrt(u1 ** 2 - SQRADCURVMUL2) - RADCURVMUL4
        if u >= 0:
            alpha = np.arctan2(RADCURVMUL2, (u + RADCURVMUL4))
            t = wrapToPi(MPIDIV2 + theta + alpha)
            v = wrapToPi(t - phi)

            length_rs = RADCURV * (t + MPI + v) + u
            if t >= 0 and u >= 0 and v >= 0:
                return True, t, u, v, length_rs

    return False, 0.0, 0.0, 0.0, INFINITY


def CCSCC(x, y, phi, paths):
    """
        Formula 8.11: C|C[pi/2]SC[pi/2]|C
    """
    sphi, cphi, ap, am, b1, b2 = curve_params(phi)
    ########## ---------------------------------- C | C2 S C2 | C   Curves --------------------------------------------

    #  num 33     : ["L+", "R-", "S-", "L-", "R+"], arguments -- c_c2sc2_c(x,y,phi,ap,b2)       -> given
    flag, t, u, v, length_rs = c_c2sc2_c(x, y, phi, ap, b2)
    if flag:
        num = 33
        cytpe = path_dict[num]
        paths = set_path(paths, [t, MPIDIV2, u, MPIDIV2, v], cytpe, length_rs, num)

    # num 34     : ["R+", "L-", "S-", "R-", "L+"], arguments -- c_c2sc2_c(x,-y,-phi,am,b2)    -> reflection of num 33
    flag, t, u, v, length_rs = c_c2sc2_c(x, -y, -phi, am, b2)
    if flag:
        num = 34
        cytpe = path_dict[num]
        paths = set_path(paths, [t, MPIDIV2, u, MPIDIV2, v], cytpe, length_rs, num)

    # num 35     : ["L-", "R+", "S+", "L+", "R-"], arguments -- c_c2sc2_c(-x,y,-phi,am,b2)    -> time flip of num 33
    flag, t, u, v, length_rs = c_c2sc2_c(-x, y, -phi, am, b2)
    if flag:
        num = 35
        cytpe = path_dict[num]
        paths = set_path(paths, [-t, -MPIDIV2, -u, -MPIDIV2, -v], cytpe, length_rs, num)

    # num 36     : ["R-", "L+", "S+", "R+", "L-"], arguments -- c_c2sc2_c(-x,-y,phi,ap,b2)      -> time flip of num 34
    flag, t, u, v, length_rs = c_c2sc2_c(-x, -y, phi, ap, b2)
    if flag:
        num = 36
        cytpe = path_dict[num]
        paths = set_path(paths, [-t, -MPIDIV2, -u, -MPIDIV2, -v], cytpe, length_rs, num)

    return paths


### -------------------------------------- PATH PLANNING STARTS HERE ---------------------------------------------
def calc_shortest_path_length(sx, sy, syaw, gx, gy, gyaw):
    '''
        Finding the shortest path
    '''

    q0 = [sx, sy, syaw]  # starting node
    q1 = [gx, gy, gyaw]  # ending node
    paths = generate_path(q0, q1)

    minL = np.Inf
    best_path_index = -1
    for i, _ in enumerate(paths):
        if paths[i].total_length <= minL:
            minL = paths[i].total_length
            best_path_index = i

    # bpath = paths[best_path_index]

    return minL


def generate_coordinates(path, q0):
    '''
        Given the path object, returns px, py, pyaw and directions with the interpolated points
        ds = rho * dheta --> dheta = ds / rho
    :param path:
    :return:
    '''

    rho = RADCURV
    ds = MOTION_RESOLUTION

    mode = path.ctypes
    clen = path.lengths

    # angle resolution
    dth = ds / rho

    # path always start at 0, 0, 0
    x0, y0, yaw0 = 0.0, 0.0, 0.0

    px, py, pyaw, directions = [], [], [], []

    for maneuver, ll in zip(mode, clen):
        '''
            First compute the pure maneuver, than stitch 
        '''
        if maneuver == 'L+':
            ll_pos = np.abs(ll)
            t = np.arange(0, ll_pos, dth)
            t[-1] = ll_pos

            x = rho * np.sin(t)
            y = rho * (1 - np.cos(t))
            yaw = t
            dirs = np.zeros(t.shape)

            # Rotate to the world coordinates
            xb, yb = rotate_bw(x, y, yaw0)

            # pxsec = (xb + x0).tolist()
            # pysec = (yb + y0).tolist()
            # pyawsec = (yaw + yaw0).tolist()
            # pdirsec = (dirs + 1).tolist()
            #
            # ## Patch the sections
            # path.man_pathches += [pxsec, pysec, pyawsec, pdirsec]
            #
            # px += pxsec
            # py += pysec
            # pyaw += pyawsec
            # directions += pdirsec

            px += (xb + x0).tolist()
            py += (yb + y0).tolist()
            pyaw += (yaw + yaw0).tolist()
            directions += (dirs + 1).tolist()

            x0, y0, yaw0 = px[-1], py[-1], pyaw[-1]
            # plot_xy(px, py)

        if maneuver == 'S+':
            ll_pos = np.abs(ll)
            s = np.arange(0, ll_pos, ds)
            s[-1] = ll_pos

            x = s
            y = np.zeros(x.shape[0])
            yaw = np.zeros(x.shape[0])
            dirs = np.zeros(s.shape)

            xb, yb = rotate_bw(x, y, yaw0)

            px += (xb + x0).tolist()
            py += (yb + y0).tolist()
            pyaw += (yaw + yaw0).tolist()
            directions += (dirs + 1).tolist()

            x0, y0, yaw0 = px[-1], py[-1], pyaw[-1]
            # plot_xy(px, py)

        if maneuver == 'R+':
            ll_pos = np.abs(ll)
            t = np.arange(0, ll_pos, dth)

            if len(t):
                t[-1] = ll_pos
            else:
                t = np.array([0.0])

            x = rho * np.sin(t)
            y = - rho * (1 - np.cos(t))
            yaw = -t
            dirs = np.zeros(t.shape)

            # Rotate to the world coordinates
            xb, yb = rotate_bw(x, y, yaw0)

            px += (xb + x0).tolist()
            py += (yb + y0).tolist()
            pyaw += (yaw + yaw0).tolist()
            directions += (dirs + 1).tolist()

            x0, y0, yaw0 = px[-1], py[-1], pyaw[-1]
            # plot_xy(px, py)

        if maneuver == 'L-':
            ll_pos = np.abs(ll)
            t = np.arange(0, ll_pos, dth)

            if len(t):
                t[-1] = ll_pos
            else:
                t = np.array([0.0])

            x = -rho * np.sin(t)
            y = rho * (1 - np.cos(t))
            yaw = -t
            dirs = np.zeros(t.shape)

            # Rotate to the world coordinates
            xb, yb = rotate_bw(x, y, yaw0)

            px += (xb + x0).tolist()
            py += (yb + y0).tolist()
            pyaw += (yaw + yaw0).tolist()
            directions += (dirs - 1).tolist()

            x0, y0, yaw0 = px[-1], py[-1], pyaw[-1]
            # plot_xy(px, py)

        if maneuver == 'R-':
            '''
                We find the maneuver motion in positive discretization, then we adjust the sign of the motion in 
                the 
                equations 
            '''
            ll_pos = np.abs(ll)
            t = np.arange(0, ll_pos, dth)

            if len(t):
                t[-1] = ll_pos
            else:
                t = np.array([0.0])

            x = -rho * np.sin(t)
            y = -rho * (1 - np.cos(t))
            yaw = t
            dirs = np.zeros(t.shape)

            # Rotate to the world coordinates
            xb, yb = rotate_bw(x, y, yaw0)

            px += (xb + x0).tolist()
            py += (yb + y0).tolist()
            pyaw += (yaw + yaw0).tolist()
            directions += (dirs - 1).tolist()

            x0, y0, yaw0 = px[-1], py[-1], pyaw[-1]
            # plot_xy(px, py)

        if maneuver == 'S-':
            ll_pos = np.abs(ll)
            s = np.arange(0, ll_pos, ds)
            s[-1] = ll_pos

            x = -s
            y = np.zeros(x.shape[0])
            yaw = np.zeros(x.shape[0])
            dirs = np.zeros(s.shape)

            xb, yb = rotate_bw(x, y, yaw0)

            px += (xb + x0).tolist()
            py += (yb + y0).tolist()
            pyaw += (yaw + yaw0).tolist()
            directions += (dirs - 1).tolist()

            x0, y0, yaw0 = px[-1], py[-1], pyaw[-1]
            # plot_xy(px, py)

    '''
        In the codes above we computed the paths starting from 0, 0, 0
        To return a path, we need to change the path to the q0
    '''
    xq0, yq0, yawq0 = q0

    pxnp = np.array(px)
    pynp = np.array(py)
    pyawnp = np.array(pyaw)

    px, py = rotate_bw(pxnp, pynp, yawq0)
    px += xq0
    py += yq0
    pyaw = wrapToPi(pyawnp + yawq0)

    # check_yaw(px, py, pyaw)

    # SET PATH COORDINATES
    path.set_coords(px, py, pyaw)
    path.set_directions(directions)

    return px, py, pyaw, mode


def trapezoid_fit(joned_paths):
    '''
        Aim is to build piecewise a linear function
    '''

    conditions = []  #
    acc_profile = 0

    num_of_paths = len(joned_paths)
    if num_of_paths == 1:
        patch = joned_paths[0]
        s0 = patch['upper_lower'][0]
        sf = patch['upper_lower'][1]

        # compute total distance for accelerating and decelaring along the path
        t_acc_maxv = VMAX / ACCMAX  # time required to reach to the maximum velocity
        s_acc_2 = 2 * (1 / 2) * ACCMAX * t_acc_maxv ** 2  # the distance travelled during acc and deceleration

        length = patch['length']
        direction = patch['direction']

        if s_acc_2 <= length:
            cond = [0, t_acc_maxv, ACCMAX * direction]
            conditions.append(cond)

            # compute constant velocity and zero acc region
            s_remaining = length - s_acc_2
            t_constant = s_remaining / VMAX

            cond = [t_acc_maxv, t_acc_maxv + t_constant, 0]
            conditions.append(cond)

            # compute the deceleration conditions
            cond = [t_acc_maxv + t_constant, 2 * t_acc_maxv + t_constant, -ACCMAX * direction]
            conditions.append(cond)

        else:
            t_acc = np.sqrt(length / ACCMAX)
            cond = [0, t_acc, ACCMAX * direction]
            conditions.append(cond)

            cond = [t_acc, 2 * t_acc, -ACCMAX * direction]
            conditions.append(cond)


    else:

        t_last = 0
        for i in range(num_of_paths):

            patch = joned_paths[i]
            s0 = patch['upper_lower'][0]
            sf = patch['upper_lower'][1]

            length = patch['length']
            direction = patch['direction']

            # compute total distance for accelerating and decelaring along the path
            t_acc_maxv = VMAX / ACCMAX  # time required to reach to the maximum velocity
            s_acc_2 = 2 * (1 / 2) * ACCMAX * t_acc_maxv ** 2  # the distance travelled during acc and deceleration

            if s_acc_2 <= length:
                cond = [t_last, t_last + t_acc_maxv, ACCMAX * direction]
                conditions.append(cond)

                # compute constant velocity and zero acc region
                s_remaining = length - s_acc_2
                t_constant = s_remaining / VMAX

                cond = [t_last + t_acc_maxv, t_last + t_acc_maxv + t_constant, 0]
                conditions.append(cond)

                # compute the deceleration conditions
                cond = [t_last + t_acc_maxv + t_constant, t_last + 2 * t_acc_maxv + t_constant, -ACCMAX * direction]
                conditions.append(cond)
                a = 1

            else:
                t_acc = np.sqrt(length / ACCMAX)
                cond = [t_last, t_last + t_acc, ACCMAX * direction]
                conditions.append(cond)

                cond = [t_last + t_acc, t_last + 2 * t_acc, -ACCMAX * direction]
                conditions.append(cond)

            # Update last time point
            t_last = cond[1]

    fa = piecewise_func(conditions)
    t_last = conditions[-1][1]

    # # plot algorithm outcome
    # t_last = cond[1]
    # t_vector = np.linspace(0, t_last, 100)
    # plt.plot(t_vector, fa(t_vector))
    # plt.show()

    return fa, t_last


def piecewise_func(conditions):
    '''
        Creates a piewise function given the intervals and corresponding values

    '''
    vals = [vallist[-1] for vallist in conditions] + [0]
    fa = lambda x: np.piecewise(x, [(x >= suplow[0]) & (x < suplow[1]) for suplow in conditions], vals)

    return fa


def append_control(path):
    """
        Given the path objects, we compute the controls based on the geometry of the path
        we attach path two profiles
            - steering profile which is an numpy array 4 by 8, where the rows are for time, acc, velocity and distance
            - steering profile is again a numpy array which time, steering and distance travalled


    :param path:
    :return: the path whose control signals are computed and arc-length parametrized
    """

    '''
        The lengths for the curves are the angles in radians not in meters, we first compute the length of the curve 
        ll = angle * radius where radius is the curvature radius as the vehicle take the turn using the max curvature
        
        the results are first stored in the lengths_expanded then the lenghts are assigned back again. 
    '''

    lengths = path.lengths
    ctypes = path.ctypes
    lengths_expanded = [lengths[i] * RADCURV if (ctypes[i].startswith('L') or ctypes[i].startswith('R')) else lengths[i]
                        for i in range(len(lengths))]

    lengths = np.abs(lengths_expanded)
    path.lengths_expanded = lengths
    lengths_cum = np.hstack((0, lengths)).cumsum()

    '''
        Assert lengths == total_length
    '''
    assert round(lengths.sum(), 1) == round(path.total_length, 1)

    '''
        We fit trapezoidal acceleration profile to the section. This requires to split the path into continuous + or - 
        sections.        
        
    '''

    signs = [1 if list(k)[1] == '+' else -1 for k in ctypes]

    '''
        Count how many sign changes in the list
    '''

    joined_patches = {}  # keep the length of the joined patches

    # # Get the values of the first maneuver
    # length_0 = lengths[0]
    # sign_0 = signs[0]
    # mane_0 = list(ctypes[0])[0]
    # steering_value_0 = path.steering_map[mane_0]  # that maps the maneuver to the maximum minimum steering

    # KEEP - THE FOLLOWING VARIABLES FOR EACH PATH
    '''
        length, 
        [path_start_s0, path_end_s1], 
        path_direction, 
        path_steering_value
    
    '''
    patch_attributes = {}
    conds_steer = []

    for k in range(len(lengths)):
        sign_of_path = signs[k]
        length = lengths[k]
        s0 = lengths_cum[k]
        s1 = lengths_cum[k + 1]

        steering_value = path.steering_map[list(ctypes[k])[0]]

        ## CREATE A PIECEWISE FUNCTION

        patch_attributes[k] = {'direction': sign_of_path, 'length': length, 'upper_lower': [s0, s1],
                               'steering_val': [steering_value]}

        conds_steer.append([s0, s1, steering_value])

        # fa = piecewise_func(s0, s1, steering_value)
        # funcs_steer.append(fa)  # one can append a lambda

    # STEERING PIECEWISE AFFINE FUNCTION OF THE PATH
    '''
        steering = fa(s=distance)
         
    '''
    # STORE STEERING CONDITIONS to transform the boundaries from dist to time
    path.steering_conditions = conds_steer
    fa = piecewise_func(conditions=conds_steer)
    path.final_distance = lengths_cum[-1]
    path.controls['steering_func'] = fa

    # to test
    # s_vector = np.linspace(0, lengths_cum[-1], 100)
    # plt.plot(s_vector, fa(s_vector))
    # plt.show()

    # NOW JOIN THE PATHES 
    joined_patches[0] = patch_attributes[0]
    current_dict_ind = 0

    for i in range(len(signs) - 1):
        if signs[i + 1] == signs[i]:
            '''
                if the next curve has the same sign with the current curve, append the next curve to this patche 
                else, create the reverse direction patche in the next condition.  
                - sign changes, add new slot in the joined patches
            '''

            # update the length of the joined paths
            joined_patches[current_dict_ind]['length'] += patch_attributes[i + 1]['length']

            # update the upper bound of the joined path
            joined_patches[current_dict_ind]['upper_lower'][-1] = patch_attributes[i + 1]['upper_lower'][-1]

            # add the steeering value of the joined patche
            joined_patches[current_dict_ind]['steering_val'].append(patch_attributes[i + 1]['steering_val'][0])


        else:
            '''
                The next path section has different direction than the current one, process it in the different 
                patche section
            '''
            current_dict_ind += 1
            joined_patches[current_dict_ind] = patch_attributes[i + 1]

        ''' 
            We partitioned and joined the paths wrt their continuity, for each path partition, now we can fit a 
            trapezoid speed profile
        '''

    # DEFINE ACCELERATION PROFILE FOR THE JOINED CURVES
    fa1, final_time = trapezoid_fit(joined_patches)

    path.controls['acceleration_func'] = fa1
    path.final_time = final_time

    # t_vector = np.linspace(0, lengths_cum[-1] / VMAX, 100)
    # plt.plot(t_vector, fa1(t_vector))
    # plt.show()


def apply_control(path):
    '''
        Apply control the path
        - acceleration and steering are stored in the path object
        - acc is a function of time,
        - steering is a function of distance

        Therefore, need to transform distance <--> time functions

    :param path:
    :return:
    '''

    acc_func = path.controls['acceleration_func']
    str_fun = path.controls['steering_func']

    final_time = path.final_time
    final_distance = path.final_distance

    t_sim = np.linspace(0, final_time, 100)
    dt = t_sim[1] - t_sim[0]
    acc_values = acc_func(t_sim)

    # PLOT ACC
    # plt.plot(t_sim, acc_values)
    # plt.show()

    velocity = np.cumsum(dt * acc_values)

    # # PLOT Velocity
    # plt.plot(t_sim, velocity)
    # plt.grid()
    # plt.show()

    ## Compute Distance
    dist_travelled = np.cumsum(np.abs(velocity * dt))

    ## Compute steering with respect to time

    steering_conditions = path.steering_conditions  # dependent on the distance [s0, s1, value]
    steering_conditions_np = np.array(steering_conditions)

    # transform distance to the time
    steering_conditions_np[:, (0, 1)] = np.interp(steering_conditions_np[:, (0, 1)], dist_travelled, t_sim)

    # re-form the steering function with respect to time
    path.controls['steering_func_time'] = piecewise_func(steering_conditions_np.tolist())
    path.controls['velocity'] = velocity

    ## plot steering signal as a function of time

    # steering_values = path.controls['steering_func_time'](t_sim)
    # plot steering values
    # plt.plot(t_sim, steering_values)
    # plt.show()

    # print(path.ctypes)


def reeds_shepp_path_planning(sx, sy, syaw, gx, gy, gyaw):
    '''
        Given start and end points, returns the best - optimal path
    :param sx:
    :param sy:
    :param syaw:
    :param gx:
    :param gy:
    :param gyaw:
    :return:
    '''
    q0 = [sx, sy, syaw]
    q1 = [gx, gy, gyaw]

    paths = generate_path(q0, q1)

    if not paths:
        #  print("No path")
        #  print(sx, sy, syaw, gx, gy, gyaw)
        return False, None, None, None, None, None

    '''
        Sort the paths list inplace wrt their lengths.
        The first
    '''
    paths.sort(key=lambda x: x.total_length, reverse=False)

    # minL = float("Inf")
    # best_path_index = -1
    # for i, _ in enumerate(paths):
    #     total_length = paths[i].total_length
    #
    #     # print the length of the path, as there are almost equal length paths which are needed to be selected in
    #     # between
    #     print(f'total_length {total_length} of maneuver {paths[i].ctypes}')
    #
    #     if total_length <= minL:
    #         minL = total_length
    #         best_path_index = i
    #
    # bpath = paths[best_path_index]

    '''
        We have sorted the paths with their increasing total_length, the shortest length path is at the beginning. In 
        the following conditions we check if the two minimum length paths are close to each other such as 20cm apart. 
        If so we check their maneuver lengths and choose the one with minimum number of maneuver. 

    '''
    condition1 = paths[1].total_length - paths[0].total_length < 0.2
    condition2 = len(paths[1].ctypes) < len(paths[0].ctypes)

    if condition1 and condition2:
        bpath = paths[1]

    else:
        bpath = paths[0]

    ''' 
        Interpolate the curves and append on the object
    '''

    generate_coordinates(bpath, q0)

    """
        Append Controls to the Path 
        Apply Controls to the Path to test
    """

    append_control(bpath)
    apply_control(bpath)

    # print('\n')
    # ## test path
    # px, py, pyaw, mode = generate_coordinates(bpath, q0)
    # print(mode)
    # plt.cla()
    #
    # #  plotting
    # plot_arrow(q0[0], q0[1], q0[2], label='start')
    # plot_arrow(q1[0], q1[1], q1[2], label='end')
    # plt.plot(px, py, label="final course " + str(mode))
    #
    # plt.legend()
    # plt.grid(True)
    # plt.axis("equal")
    #
    # # plt.xlim(min(sx, gx) - 2, max(sx, gx) + 2)
    # # plt.ylim(min(sy, gy) - 2, max(sy, gy) + 2)
    # plt.pause(0.1)
    # check_yaw(px, py, pyaw)

    return True, bpath.x, bpath.y, bpath.yaw, bpath.ctypes, bpath.total_length, bpath.directions, bpath.man_patches


def reeds_shepp_path_planningq(q0, q1):
    '''
        Given start and end points, returns the best - optimal path
        # q0 = [sx, sy, syaw]
        # q1 = [gx, gy, gyaw]
    :return:
    '''
    # q0 = [sx, sy, syaw]
    # q1 = [gx, gy, gyaw]

    paths = generate_path(q0, q1)

    if not paths:
        #  print("No path")
        #  print(sx, sy, syaw, gx, gy, gyaw)
        return False, None, None, None, None, None

    '''
        Sort the paths list inplace wrt their lengths.
        Take the first one as the shortest path
    '''
    paths.sort(key=lambda x: x.total_length, reverse=False)

    # for path in paths:
    #     print(f'total_length {path.total_length} of maneuver {path.ctypes}')

    # minL = float("Inf")
    # best_path_index = -1
    # for i, _ in enumerate(paths):
    #     total_length = paths[i].total_length
    #
    #     # print the length of the path, as there are almost equal length paths which are needed to be selected in
    #     # between
    #     print(f'total_length {total_length} of maneuver {paths[i].ctypes}')
    #
    #     if total_length <= minL:
    #         minL = total_length
    #         best_path_index = i
    #
    # bpath = paths[best_path_index]

    '''
        We have sorted the paths with their increasing total_length, the shortest length path is at the beginning. In 
        the following conditions we check if the two minimum length paths are close to each other such as 20cm apart. 
        If so we check their maneuver lengths and choose the one with minimum number of maneuver. 

    '''
    condition1 = paths[1].total_length - paths[0].total_length < 0.2
    condition2 = len(paths[1].ctypes) < len(paths[0].ctypes)

    if condition1 and condition2:
        bpath = paths[1]

    else:
        bpath = paths[0]

    ''' 
        Interpolate the curves and append on the object
    '''

    generate_coordinates(bpath, q0)

    """
        Append Controls to the Path 
        Apply Controls to the Path to test
    """

    append_control(bpath)
    apply_control(bpath)

    # print('\n')
    # ## test path
    # px, py, pyaw, mode = generate_coordinates(bpath, q0)
    # print(mode)
    # plt.cla()
    #
    # #  plotting
    # plot_arrow(q0[0], q0[1], q0[2], label='start')
    # plot_arrow(q1[0], q1[1], q1[2], label='end')
    # plt.plot(px, py, label="final course " + str(mode))
    #
    # plt.legend()
    # plt.grid(True)
    # plt.axis("equal")
    #
    # # plt.xlim(min(sx, gx) - 2, max(sx, gx) + 2)
    # # plt.ylim(min(sy, gy) - 2, max(sy, gy) + 2)
    # plt.pause(0.1)
    # check_yaw(px, py, pyaw)

    return True, bpath.x, bpath.y, bpath.yaw, bpath.ctypes, bpath.total_length, bpath.directions, bpath


### ---------------------------------------  Test---------------------------------------------


if __name__ == '__main__':

    NTEST = 15
    # file_name = '../RS_paper_figs/rs_paths_paper.pickle'

    for i in range(NTEST):
        start_x = (np.random.rand() - 0.5) * np.random.randint(2)  # [m]
        start_y = (np.random.rand() - 0.5) * np.random.randint(2)  # [m]
        start_yaw = np.deg2rad((np.random.rand() - 0.5) * 180.0)  # [rad]

        end_x = (np.random.rand() + 0.5) * np.random.randint(20)  # [m]
        end_y = (np.random.rand() + 0.5) * np.random.randint(20)  # [m]
        end_yaw = np.deg2rad((np.random.rand() + 0.5) * 180.0)  # [rad]

        # reeds_shepp_path_planning(start_x, start_y, start_yaw, end_x, end_y, end_yaw)

        flag, px, py, pyaw, mode, total_length, directions, pathces = reeds_shepp_path_planning(start_x, start_y,
                                                                                                start_yaw,
                                                                                                end_x, end_y, end_yaw)

        logs_pickle[i] = dict()
        logs_pickle[i]['px'] = px
        logs_pickle[i]['py'] = py
        logs_pickle[i]['pyaw'] = pyaw
        logs_pickle[i]['mode'] = mode
        logs_pickle[i]['patches'] = pathces

        if show_animation and flag:  # pragma: no cover
            plt.cla()
            plt.plot(px, py, label="final course " + str(mode))
            plt.xlim(-25, 25)
            plt.ylim(-25, 25)

            #  plotting
            plot_arrow(start_x, start_y, start_yaw, label='start')
            plot_arrow(end_x, end_y, end_yaw, label='end')

            plt.legend()
            plt.grid(True)
            plt.axis("equal")
            plt.pause(0.1)

            #  plt.show()

        # with open(file_name, 'wb') as handle:
        #     pickle.dump(logs_pickle, handle, protocol=pickle.HIGHEST_PROTOCOL)

        print("Test done")
