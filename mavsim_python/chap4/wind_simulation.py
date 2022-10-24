"""
Class to determine wind velocity at any given moment,
calculates a steady wind speed and uses a stochastic
process to represent wind gusts. (Follows section 4.4 in uav book)
"""
import sys
sys.path.append('..')
from tools.transfer_function import transferFunction
import numpy as np
import math
from control.matlab import *


class WindSimulation:
    def __init__(self, Ts):
        # steady state wind defined in the inertial frame
        self._steady_state = np.array([[0., 0., 0.]]).T
        #self._steady_state = np.array([[0., 5., 0.]]).T

        #   Dryden gust model parameters (section 4.4 UAV book)
        Va = 10# must set Va to a constant value
        Lu = 200.
        Lv = Lu
        Lw = 50.
        gust_flag = True
        if gust_flag==True:
            sigma_u = 1.06
            sigma_v = sigma_u
            sigma_w = 1.7
        else:
            sigma_u = 1.06
            sigma_v = sigma_u
            sigma_w = 1.7
        s = tf('s')
        # Dryden transfer functions (section 4.4 UAV book)
        self._Ts = Ts
        self.u_w = tf(sigma_u*math.sqrt((2*Va)/(np.pi*Lu))/(s+Va/Lu),dt=self._Ts)
        self.v_w = tf(sigma_v*math.sqrt((3*Va)/(np.pi*Lv))*(s+(Va/Lv*3**.5))/((s+Va/Lv)**2), dt=Ts)
        self.w_w = tf(sigma_w*math.sqrt((3*Va)/(np.pi*Lw))*(s+(Va/Lw*3**.5))/((s+Va/Lw)**2), dt=Ts)

    def update(self):
        # returns a six vector.
        #   The first three elements are the steady state wind in the inertial frame
        #   The second three elements are the gust in the body frame
        #gust = np.array([[self.u_w.update(np.random.randn())],
        #                 [self.v_w.update(np.random.randn())],
        #                 [self.w_w.update(np.random.randn())]])
        gust = np.array([[0.],[0.],[0.]])
        return np.concatenate(( self._steady_state, gust ))

