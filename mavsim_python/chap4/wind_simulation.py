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
        self._steady_state = np.array([[0., -10., 0.]]).T
        #self._steady_state = np.array([[0., 5., 0.]]).T

        #   Dryden gust model parameters (section 4.4 UAV book)
        Va = 25# must set Va to a constant value
        Lu = 200.
        Lv = Lu
        Lw = 50.
        gust_flag = True
        if gust_flag==True:
            sigma_u = 1.06
            sigma_v = sigma_u
            sigma_w = 1.7
        else:
            sigma_u = 0
            sigma_v = 0
            sigma_w = 0
        s = tf('s')
        # Dryden transfer functions (section 4.4 UAV book)
        self.u_w = transferFunction(num=np.array([[0, 0, sigma_u*math.sqrt(2*Va/Lu)]]),den=np.array([[0, 1, Va/Lu]]), Ts=Ts)
        self.v_w = transferFunction(num=np.array([[0, sigma_v*math.sqrt(3*Va/Lv), sigma_v*math.sqrt(3*Va/Lv)*Va/(math.sqrt(3)*Lv)]]),den=np.array([[1, 2*Va/Lv, (Va/Lv)**2]]), Ts=Ts)
        self.w_w = transferFunction(num=np.array([[0, sigma_w*math.sqrt(3*Va/Lw), sigma_w*math.sqrt(3*Va/Lw)*Va/(math.sqrt(3)*Lw)]]),den=np.array([[1, 2*Va/Lw, (Va/Lw)**2]]), Ts=Ts)
        self._Ts = Ts

    def update(self):
        # returns a six vector.
        #   The first three elements are the steady state wind in the inertial frame
        #   The second three elements are the gust in the body frame
        #gust = np.array([[self.u_w.update(np.random.randn())],
        #                 [self.v_w.update(np.random.randn())],
        #                 [self.w_w.update(np.random.randn())]])
        gust = np.array([[0.],[0.],[0.]])
        return np.concatenate(( self._steady_state, gust ))

