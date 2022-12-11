"""
observer
    - Beard & McLain, PUP, 2012
    - Last Update:
        3/2/2019 - RWB
"""
import sys
import numpy as np
from scipy import stats
sys.path.append('..')
import parameters.control_parameters as CTRL
import parameters.simulation_parameters as SIM
import parameters.sensor_parameters as SENSOR
from tools.wrap import wrap
from message_types.msg_state import MsgState
from message_types.msg_sensors import MsgSensors
import parameters.aerosonde_parameters as MAV
from math import tan
from math import cos
from math import sin

class Observer:
    def __init__(self, ts_control, initial_state = MsgState(), initial_measurements = MsgSensors()):
        # initialized estimated state message
        self.estimated_state = initial_state
        # use alpha filters to low pass filter gyros and accels
        # alpha = Ts/(Ts + tau) where tau is the LPF time constant
        self.lpf_gyro_x = AlphaFilter(alpha=0.7, y0=initial_measurements.gyro_x)
        self.lpf_gyro_y = AlphaFilter(alpha=0.7, y0=initial_measurements.gyro_y)
        self.lpf_gyro_z = AlphaFilter(alpha=0.7, y0=initial_measurements.gyro_z)
        self.lpf_accel_x = AlphaFilter(alpha=0.7, y0=initial_measurements.accel_x)
        self.lpf_accel_y = AlphaFilter(alpha=0.7, y0=initial_measurements.accel_y)
        self.lpf_accel_z = AlphaFilter(alpha=0.7, y0=initial_measurements.accel_z)
        # use alpha filters to low pass filter absolute and differential pressure
        self.lpf_abs = AlphaFilter(alpha=0.9, y0=initial_measurements.abs_pressure)
        self.lpf_diff = AlphaFilter(alpha=0.7, y0=initial_measurements.diff_pressure)
        # ekf for phi and theta
        self.attitude_ekf = EkfAttitude()
        # ekf for pn, pe, Vg, chi, wn, we, psi
        self.position_ekf = EkfPosition()


    def update(self, measurement):

        # estimates for p, q, r are low pass filter of gyro minus bias estimate
        self.estimated_state.p = self.lpf_gyro_x.update(measurement.gyro_x) - SENSOR.gyro_x_bias
        self.estimated_state.q = self.lpf_gyro_y.update(measurement.gyro_y) - SENSOR.gyro_y_bias
        self.estimated_state.r = self.lpf_gyro_z.update(measurement.gyro_x) - SENSOR.gyro_z_bias

        # invert sensor model to get altitude and airspeed
        self.estimated_state.altitude = self.lpf_abs.update(measurement.static_pressure)/(MAV.rho*MAV.gravity)
        self.estimated_state.Va = np.sqrt(2/MAV.rho*self.lpf_diff.update(measurement.diff_pressure))

        # estimate phi and theta with simple ekf
        self.attitude_ekf.update(measurement, self.estimated_state)

        # estimate pn, pe, Vg, chi, wn, we, psi
        self.position_ekf.update(measurement, self.estimated_state)

        # not estimating these
        self.estimated_state.alpha = self.estimated_state.theta
        self.estimated_state.beta = 0.0
        self.estimated_state.bx = 0.0
        self.estimated_state.by = 0.0
        self.estimated_state.bz = 0.0
        return self.estimated_state


class AlphaFilter:
    # alpha filter implements a simple low pass filter
    # y[k] = alpha * y[k-1] + (1-alpha) * u[k]
    def __init__(self, alpha=0.5, y0=0.0):
        self.alpha = alpha  # filter parameter
        self.y = y0  # initial condition

    def update(self, u):
        self.y = self.alpha*self.y+(1-self.aplha)*u
        return self.y


class EkfAttitude:
    # implement continous-discrete EKF to estimate roll and pitch angles
    def __init__(self):
        self.Q = np.diag((1e-9, 1e-9))
        self.Q_gyro = np.eye(3) * SENSOR.gyro_sigma ** 2
        self.R_accel = np.eye(3) * SENSOR.accel_sigma ** 2
        self.N = 5  # number of prediction step per sample
        self.xhat = np.vstack((MAV.phi0, MAV.theta0)) # initial state: phi, theta
        self.P = np.eye(2)
        self.Ts = SIM.ts_control/self.N
        self.gate_threshold = stats.chi2.isf()

    def update(self, measurement, state):
        self.propagate_model(measurement, state)
        self.measurement_update(measurement, state)
        state.phi = self.xhat.item(0)
        state.theta = self.xhat.item(1)

    def f(self, x, measurement, state):
        # system dynamics for propagation model: xdot = f(x, u)
        p = measurement.gyro_x
        q = measurement.gyro_y
        r = measurement.gyro_z
        phi = x.item(0)
        theta = x.item(1)
        one = p+q*sin(phi)*tan(theta)+r*cos(phi)*tan(theta)
        two = q*cos(phi)-r*sin(phi)
        f_ = np.vstack((one,two))
        return f_

    def h(self, x, measurement, state):
        # measurement model y
        p = measurement.gyro_x
        q = measurement.gyro_y
        r = measurement.gyro_z
        phi = x.item(0)
        theta = x.item(1)
        Va = state.Va
        one = q * Va * sin(theta) + MAV.gravity * sin((theta))
        two = r * Va * cos(theta) - p * Va * sin(theta) - MAV.gravity * cos(theta) * sin(phi)
        three = -q * Va * cos(theta) - MAV.gravity * cos(theta) * cos(phi)
        h_ = np.vstack((one,two,three))
        return h_

    def propagate_model(self, measurement, state):
        # model propagation
        for i in range(0, self.N):

            # propagate model
            self.xhat = self.xhat+ self.f(self.xhat, state, measurement)
            # compute Jacobian
            A = jacobian(self.f, self.xhat, state, measurement)
            #testing again
            # compute G matrix for gyro noise
            phi = state.psi
            theta = state.theta
            G = np.array([
                [1.0, sin(phi) * tan(theta), cos(phi) * tan(theta), 0],
                [0.0, cos(phi), -sin(phi), 0]
            ])
            # convert to discrete time models
            A_d = np.eye(len(A))+A*self.Ts+(A@A)*(self.Ts**2)/2.0
            # update P with discrete time model
            G_d = G*self.Ts
            self.P = A_d@self.P@A_d.transpose()+self.Ts**2*self.Q+G_d @ self.Q_gyro @ G_d.T

    def measurement_update(self, measurement, state):
        # measurement updates
        h = self.h(self.xhat, measurement, state)
        C = jacobian(self.h, self.xhat, measurement, state)
        y = np.array([[measurement.accel_x, measurement.accel_y, measurement.accel_z]]).T
        S_inv = np.linalg.inv(self.R_accel+C@self.P@C.T)
        if (y-h).T @ S_inv @ (y-h) < self.gate_threshold:
            L = self.P@C.T@S_inv
            tmp = np.eye(2)-L@C
            self.P = tmp@self.P@tmp.T+L@self.R_accel@L.T
            self.xhat =self.xhat+L@(y-h)
            #print('updating')


class EkfPosition:
    # implement continous-discrete EKF to estimate pn, pe, Vg, chi, wn, we, psi
    def __init__(self):
        self.Q = 
        self.R_gps = 
        self.R_pseudo = 
        self.N =   # number of prediction step per sample
        self.Ts = (SIM.ts_control / self.N)
        self.xhat = 
        self.P = 
        self.gps_n_old = 9999
        self.gps_e_old = 9999
        self.gps_Vg_old = 9999
        self.gps_course_old = 9999
        self.pseudo_threshold = stats.chi2.isf()
        self.gps_threshold = 100000 # don't gate GPS

    def update(self, measurement, state):
        self.propagate_model(measurement, state)
        self.measurement_update(measurement, state)
        state.north = self.xhat.item(0)
        state.east = self.xhat.item(1)
        state.Vg = self.xhat.item(2)
        state.chi = self.xhat.item(3)
        state.wn = self.xhat.item(4)
        state.we = self.xhat.item(5)
        state.psi = self.xhat.item(6)

    def f(self, x, measurement, state):
        # system dynamics for propagation model: xdot = f(x, u)
        pn = x.item(0)
        pe = x.item(1)
        Vg = x.item(2)
        chi = x.item(3)
        wn = x.item(4)
        we = x.item(5)
        psi = x.item(6)
        p = measurement.gyro_x
        r = measurement.gyro_y
        q = measurement.gyro_z
        theta = state.theta
        phi = state.phi
        pndot = Vg * cos(chi)
        pedot = Vg * sin(chi)
        Va = state.Va
        chidot = MAV.gravity/Vg*tan(phi)*cos(chi-psi)
        psidot = q * sin(phi) / cos(theta) + r * cos(phi) / cos(theta)
        Vgdot = ((Va*cos(psi) + wn) * (-Va*psidot*sin(psi)) + (Va*sin(psi) + we)*(Va*psidot*cos(psi)))*1/Vg
        f_ = np.vstack((pndot, pedot, Vgdot, chidot, 0, 0, psidot))
        return f_

    def h_gps(self, x, measurement, state):
        # measurement model for gps measurements
        pn = x.item(0)
        pe = x.item(1)
        Vg = x.item(2)
        chi = x.item(3)
        h_ = np.vstack((pn,pe,Vg,chi))
        return h_

    def h_pseudo(self, x, measurement, state):
        # measurement model for wind triangale pseudo measurement
        Vg = x.item(2)
        chi = x.item(3)
        wn = x.item(4)
        we = x.item(5)
        psi = x.item(6)
        Va = state.Va
        h_ = np.vstack((Va*cos(psi)+wn - (Vg*cos(chi))), Va*sin(psi)+we - (Vg*sin(chi)))
        return h_

    def propagate_model(self, measurement, state):
        # model propagation
        for i in range(0, self.N):
            # propagate model
            self.xhat = 
            # compute Jacobian
            A = #jacobian()
            # convert to discrete time models
            A_d = 
            # update P with discrete time model
            self.P = 

    def measurement_update(self, measurement, state):
        # always update based on wind triangle pseudu measurement
        h = self.h_pseudo(self.xhat, measurement, state)
        C = jacobian(self.h_pseudo, self.xhat, measurement, state)
        y = np.array([[0, 0]]).T
        S_inv = 
        if (y-h).T @ S_inv @ (y-h) < self.pseudo_threshold:
            L = 
            self.P = 
            self.xhat = 

        # only update GPS when one of the signals changes
        if (measurement.gps_n != self.gps_n_old) \
            or (measurement.gps_e != self.gps_e_old) \
            or (measurement.gps_Vg != self.gps_Vg_old) \
            or (measurement.gps_course != self.gps_course_old):

            h = self.h_gps(self.xhat, measurement, state)
            C = jacobian(self.h_gps, self.xhat, measurement, state)
            y_chi = wrap(measurement.gps_course, h[3, 0])
            y = np.array([[measurement.gps_n,
                           measurement.gps_e,
                           measurement.gps_Vg,
                           y_chi]]).T
            S_inv = 
            if (y-h).T @ S_inv @ (y-h) < self.gps_threshold:
                L =
                self.xhat = 
                self.P = 

            # update stored GPS signals
            self.gps_n_old = measurement.gps_n
            self.gps_e_old = measurement.gps_e
            self.gps_Vg_old = measurement.gps_Vg
            self.gps_course_old = measurement.gps_course


def jacobian(fun, x, measurement, state):
    # compute jacobian of fun with respect to x
    f = fun(x, measurement, state)
    m = f.shape[0]
    n = x.shape[0]
    eps = 0.0001  # deviation
    J = np.zeros((m, n))
    for i in range(0, n):
        x_eps = np.copy(x)
        x_eps[i][0] += eps
        f_eps = fun(x_eps, measurement, state)
        df = (f_eps - f) / eps
        J[:, i] = df[:, 0]
    return J
