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
        self.estimated_state.altitude = self.lpf_static.update(measurement.static_pressure)/(MAV.rho*MAV.gravity)
        self.estimated_state.Va = 

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
        self.Q = 
        self.Q_gyro =
        self.R_accel = 
        self.N =   # number of prediction step per sample
        self.xhat =  # initial state: phi, theta
        self.P = 
        self.Ts = 
        self.gate_threshold = #stats.chi2.isf()

    def update(self, measurement, state):
        self.propagate_model(measurement, state)
        self.measurement_update(measurement, state)
        state.phi = self.xhat.item(0)
        state.theta = self.xhat.item(1)

    def f(self, x, measurement, state):
        # system dynamics for propagation model: xdot = f(x, u)
        p = measurement[0]
        q = measurement[1]
        r = mesurement[2]
        G = np.array([[1, np.sin(x[0])*np.tan(x[1]), np.cos(x[0])*np.tan(x[1]), 0], [0, np.cos(x[0]), -np.sin(x[0]), 0]])
        f_ = np.array([[p+q*np.sin(x[0])*np.tan(x[1])+r*np.cos(x[0])*np.tan(x[1])], [q*np.cos(x[0])-r*np.sin(x[0])]])
        return f_

    def h(self, x, measurement, state):
        # measurement model y
        p = measurement[0]
        q = mesurement[1]
        r = measurement[2]
        Va = measurement[3]
        h_ = np.array([[q*Va*np.sin(x[1]) + g*np.sin(x[1])], [r*Va*np.cos(x[1])-p*Va*np.sin(x[1])-g*np.cos(x[1])*np.sin(x[0])], [-q*Va*np.cos(x[1])-g*np.cos(x[1])*np.cos(x[0])]])
        return h_

    def propagate_model(self, measurement, state):
        # model propagation
        for i in range(0, self.N):

            # propagate model
            self.xhat = 
            # compute Jacobian
            A = #jacobian()

            # compute G matrix for gyro noise
            G = 
            # convert to discrete time models
            A_d = 
            # update P with discrete time model
            self.P = 

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
        pn = x[0]
        pe = x[1]
        Vg = x[2]
        chi = x[3]
        wn = x[4]
        we = x[5]
        psi = x[6]
        Va = measurement[0]
        q = measurement[1]
        r = measurement[2]
        phi = measurement[3]
        theta = measurement[4]
        wndot = 0
        wedot = 0
        Vadot = 0
        psidot = q*(np.sin(phi)/np.cos(theta)) + r*(np.cos(phi)/np.cos(theta))
        Vgdot = ((Va*np.cos(psi)+wn)*(Vadot*np.cos(psi)-Va*psidot*np.sin(psi)+wndot)+(Va*np.sin(psi)+we)*(Vadot*np.sin(psi)+Va*psidot*np.cos(psi)+wedot))/Vg
        f_ = np.array([[Vg*np.cos(chi)], [Vg*np.sin(chi)], [Vgdot], [(g/Vg)*np.tan(phi)*np.cos(chi-psi)], [0], [0], [psidot]])
        return f_

    def h_gps(self, x, measurement, state):
        # measurement model for gps measurements
        pn = x[0]
        pe = x[1]
        Vg = x[2]
        chi = x[3]
        Va = measurement[0]
        wn = x[4]
        we = x[5]
        psi = x[6]
        h_ = np.array([[pn], [pe], [Vg], [chi], [Va*np.cos(psi)+wn-Vg*np.cos(chi)], [Va*np.sin(psi)+we-Vg*np.sin(chi)]])
        return h_

    def h_pseudo(self, x, measurement, state):
        # measurement model for wind triangale pseudo measurement
        Vg = 
        chi = 
        wn = 
        we = 
        psi = 
        h_ = 
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
