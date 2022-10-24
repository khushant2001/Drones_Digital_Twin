"""
mavDynamics 
    - this file implements the dynamic equations of motion for MAV
    - use unit quaternion for the attitude state
    
part of mavPySim 
    - Beard & McLain, PUP, 2012
    - Update history:  
        12/20/2018 - RWB
"""
import sys
sys.path.append('..')
import numpy as np
import math

# load message types
from message_types.msg_state import MsgState

import parameters.aerosonde_parameters as MAV
from tools.rotations import Quaternion2Rotation, Quaternion2Euler

rho = MAV.rho
S = MAV.S_wing
e = MAV.e
g = MAV.gravity
mass = MAV.mass
b = MAV.b
c = MAV.c

class MavDynamics:
    def __init__(self, Ts):
        self._ts_simulation = Ts
        # set initial states based on parameter file
        # _state is the 13x1 internal state of the aircraft that is being propagated:
        # _state = [pn, pe, pd, u, v, w, e0, e1, e2, e3, p, q, r]
        # We will also need a variety of other elements that are functions of the _state and the wind.
        # self.true_state is a 19x1 vector that is estimated and used by the autopilot to control the aircraft:
        # true_state = [pn, pe, h, Va, alpha, beta, phi, theta, chi, p, q, r, Vg, wn, we, psi, gyro_bx, gyro_by, gyro_bz]
        self._state = np.array([[MAV.north0],  # (0)
                               [MAV.east0],   # (1)
                               [MAV.down0],   # (2)
                               [MAV.u0],    # (3)
                               [MAV.v0],    # (4)
                               [MAV.w0],    # (5)
                               [MAV.e0],    # (6)
                               [MAV.e1],    # (7)
                               [MAV.e2],    # (8)
                               [MAV.e3],    # (9)
                               [MAV.p0],    # (10)
                               [MAV.q0],    # (11)
                               [MAV.r0]])   # (12)
        # store wind data for fast recall since it is used at various points in simulation
        self._wind = np.array([[0.], [0.], [0.]])  # wind in NED frame in meters/sec
        self._update_velocity_data()
        # store forces to avoid recalculation in the sensors function
        self._forces = np.array([[0.], [0.], [0.]])
        self._Va = MAV.u0
        self._alpha = 0
        self._beta = 0
        # initialize true_state message
        self.true_state = MsgState()

    ###################################
    # public functions
    def update(self, delta, wind):
        """
            Integrate the differential equations defining dynamics, update sensors
            delta = (delta_a, delta_e, delta_r, delta_t) are the control inputs
            wind is the wind vector in inertial coordinates
            Ts is the time step between function calls.
        """
        # get forces and moments acting on rigid bod
        forces_moments = self._forces_moments(delta)

        # Integrate ODE using Runge-Kutta RK4 algorithm
        time_step = self._ts_simulation
        k1 = self._derivatives(self._state, forces_moments)
        k2 = self._derivatives(self._state + time_step/2.*k1, forces_moments)
        k3 = self._derivatives(self._state + time_step/2.*k2, forces_moments)
        k4 = self._derivatives(self._state + time_step*k3, forces_moments)
        self._state += time_step/6 * (k1 + 2*k2 + 2*k3 + k4)

        # normalize the quaternion
        e0 = self._state.item(6)
        e1 = self._state.item(7)
        e2 = self._state.item(8)
        e3 = self._state.item(9)
        normE = np.sqrt(e0**2+e1**2+e2**2+e3**2)
        self._state[6][0] = self._state.item(6)/normE
        self._state[7][0] = self._state.item(7)/normE
        self._state[8][0] = self._state.item(8)/normE
        self._state[9][0] = self._state.item(9)/normE

        # update the airspeed, angle of attack, and side slip angles using new state
        self._update_velocity_data(wind)

        # update the message class for the true state
        self._update_true_state()

    def external_set_state(self, new_state):
        self._state = new_state

    ###################################
    # private functions
    def _derivatives(self, state, forces_moments):
        """
        for the dynamics xdot = f(x, u), returns f(x, u)
        """
        # extract the states
        # north = state.item(0)
        # east = state.item(1)
        # down = state.item(2)
        u = state.item(3)
        v = state.item(4)
        w = state.item(5)
        e0 = state.item(6)
        e1 = state.item(7)
        e2 = state.item(8)
        e3 = state.item(9)
        p = state.item(10)
        q = state.item(11)
        r = state.item(12)
        #   extract forces/moments
        fx = forces_moments.item(0)
        fy = forces_moments.item(1)
        fz = forces_moments.item(2)
        l = forces_moments.item(3)
        m = forces_moments.item(4)
        n = forces_moments.item(5)

        # position kinematics
        pos_dot = np.array([
            [(e1 ** 2) + (e0 ** 2) - (e2 ** 2) - (e3 ** 2), 2 * ((e1 * e2) - (e3 * e0)), 2 * ((e1 * e3) + (e2 * e0))],
            [2 * ((e1 * e2) + (e3 * e0)), (e2 ** 2) + (e0 ** 2) - (e1 ** 2) - (e3 ** 2), 2 * ((e2 * e3) - (e1 * e0))],
            [2 * ((e1 * e3) - (e2 * e0)), 2 * ((e2 * e3) + (e1 * e0)), (e3 ** 2) + (e0 ** 2) - (e1 ** 2) - (e2 ** 2)]
        ]) @ np.array([[u], [v], [w]])

        north_dot = pos_dot.item(0)
        east_dot = pos_dot.item(1)
        down_dot = pos_dot.item(2)

        # position dynamics
        vel = np.array([[r * v - q * w],
                        [p * w - r * u],
                        [q * u - p * v]]) + (1 / MAV.mass) * np.array([[fx], [fy], [fz]])
        u_dot = vel.item(0)
        v_dot = vel.item(1)
        w_dot = vel.item(2)

        # rotational kinematics
        e_vel = .5 * np.array([
            [0, -p, -q, -r],
            [p, 0, r, -q],
            [q, -r, 0, p],
            [r, q, -p, 0]
        ]) @ np.array([
            [e0],
            [e1],
            [e2],
            [e3]
        ])
        e0_dot = e_vel.item(0)
        e1_dot = e_vel.item(1)
        e2_dot = e_vel.item(2)
        e3_dot = e_vel.item(3)

        jx = MAV.Jx
        jy = MAV.Jy
        jz = MAV.Jz
        jxy = 0.
        jyx = 0.
        jxz = MAV.Jxz
        jzx = jxz
        jyz = 0.
        jzy = 0.
        r0 = jx * jz - jxz ** 2
        r1 = (jxz * (jx - jy + jz)) / r0
        r2 = (jz * (jz - jy) + jxz ** 2) / r0
        r3 = jz / r0
        r4 = jxz / r0
        r5 = (jz - jx) / jy
        r6 = jxz / jy
        r7 = ((jx - jy) * jx + jxz ** 2) / r0
        r8 = jx / r0

        # rotatonal dynamics
        p_dot = r1 * p * q - r2 * q * r + r3 * l + r4 * n
        q_dot = r5 * p * r - r6 * (p ** 2 - r ** 2) + m / jy
        r_dot = r7 * p * q - r1 * q * r + r4 * l + r8 * n

        # collect the derivative of the states
        x_dot = np.array([[north_dot, east_dot, down_dot, u_dot, v_dot, w_dot,
                           e0_dot, e1_dot, e2_dot, e3_dot, p_dot, q_dot, r_dot]]).T
        return x_dot

    def _update_velocity_data(self, wind=np.zeros((6,1))):
        steady_state = wind[0:3]
        gust = wind[3:6]
        e0 = self._state.item(6)
        e1 = self._state.item(7)
        e2 = self._state.item(8)
        e3 = self._state.item(9)
        # convert wind vector from world to body frame and add gust
        wind_body_frame = Quaternion2Rotation(self._state[6:10])@steady_state + gust
        # velocity vector relative to the airmass
        v_air = self._state[3:6] - wind_body_frame
        ur = v_air.item(0)
        vr = v_air.item(1)
        wr = v_air.item(2)

        # compute airspeed
        self._Va = math.sqrt(ur**2+vr**2+wr**2)
        # compute angle of attack
        if ur == 0:
            self._alpha = 0
        else:
            self._alpha = np.arctan(wr/ur)
        # compute sideslip angle
        if self._Va == 0:
            self._beta = 0
        else:
            self._beta = np.arcsin(vr / self._Va)

    def _forces_moments(self, delta):
        """
        return the forces on the UAV based on the state, wind, and control surfaces
        :param delta: np.matrix(delta_a, delta_e, delta_r, delta_t)
        :return: Forces and Moments on the UAV np.matrix(Fx, Fy, Fz, Ml, Mn, Mm)
        """
        phi, theta, psi = Quaternion2Euler(self._state[6:10])
        p = self._state.item(10)
        q = self._state.item(11)
        r = self._state.item(12)

        # Constants required
        clo = MAV.C_L_0
        cdo = MAV.C_D_0
        cmo = MAV.C_m_0
        cl_alpha = MAV.C_L_alpha
        cd_alpha = MAV.C_D_alpha
        cm_alpha = MAV.C_m_alpha
        clq = MAV.C_L_q
        cdq = MAV.C_D_q
        cmq = MAV.C_m_q
        cl_elevator = MAV.C_L_delta_e
        cd_elevator = MAV.C_D_delta_e
        cm_elevator = MAV.C_m_delta_e
        M = MAV.M
        alpha_0 = MAV.alpha0
        epsilon = MAV.epsilon
        cdp = MAV.C_D_p
        cyo = MAV.C_Y_0
        clo = MAV.C_ell_0
        cno = MAV.C_n_0
        cy_beta = -.83
        cl_beta = MAV.C_ell_beta
        cn_beta = MAV.C_n_beta
        cyp = MAV.C_Y_p
        clp = MAV.C_ell_p
        cnp = MAV.C_n_p
        cyr = MAV.C_Y_r
        clr = MAV.C_ell_r
        cnr = MAV.C_n_r
        cy_aileron = MAV.C_Y_delta_a
        cl_aileron = MAV.C_ell_delta_a
        cn_aileron = MAV.C_n_delta_a
        cy_rudder = MAV.C_Y_delta_r
        cl_rudder = MAV.C_ell_delta_r
        cn_rudder = MAV.C_n_delta_r
        # compute gravitaional forces
        e0 = self._state.item(6)
        e1 = self._state.item(7)
        e2 = self._state.item(8)
        e3 = self._state.item(9)
        f_g = Quaternion2Rotation(self._state[6:10]) @ np.array([[0.], [0.], [MAV.mass * MAV.gravity]])
        c_2Va = c / (2 * self._Va)
        q = self._state.item(11)
        e_neg_M = math.exp(-M * (self._alpha - MAV.alpha0))
        e_pos_M = math.exp(M * (self._alpha + MAV.alpha0))
        b_2Va = MAV.b / (2 * self._Va)
        sigma = (1 + e_neg_M + e_pos_M) / ((1 + e_neg_M) * (1 + e_pos_M))
        CL = (1 - sigma) * (MAV.C_L_0 + MAV.C_L_alpha * self._alpha) + sigma * (
                    2 * np.sign(self._alpha) * (np.sin(self._alpha) ** 2) * np.cos(self._alpha))
        CD = MAV.C_D_p + ((MAV.C_L_0 + MAV.C_L_alpha * self._alpha) ** 2) / (np.pi * MAV.e * MAV.AR)
        # formulas
        cx_alpha = -CD * np.cos(self._alpha) + CL * np.sin(self._alpha)
        cxq_alpha = -cdq * np.cos(self._alpha) + clq * np.sin(self._alpha)
        cx_elevator_alpha = -cd_elevator * np.cos(self._alpha) + cl_elevator * np.sin(self._alpha)
        cz_alpha = -CD * np.sin(self._alpha) - CL * np.cos(self._alpha)
        czq_alpha = -cdq * np.sin(self._alpha) - clq * np.cos(self._alpha)
        cz_elevator_alpha = -cd_elevator * np.sin(self._alpha) - cl_elevator * np.cos(self._alpha)
        # compute Lift and Drag coefficients


        # compute propeller thrust and torque
        thrust_prop, torque_prop = self._motor_thrust_torque(self._Va, delta.throttle)
        thrust_prop11 = np.array([
            [thrust_prop],
            [0],
            [0]
        ])

        # compute longitudinal forces in body frame
        fx = 0.5 * rho * (self._Va ** 2) * S*(-CD * np.cos(self._alpha) - CL * np.sin(self._alpha) - c_2Va * q * (cdq * np.cos(self._alpha) + clq * np.sin(self._alpha)) - delta.elevator*(cd_elevator*np.cos(self._alpha)) - cl_elevator*np.sin(self._alpha)) +f_g.item(0) + thrust_prop11.item(0)
        fz = 0.5 * rho * (self._Va ** 2) * S*(-CD * np.sin(self._alpha) - CL * np.cos(self._alpha) - (cdp * np.sin(self._alpha)) - clq * np.cos(self._alpha) * c_2Va * q + (-cd_elevator*np.sin(self._alpha) - cl_elevator*np.cos(self._alpha)*delta.elevator)) +f_g.item(2) +thrust_prop11.item(2)

        # compute lateral forces in body frame
        fy = 0.5 * rho * (self._Va ** 2) * S*(cyo + cy_beta * self._beta + cyp * b_2Va * p + cyr * b_2Va * r + cy_aileron * delta.aileron + cy_rudder * delta.rudder) + f_g.item(1) +thrust_prop11.item(1)

        # compute logitudinal torque in body frame
        torque_prop11 = np.array([
            [torque_prop],
            [0],
            [0]
        ])
        final_moments = torque_prop11 + .5 * self._Va ** 2 * rho * S * np.array([
            [b * (clo + cl_beta * self._beta + (clp * p * b / (2 * self._Va)) + (clr * b * r / (2 * self._Va)))],
            [c * (cmo + cm_alpha * self._alpha + (cmq * c * q / (2 * self._Va)))],
            [b * (cno + cn_beta * self._beta + (cnp * b * p / (2 * self._Va)) + (cnr * b * r / (2 * self._Va)))]
        ]) + np.array([
            [b * (cl_aileron * delta.aileron + cl_rudder * delta.rudder)],
            [c * cm_elevator * delta.elevator],
            [b * (cn_aileron * delta.aileron + cn_rudder * delta.rudder)]
        ]) * .5 * rho * self._Va ** 2 * S
        My = final_moments.item(1)
        # compute lateral torques in body frame
        Mx = final_moments.item(0)
        Mz = final_moments.item(2)

        self._forces[0] = fx
        self._forces[1] = fy
        self._forces[2] = fz
        return np.array([[fx, fy, fz, Mx, My, Mz]]).T


    def _motor_thrust_torque(self, Va, delta_t):
        # compute thrust and torque due to propeller  (See addendum by McLain)
        # map delta_t throttle command(0 to 1) into motor input voltage

        v_max = 3.7 * MAV.ncells
        d_prop = MAV.D_prop
        kv = 60./(2*np.pi*MAV.KVstar)
        kq = kv
        R = MAV.R_motor
        io = MAV.i0
        cq2 = MAV.C_Q2
        cq1 = MAV.C_Q1
        cqo = MAV.C_Q0
        ct2 = MAV.C_T2
        ct1 = MAV.C_T1
        ct0 = MAV.C_T0
        V_in = delta_t * v_max

        # Angular speed of propeller
        a = rho * d_prop ** 5 * cqo / ((2 * np.pi) ** 2)
        b = (rho * d_prop ** 4 * cq1 * self._Va / (2 * np.pi)) + (kq * kv / R)
        c = (kq * io) - (kq * V_in / R) + (rho * d_prop ** 3 * cq2 * self._Va ** 2)
        Omega_p = (-b + np.sqrt(b ** 2 - 4 * a * c)) / (2.0 * a)
        J = 2 * np.pi * self._Va / (Omega_p * d_prop)
        # thrust and torque due to propeller
        ct_func = ct2 * J ** 2 + ct1 * J + ct0
        cq_func = cq2 * J ** 2 + cq1 * J + cqo
        # thrust and torque due to propeller
        thrust_prop = rho * (Omega_p / (2 * np.pi)) ** 2 * d_prop ** 4 * ct_func
        #torque_prop = kq * ((1 / R) * (V_in - kv * Omega_p) - io)
        torque_prop = rho * (d_prop ** 5 / (2 * np.pi)**2) * Omega_p * cq_func
        return thrust_prop, torque_prop

    def _update_true_state(self):
        # update the class structure for the true state:
        #   [pn, pe, h, Va, alpha, beta, phi, theta, chi, p, q, r, Vg, wn, we, psi, gyro_bx, gyro_by, gyro_bz]
        phi, theta, psi = Quaternion2Euler(self._state[6:10])
        pdot = Quaternion2Rotation(self._state[6:10]) @ self._state[3:6]
        self.true_state.north = self._state.item(0)
        self.true_state.east = self._state.item(1)
        self.true_state.altitude = -self._state.item(2)
        self.true_state.Va = self._Va
        self.true_state.alpha = self._alpha
        self.true_state.beta = self._beta
        self.true_state.phi = phi
        self.true_state.theta = theta
        self.true_state.psi = psi
        self.true_state.Vg = np.linalg.norm(pdot)
        self.true_state.gamma = np.arcsin(pdot.item(2) / self.true_state.Vg)
        self.true_state.chi = np.arctan2(pdot.item(1), pdot.item(0))
        self.true_state.p = self._state.item(10)
        self.true_state.q = self._state.item(11)
        self.true_state.r = self._state.item(12)
        self.true_state.wn = self._wind.item(0)
        self.true_state.we = self._wind.item(1)
