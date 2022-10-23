# -*- coding: utf-8 -*-
"""
Created on Sat Oct 22 19:06:41 2022

@author: khush
"""

"""
mavDynamics 
    - this file implements the dynamic equations of motion for MAV
    - use unit quaternion for the attitude state
    
part of mavPySim 
    - Beard & McLain, PUP, 2012
    - Update history:  
        12/20/2018 - RWB
"""

#import sys
#sys.path.append('..')
import numpy as np
import math
import matplotlib.pyplot as plt

# load message types
#from message_types.msg_state import MsgState

#import parameters.aerosonde_parameters as MAV
#from tools.rotations import Quaternion2Rotation, Quaternion2Euler
m = 11
S = .55
b = 2.9
c = .19
rho = 1.268
e = .9
aileron = 0 
elevator = 0
rudder = 0
throttle = 0
g = 9.81
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
            [(e1**2)+(e0**2) - (e2**2) - (e3**2), 2*((e1*e2)-(e3*e0)), 2*((e1*e3)+(e2*e0))],
            [2*((e1*e2)+(e3*e0)), (e2**2)+(e0**2)-(e1**2)-(e3**2), 2*((e2*e3)-(e1*e0))],
            [2*((e1*e3)-(e2*e0)), 2*((e2*e3)+(e1*e0)), (e3**2)+(e0**2)-(e1**2)-(e2**2)]
        ])@np.array([[u],[v],[w]])
         
        north_dot = pos_dot[0]
        east_dot = pos_dot[1]
        down_dot = pos_dot[2]

        # position dynamics
        vel = np.array([
                          [(r*v) - (q*w)],
                          [(p*w) - (r*u)],
                          [(q*u) - (p*v)]
                          ]) + (1/mass)*np.array([[fx], [fy], [fz]])
        u_dot = vel[0]
        v_dot = vel[1]
        w_dot = vel[2]

        # rotational kinematics
        e_vel = .5*np.array([
            [0.,-p,-q,-r],
            [p,0.,r,-q],
            [q,-r, 0., p],
            [r,q,-p,0.]
        ])@np.array([
            [e0],
            [e1],
            [e2],
            [e3]
        ])
        e0_dot = e_vel[0][0]
        e1_dot = e_vel[1][0]
        e2_dot = e_vel[2][0]
        e3_dot = e_vel[3][0]

        # rotational dynamics
        p_dot = r1*p*q-r2*q*r+r3*l+r4*n 
        q_dot = r5*p*r-r6*(p**2-r**2)+m/jy
        r_dot = r7*p*q-r1*q*r+r4*l+r8*n 

        # collect the derivative of the states
        x_dot = np.array([[north_dot, east_dot, down_dot, u_dot, v_dot, w_dot,
                           e0_dot, e1_dot, e2_dot, e3_dot, p_dot, q_dot, r_dot]]).T
        return x_dot

    def _update_velocity_data(self, wind=np.zeros((6,1))):
        steady_state = wind[0:3]
        gust = wind[3:6]
        # convert wind vector from world to body frame and add gust
        wind_body_frame = 
        # velocity vector relative to the airmass
        v_air = vel - wind_body_frame
        ur = v_air[0]
        vr = v_air[1]
        wr = v_air[2]
        # compute airspeed
        self._Va = math.sqrt(ur**2+vr**2+wr**2)
        # compute angle of attack
        if ur == 0:
            self._alpha = 90 
        else:
            self._alpha = np.atan(wr/ur)
        # compute sideslip angle
        if self._Va == 0:
            self._beta = 0
        else:
            self._beta = np.asin(vr/self._Va)

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
        
        #Constants required
        clo = .23
        cdo = .043
        cmo = .0135
        cl_alpha = 5.61
        cd_alpha = .03
        cm_alpha = -2.74
        clq = 7.95
        cdq = 0
        cmq = -38.21
        cl_elevator = .13
        cd_elevator = .0135
        cm_elevator = -.99
        m = 50
        alpha_0 = .47
        epsilon = .16
        cdp = 0
        cyo = 0
        clo = 0
        cno = 0
        cy_beta = -.83
        cl_beta = -.13
        cn_beta = .073
        cyp = 0
        clp = -.51
        cnp = -.069
        cyr = 0
        clr = .25
        cnr = -.095
        cy_aileron = .075
        cl_aileron = .17
        cn_aileron = -.011
        cy_rudder = .19
        cl_rudder = .0024
        cn_rudder = -.069
        cl_alpha_func = clo + cl_alpha*self.alpha
        cd_alpha_func = cdo +cd_alpha*self.alpha
        
        #Motor constants
        v_max = 44.4
        d_prop = .5
        kv = .0659
        kq = .0659
        R = .042
        io = 1.5
        cq2 = -.01664
        cq1 = .004970
        cqo = .005230
        ct2 = -.1079
        ct1 = -.06044
        ct0 = .09357
        # compute gravitaional forces
        f_g = np.array([
           [mass*g * 2 * (e1 * e3 - e2 * e0)],
           [mass*g * 2 * (e2 * e3 + e1 * e0)],
           [(e3**2 + e0**2 - e1**2 - e2**2) * mass*g]
           ])
        
        #formulas
        cx_alpha = -cd_alpha_func*np.cos(self.alpha)+cl_alpha_func*np.sin(self.alpha)
        cxq_alpha = -cdq*np.cos(self.alpha)+clq*np.sin(self.alpha)
        cx_elevator_alpha = -cd_elevator*np.cos(self.alpha)+cl_elevator*np.sin(self.alpha)
        cz_alpha = -cd_alpha_func*np.sin(self.alpha) - cl_alpha_func*np.cos(self.alpha)
        czq_alpha = -cdq*np.sin(self.alpha) - clq*np.cos(self.alpha)
        cz_elevator_alpha = -cd_elevator*np.sin(self.alpha) - cl_elevator*np.cos(self.alpha)
        # compute Lift and Drag coefficients
        lift = .5*rho*self._Va**2*S+np.array([
               [cx_alpha+cxq_alpha*c*q/(2.*self._Va)],
               [cyo + cy_beta*self.beta + cyp*b*p*(1/2*self._Va)+cyr*b*r/(2.*self._Va)],
               [cz_alpha+czq_alpha*c*q/(2*self._Va)]
            ])+
            np.array([
                [cx_elevator_alpha],
                [cy_aileron*aileron+cy_rudder*rudder],
                [cz_elevator_alpha*elevator]
                ])*.5*self._Va**2*rho*S

        #compute propeller thrust and torque
        thrust_prop, torque_prop = self._motor_thrust_torque(self._Va, throttle)
        thrust_prop11 = np.array([
            [thrust_prop],
            [0],
            [0]
            ])
        total_forces = f_g+lift+thrust_prop11
        # compute longitudinal forces in body frame
        fx = total_forces[0]
        fz = total_forces[2]

        # compute lateral forces in body frame
        fy = total_forces[1]

        # compute logitudinal torque in body frame
        torque_prop11 = np.array([
            [torque_prop],
            [0],
            [0]
            ])
        final_moments = torque_prop11 + .5*self._Va**2*rho*S*np.array([
            [b*(clo+cl_beta*self.beta+clp*p*b/(2*self._Va) + clr*b*r/(2*self._Va)],
            [c*(cmo + cm_alpha*self.alpha+cmq*c*q/(2*self._Va)],
            [b*(cno + cn_beta*self.beta+cnp*b*p/(2*self._Va))+ cnr*b*r/(2*self._Va)]
            ])+
            np.array([
               [b*(cl_aileron*aileron + cl_rudder*rudder)],
               [c*cm_elevator*elevator],
               [b*(cn_aileron*aileron+cn_rudder*rudder]
                ])*.5*rho*self._Va**2*S
        My = final_moments[1]
        # compute lateral torques in body frame
        Mx = final_moments[0]
        Mz = final_moments[2]

        self._forces[0] = fx
        self._forces[1] = fy
        self._forces[2] = fz
        return np.array([[fx, fy, fz, Mx, My, Mz]]).T

    def _motor_thrust_torque(self, Va, delta_t):
        # compute thrust and torque due to propeller  (See addendum by McLain)
        # map delta_t throttle command(0 to 1) into motor input voltage
        V_in = throttle*v_max

        # Angular speed of propeller
        a = rho*d_prop**5*cqo/((2*np.pi)**2)
        b = (rho*d_prop**4*cq1*self._Va/(2*np.pi))+kq*kv/R
        c = (kq*io)-(kq*V_in/R)+(rho*d_prop**3*cq2*self._Va**2)
        
        Omega_p = np.roots(a,b,c)
        J = 2*np.pi*self._Va/(Omega_p*d_prop)
        # thrust and torque due to propeller
        ct_func = ct2*J**2 + ct1*J + ct0
        cq_func = cq2*J**2 + cq1*J + cqo
        thrust_prop = rho*(Omega_p/(2*np.pi))**2*d_prop**4*ct_func
        torque_prop = kq*((1/R)*(V_in - kv*Omega_p) - io) 
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