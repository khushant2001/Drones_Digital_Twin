import numpy as np
from numpy.random import standard_normal
import matplotlib.pyplot as plt
from control.matlab import *  # MATLAB-like control toolbox functionality

# DC motor parameters 
K = 4.9 # rad/s/V 
tau = 0.085 # 1/s

# Process and measurement noise variance
q  = 10**2 # covariance of process noise
Q = (q**2)*np.eye(2)

r = (1/100)*q
R  = r**2 # 1 standard dev.   covariance of measurement noise
#Making R >>Q will give a smooth graph, but Q >> R has very sharp jagged lines

# Time between measurements, Sampling time, # iters.
Tm = 0.2 # Time between measurements
nm = 150  # number of measurements
ns = 5  
Ts = Tm/ns # sampling time

# Initial conditions
x0  = np.array([0,0]) # true initial state
yhat0  = x0 + standard_normal() * q # initial measurement
xhat0 = yhat0 # initial estimated state

# Setup KF system
A = np.array([
    [0,1],
    [0,-1/tau]
    ])
B = np.array([
    [0],
    [K/tau]
    ])
C = np.array([1,0])
D = 0
sys_kf = ss(A, B, C, D)
# Setup covariance matrices
Ad =  np.exp(A*Ts) #1 + A*Ts + A**2 * Ts**2/2
P0 = R

# Setup "true" system
Bn = np.hstack((B, np.eye(2)))
Dn = np.array([0, 0, 0])
sys_true = ss(A, Bn, C, Dn)

# Init arrays for simulation
t_history = np.array([0])
x_true_history = np.array([0])
x_kf_history = np.array([xhat0[0]])
t_kf_plus = np.array([])
x_kf_plus = np.array([])
P_history = np.array([])
T1 = 0
T2 = 0
P = P0
for i in range(nm): # for all measurements
    T1 = T2
    T2 += Tm
    # Simulate "true" system for T1 to T2
    t = np.linspace(T1, T2, ns)
    u = np.sin(t) # np.ones(ns) 
    
    zeta1 = standard_normal(size = (ns,2)) @ np.sqrt(Q) #5 by 2 matrix
    zeta = zeta1[:,0]
    un = np.hstack((u, zeta, zeta1[:,1])).reshape(5,3)
    y, _, x = lsim(sys_true, un, t, x0)
    #  add noise to sensor measurement
    yn = y[-1] + standard_normal() * np.sqrt(R)
    #  update initial condition
    x0 = x[-1]
    tn = T2
    # store values
    t_history = np.hstack((t_history, np.squeeze(t)))
    x_true_history = np.append(x_true_history, np.squeeze(x)[:,0])#np.hstack((x_true_history, np.squeeze(x)[:,0])) #ERROR!!!!!
    #print(x_true_history)
    #print()
    # Kalman filter
    #  Prediction
    
    y, _, x = lsim(sys_kf, u, t, xhat0)
    for i in range(ns-1):
        P = Ad * P * Ad + Ts**2 * Q #could be a problem
    # Store vals
    x_kf_history = np.append(x_kf_history, np.squeeze(x)[:,0]) #ERROR!!!!
    #  Correction
    L = P @ C.transpose() / (R + C@P@C.transpose())
    x_plus = x[-1] + L*(yn - y[-1])
    P = (np.eye(2) - L*C)*P*(np.eye(2)-L*C).transpose() + L*R*L.transpose() #R is the covariance
    # Update IC
    xhat0 = x_plus #prediction of initial condition
    t_kf_plus = np.hstack((t_kf_plus, T2))
    x_kf_plus = np.hstack((x_kf_plus, x_plus[0]))


fig, ax = plt.subplots()
ax.plot(t_history, x_true_history, label='true')
ax.plot(t_history, x_kf_history, label='KF')
ax.plot(t_kf_plus, x_kf_plus, 'rs',label='correction')
ax.set_xlabel('t [s]')
ax.set_ylabel(r'$\omega$ [rad/s]')
ax.legend()
ax.grid()
ax.set_title('Kalaman filter')
plt.show()