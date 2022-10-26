import sys
sys.path.append('..')
from control.matlab import *
import numpy as np
import math
import matplotlib.pyplot as plt
from numpy import linalg as LA
from chap5.trim import compute_trim
import chap5.model_coef as MOD

alpha1 = MOD.a_phi1
alpha2 = MOD.a_phi2
K = 0.5

A = MOD.A_lon
B = MOD.B_lon

#Defining Transfer Function Inputs for Roll
s = tf('s')
d_phi = 1/s
alpha = 1/s

G0 = alpha2 / (s*(s + alpha1))
#G1 = G0 * (alpha + (d_phi / alpha2))
C = K
L = G0 * C
x = feedback(L,1)
y,t = step(x)

plt.plot(t,y)
plt.title('step input')
print(stepinfo(x))