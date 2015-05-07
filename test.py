# -*- coding: utf-8 -*-
"""
Created on Wed Mar 26 14:59:31 2014

@author: hippo
"""

import numpy as np
import matplotlib.pyplot as plt
import os

pwd = os.getcwd()
N = 16
Nstr = str(N)

time = 0
delay = 0
ele = 0

stime = str(time)
sdelay = str(delay)
sele = str(ele)
stotal = str(time + delay * 10)

Jac = np.loadtxt(pwd+"/"+Nstr+"x"+Nstr+"_Jac/jac_t_"+stime+".dat")
Jac_s = np.loadtxt(pwd+"/"+Nstr+"x"+Nstr+"_Jac/jac_s_t"+stime+".dat")
Jac_v = np.loadtxt(pwd+"/"+Nstr+"x"+Nstr+"_Jac/jac_v_t_"+stime+".dat")
Jac_u = np.loadtxt(pwd+"/"+Nstr+"x"+Nstr+"_Jac/jac_u_t_"+stime+".dat")

Jac_test = np.zeros((768,768))

i = 0
j = 0

for k in range(768):
    Jac_test[i][j] += Jac_s[k] * Jac_u[i][k] * Jac_v[j][k]

print Jac_test[i][j]
print Jac[i][j]