import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import os
import math

init_time = 0
nodes = 15
end_time = nodes

step = (end_time - init_time) / nodes

pwd = os.getcwd()

N =16
Nstr = str(N)
NN = N*N

P_data = []
P_mod = []
P_diff = np.zeros((nodes, 256))

def cal(time):

    time = str(time)

    P_syn = np.loadtxt(pwd+"/"+Nstr+"x"+Nstr+"/P_t"+time+".dat")
    P_mod.append(P_syn)
    P_msr = np.loadtxt(pwd+"/"+Nstr+"x"+Nstr+"_msr/P_t"+time+".dat")
    P_data.append(P_msr)
    
    u_syn = np.loadtxt(pwd+"/"+Nstr+"x"+Nstr+"/u_t"+time+".dat")
    u_msr = np.loadtxt(pwd+"/"+Nstr+"x"+Nstr+"_msr/u_t"+time+".dat")
    
    v_syn = np.loadtxt(pwd+"/"+Nstr+"x"+Nstr+"/v_t"+time+".dat")
    v_msr = np.loadtxt(pwd+"/"+Nstr+"x"+Nstr+"_msr/v_t"+time+".dat")

    P_error = sum((P_syn-P_msr)**2)
    u_error = sum((u_syn-u_msr)**2)
    v_error = sum((v_syn-v_msr)**2)

    return P_error, u_error, v_error

P_series = np.zeros(nodes)
u_series = np.zeros(nodes)
v_series = np.zeros(nodes)

P_SE = np.zeros(nodes)
u_SE = np.zeros(nodes)
v_SE = np.zeros(nodes)



for time in range(init_time, end_time, step):
    (P_series[time/step],u_series[time/step],v_series[time/step]) = cal(time)
    
P_data[0].item(0)

for i in range(0, nodes):
    P_diff[i] = P_data[i] - P_mod[i]
    
    
#t = np.arange(init_time, end_time, step)
#    
#fig, ax = plt.subplots()
#ax.plot(t, P_series, 'y', label='P_series')
#ax.plot(t, u_series, 'b', label='u_series')
#ax.plot(t, v_series, 'k', label='v_series')
#
#
#legend = ax.legend(loc='upper center', shadow=True)
#ax.set_yscale('log')
#
#plt.savefig("temp.png")
#plt.show()
