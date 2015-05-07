import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import os
import math

init_time = 0
nodes = 250
end_time = 250

step = (end_time - init_time) / nodes

pwd = os.getcwd()

N =16
Nstr = str(N)
NN = N*N

def cal(time):

    time = str(time)

    P_syn = np.loadtxt(pwd+"/"+Nstr+"x"+Nstr+"/P_t"+time+".dat")
    P_msr = np.loadtxt(pwd+"/"+Nstr+"x"+Nstr+"_msr/P_t"+time+".dat")

    u_syn = np.loadtxt(pwd+"/"+Nstr+"x"+Nstr+"/u_t"+time+".dat")
    u_msr = np.loadtxt(pwd+"/"+Nstr+"x"+Nstr+"_msr/u_t"+time+".dat")

    v_syn = np.loadtxt(pwd+"/"+Nstr+"x"+Nstr+"/v_t"+time+".dat")
    v_msr = np.loadtxt(pwd+"/"+Nstr+"x"+Nstr+"_msr/v_t"+time+".dat")

    P_error = sum((P_syn-P_msr)**2)
    #u_error = 0
    #v_error = 0
    u_error = sum((u_syn-u_msr)**2)
    v_error = sum((v_syn-v_msr)**2)


    #for i in range(NN):
        #P_error += (P_syn[i]-P_msr[i])**2
        #u_error += (u_syn[i]-u_msr[i])**2
        #v_error += (v_syn[i]-v_msr[i])**2


    return P_error, u_error, v_error

P_series = np.zeros(nodes)
u_series = np.zeros(nodes)
v_series = np.zeros(nodes)

P_SE = np.zeros(nodes)
u_SE = np.zeros(nodes)
v_SE = np.zeros(nodes)


for time in range(init_time, end_time, step):
    (P_series[time/step],u_series[time/step],v_series[time/step]) = cal(time)


#for i in xrange(0, nodes, 1):
#    for j in xrange(0, i+1, 1):
#        P_SE[i] += P_series[j] * 1.0 / (i+1)
#        u_SE[i] += u_series[j] * 1.0 / (i+1)
#        v_SE[i] += v_series[j] * 1.0 / (i+1)

t = np.arange(init_time, end_time, step)

fig, ax = plt.subplots()
ax.plot(t, P_series, 'y', label='P_series')
ax.plot(t, u_series, 'b', label='u_series')
ax.plot(t, v_series, 'k', label='v_series')

#ax.plot(t, P_SE, 'y', label='P_series')
#ax.plot(t, u_SE, 'b', label='u_series')
#ax.plot(t, v_SE, 'k', label='v_series')

legend = ax.legend(loc='upper center', shadow=True)
ax.set_yscale('log')

plt.savefig("temp.png")
#plt.show()
