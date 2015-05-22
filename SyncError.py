import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import os
import math

init_time = 0
nodes = 1000
end_time = 10000

step = (end_time - init_time) / nodes

pwd = os.getcwd()

N =16
Nstr = str(N)
NN = N*N

def cal(time):

    time = str(time)

    P_syn = np.loadtxt(pwd+"/"+Nstr+"x"+Nstr+"/P_t"+time+".dat")
    P_max = np.max(P_syn)
    P_min = np.min(P_syn)
    P_msr = np.loadtxt(pwd+"/"+Nstr+"x"+Nstr+"_msr/P_t"+time+".dat")

    u_syn = np.loadtxt(pwd+"/"+Nstr+"x"+Nstr+"/u_t"+time+".dat")
    u_max = np.max(u_syn)
    u_min = np.min(u_syn)
    u_msr = np.loadtxt(pwd+"/"+Nstr+"x"+Nstr+"_msr/u_t"+time+".dat")

    v_syn = np.loadtxt(pwd+"/"+Nstr+"x"+Nstr+"/v_t"+time+".dat")
    v_max = np.max(v_syn)
    v_min = np.min(v_syn)
    v_msr = np.loadtxt(pwd+"/"+Nstr+"x"+Nstr+"_msr/v_t"+time+".dat")

    P_error = sum(((P_syn-P_msr)/(P_max-P_min))**2)
    u_error = sum(((u_syn-u_msr)/(u_max-u_min))**2)
    v_error = sum(((v_syn-v_msr)/(v_max-v_min))**2)

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
ax.plot(t * 1.0/100, P_series, 'y', label='P_series')
ax.plot(t * 1.0/100, u_series, 'b', label='u_series')
ax.plot(t * 1.0/100, v_series, 'k', label='v_series')

#ax.plot(t, P_SE, 'y', label='P_series')
#ax.plot(t, u_SE, 'b', label='u_series')
#ax.plot(t, v_SE, 'k', label='v_series')

legend = ax.legend(loc='upper center', shadow=True)
ax.set_yscale('log')

plt.xlabel('time (in hour)')
plt.ylabel('Synchronization Error')
plt.savefig("temp.png")
plt.show()

file = open(pwd+"/data/P.txt", "w")
file.write("t   Perror \n")
for time in range(init_time, end_time, step):
    file.write(str(t[time/step] * 1.0 / 100) + "    " + str(P_series[time/step]) + "\n")
file.close()

file = open(pwd+"/data/u.txt", "w")
file.write("t   uerror \n")
for time in range(init_time, end_time, step):
    file.write(str(t[time/step] * 1.0 / 100) + "    " + str(u_series[time/step]) + "\n")
file.close()

file = open(pwd+"/data/v.txt", "w")
file.write("t   verror \n")
for time in range(init_time, end_time, step):
    file.write(str(t[time/step] * 1.0 / 100) + "    " + str(v_series[time/step]) + "\n")
file.close()
