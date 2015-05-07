import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import math

init_time = 200000
nodes = 80

pwd = os.getcwd()
print pwd

N = 16
Nstr = str(N)
NN = N*N

def cal(time):

    time = str(time)

    P_per = np.loadtxt(pwd+"/"+Nstr+"x"+Nstr+"_per/P_t"+time+".dat")
    P = np.loadtxt(pwd+"/"+Nstr+"x"+Nstr+"_un/P_t"+time+".dat")
    P_delta = np.loadtxt(pwd+"/"+Nstr+"x"+Nstr+"/P_t"+time+".dat")
    P_max = max(P)
    P_min = min(P)
    
    u_per = np.loadtxt(pwd+"/"+Nstr+"x"+Nstr+"_per/u_t"+time+".dat")
    u = np.loadtxt(pwd+"/"+Nstr+"x"+Nstr+"_un/u_t"+time+".dat")
    u_delta = np.loadtxt(pwd+"/"+Nstr+"x"+Nstr+"/u_t"+time+".dat")
    u_max = max(u)
    u_min = min(u)
    
    v_per = np.loadtxt(pwd+"/"+Nstr+"x"+Nstr+"_per/v_t"+time+".dat")
    v = np.loadtxt(pwd+"/"+Nstr+"x"+Nstr+"_un/v_t"+time+".dat")
    v_delta = np.loadtxt(pwd+"/"+Nstr+"x"+Nstr+"/v_t"+time+".dat")
    v_max = max(v)
    v_min = min(v)


    P_avg = np.mean(P)
    u_avg = np.mean(u)
    v_avg = np.mean(v)

    P_SD = np.std(P)
    u_SD = np.std(u)
    v_SD = np.std(v)
    

  #  delta = []

    var_delta = 0
    var_per = 0

    for i in range(NN):
        #var_delta += ((P_delta[i]-P[i])/P_SD)**2+((u_delta[i]-u[i])/u_SD)**2+((v_delta[i]-v[i])/v_SD)**2
        #var_per += ((P_per[i]-P[i])/P_SD)**2+((u_per[i]-u[i])/u_SD)**2+((v_per[i]-v[i])/v_SD)**2
        var_delta += ((P_delta[i]-P[i])/(P_max-P_min))**2+((u_delta[i]-u[i])/(u_max-u_min))**2+((v_delta[i]-v[i])/(v_max-v_min))**2
        var_per += ((P_per[i]-P[i])/(P_max-P_min))**2+((u_per[i]-u[i])/(u_max-u_min))**2+((v_per[i]-v[i])/(v_max-v_min))**2
    
    rms_delta= math.sqrt(1./NN * var_delta)
    rms_per= math.sqrt(1./NN * var_per)
#print var
#print rms

    return rms_delta, rms_per


series_delta = np.zeros(nodes)
series_per = np.zeros(nodes)

Lexpo = np.zeros(nodes-1)
L_sum = 0

end_time = init_time + 4000*nodes

for time in range(init_time,end_time,4000):
    (series_delta[(time-init_time)/4000],series_per[(time-init_time)/4000]) = cal(time)


for time in range(nodes-1):
    Lexpo[time]=1./40*math.log((series_delta[time+1]/series_per[time]))
    L_sum += Lexpo[time]		


t = np.arange(init_time/100, init_time / 100 + (nodes - 1) * 40 , 40 )
plt.plot(t, Lexpo, 'bo-')
plt.show()
