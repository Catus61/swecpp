import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import math

init_time = 200000
nodes = 80

pwd = os.getcwd()

N = 16
Nstr = str(N)
NN = N*N

def cal(time):

    time = str(time)

    P_per = np.loadtxt(pwd+"/"+Nstr+"x"+Nstr+"_per/P_t"+time+".dat")
    P_un = np.loadtxt(pwd+"/"+Nstr+"x"+Nstr+"_un/P_t"+time+".dat")
    P_delta = np.loadtxt(pwd+"/"+Nstr+"x"+Nstr+"/P_t"+time+".dat")
    P_max = max(P_un)
    P_min = min(P_un)
    
    u_per = np.loadtxt(pwd+"/"+Nstr+"x"+Nstr+"_per/u_t"+time+".dat")
    u_un = np.loadtxt(pwd+"/"+Nstr+"x"+Nstr+"_un/u_t"+time+".dat")
    u_delta = np.loadtxt(pwd+"/"+Nstr+"x"+Nstr+"/u_t"+time+".dat")
    u_max = max(u_un)
    u_min = min(u_un)
    
    v_per = np.loadtxt(pwd+"/"+Nstr+"x"+Nstr+"_per/v_t"+time+".dat")
    v_un = np.loadtxt(pwd+"/"+Nstr+"x"+Nstr+"_un/v_t"+time+".dat")
    v_delta = np.loadtxt(pwd+"/"+Nstr+"x"+Nstr+"/v_t"+time+".dat")
    v_max = max(v_un)
    v_min = min(v_un)


    P_avg = np.mean(P_un)
    u_avg = np.mean(u_un)
    v_avg = np.mean(v_un)


    P_SD = np.std(P_un)
    u_SD = np.std(u_un)
    v_SD = np.std(v_un)



    sum_delta = 0.0
    sum_per = 0.0
    sum_P_delta = 0.0
    sum_P_per = 0.0
    sum_u_delta= 0.0
    sum_u_per = 0.0
    sum_v_delta = 0.0
    sum_v_per= 0.0
    for i in range(NN):
        sum_P_delta += ((P_delta[i]-P_un[i])/P_SD)**2
        sum_P_per += ((P_per[i]-P_un[i])/P_SD)**2
        sum_u_delta += ((u_delta[i]-u_un[i])/u_SD)**2
        sum_u_per += ((u_per[i]-u_un[i])/u_SD)**2
        sum_v_delta += ((v_delta[i]-v_un[i])/v_SD)**2
        sum_v_per += ((v_per[i]-v_un[i])/v_SD)**2
        
        sum_delta += ((P_delta[i]-P_un[i])/P_SD)**2+((u_delta[i]-u_un[i])/u_SD)**2+((v_delta[i]-v_un[i])/v_SD)**2
        sum_per += ((P_per[i]-P_un[i])/P_SD)**2+((u_per[i]-u_un[i])/u_SD)**2+((v_per[i]-v_un[i])/v_SD)**2
        
        
        ## another method
        
        #sum_delta += ((P_delta[i]-P[i])/(P_max-P_min))**2+((u_delta[i]-u[i])/(u_max-u_min))**2+((v_delta[i]-v[i])/(v_max-v_min))**2
        #sum_per += ((P_per[i]-P[i])/(P_max-P_min))**2+((u_per[i]-u[i])/(u_max-u_min))**2+((v_per[i]-v[i])/(v_max-v_min))**2
    

    return ( math.sqrt(1./NN * sum_P_delta), math.sqrt(1./NN * sum_P_per), math.sqrt(1./NN * sum_u_delta), math.sqrt(1./NN * sum_u_per), math.sqrt(1./NN * sum_v_delta), math.sqrt(1./NN * sum_v_per), math.sqrt(1./NN * sum_delta), math.sqrt(1./NN * sum_per))


series_delta = np.zeros(nodes)
series_per = np.zeros(nodes)
series_P_delta = np.zeros(nodes)
series_P_per = np.zeros(nodes)
series_u_delta = np.zeros(nodes)
series_u_per = np.zeros(nodes)
series_v_delta = np.zeros(nodes)
series_v_per = np.zeros(nodes)


Lexpo = np.zeros(nodes-1)
P_Lexpo = np.zeros(nodes-1)
u_Lexpo = np.zeros(nodes-1)
v_Lexpo = np.zeros(nodes-1)


end_time = init_time + 4000*nodes

print cal(200000)

for time in range(init_time,end_time,4000):

    ( series_P_delta[(time-init_time)/4000], series_P_per[(time-init_time)/4000], \
    series_u_delta[(time-init_time)/4000], series_u_per[(time-init_time)/4000], \
    series_v_delta[(time-init_time)/4000], series_v_per[(time-init_time)/4000], \
    series_delta[(time-init_time)/4000], series_per[(time-init_time)/4000]   )  = cal(time)


for time in range(nodes-1):
    Lexpo[time]=1./40*math.log((series_delta[time+1]/series_per[time]))
    P_Lexpo[time]=1./40*math.log((series_P_delta[time+1]/series_P_per[time]))
    u_Lexpo[time]=1./40*math.log((series_u_delta[time+1]/series_u_per[time]))
    v_Lexpo[time]=1./40*math.log((series_v_delta[time+1]/series_v_per[time]))


t = np.arange(init_time/100, init_time / 100 + (nodes - 1) * 40 , 40 )
plt.plot(t, Lexpo, 'bo-')
plt.plot(t, P_Lexpo, 'r-')
plt.plot(t, u_Lexpo, 'g-')
plt.plot(t, v_Lexpo, 'y-')
plt.show()
