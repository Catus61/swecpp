import sys
import numpy as np
import matplotlib.pyplot as plt
import os

fig = plt.figure(figsize=(5,5),dpi=72)

if (len(sys.argv) == 1 or len(sys.argv) != 3) : # 1 is just the program name
    print "run with $: python plotSWE.py N t "
    exit()

address = os.getcwd()
time = str(sys.argv[2])

Nx = int(sys.argv[1])
Ny = Nx

P = np.loadtxt(address+"/"+str(Nx)+"x"+str(Ny)+"_msr/P_t"+time+".dat").reshape(Nx,Ny)
U = np.loadtxt(address+"/"+str(Nx)+"x"+str(Ny)+"_msr/u_t"+time+".dat").reshape(Nx,Ny)
V = np.loadtxt(address+"/"+str(Nx)+"x"+str(Ny)+"_msr/v_t"+time+".dat").reshape(Nx,Ny)

#X,Y = np.mgrid[0:16,0:16]
plt.quiver(U,V,pivot="middle")
#plt.contour(P)
axs = plt.imshow(P,interpolation='bicubic')
plt.colorbar(axs,shrink=1.0)
plt.axis([-1,Nx,-1,Ny])

plt.savefig(str(Nx)+"x"+str(Ny)+"syn_t"+time+".png",dpi=None,facecolor='w',edgecolor='w',orientation='portrait',papertype=None,format=None,transparent=False,bbox_inches=None,pad_inches=0.1,frameon=None)
#plt.show()


fig2 = plt.figure(figsize=(5,5),dpi=72)

P = np.loadtxt(address+"/"+str(Nx)+"x"+str(Ny)+"_msr/P_t"+time+".dat").reshape(Nx,Ny)
U = np.loadtxt(address+"/"+str(Nx)+"x"+str(Ny)+"_msr/u_t"+time+".dat").reshape(Nx,Ny)
V = np.loadtxt(address+"/"+str(Nx)+"x"+str(Ny)+"_msr/v_t"+time+".dat").reshape(Nx,Ny)

#X,Y = np.mgrid[0:16,0:16]
plt.quiver(U,V,pivot="middle")
#plt.contour(P)
axs = plt.imshow(P,interpolation='bicubic')
plt.colorbar(axs,shrink=1.0)
plt.axis([-1,Nx,-1,Ny])

plt.savefig(str(Nx)+"x"+str(Ny)+"msr_t"+time+".png",dpi=None,facecolor='w',edgecolor='w',orientation='portrait',papertype=None,format=None,transparent=False,bbox_inches=None,pad_inches=0.1,frameon=None)
#plt.show()


