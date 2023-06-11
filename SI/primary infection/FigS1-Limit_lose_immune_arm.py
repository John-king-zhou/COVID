import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def func(x0,time,para):
    nCoV=x0[0]
    If=x0[1]
    H=x0[2]

    k_infect=para[0]
    N1=para[1]
    d_If=para[2]
    r_H=para[3]
    d_H=para[4]
    e_ic=para[5]
    e_ik=para[6]
    e_c=para[7]
    e_k=para[8]


    dnCoVdt=N1*d_If*If-e_ic*nCoV-e_c*nCoV

    dIfdt=k_infect*nCoV*H-e_ik*If-e_k*If-d_If*If

    dHdt=r_H-k_infect*nCoV*H-d_H*H

    return np.array([dnCoVdt,dIfdt,dHdt])

k_infect=1.2e-4
N1=1500
d_If=0.4
r_H=2
d_H=0.04
e_ic = 0
e_ik = 0
e_c = 0
e_k = 2
para=[k_infect,N1,d_If,r_H,d_H,e_ic,e_ik,e_c,e_k]
init=np.array([0.01,0,r_H/d_H])

time=np.arange(0,100,0.1)
results=odeint(func,init,time,args=(para,))
fig,axes=plt.subplots(3,3,figsize=(6,6))
ax1=axes[0,:]
ax2=axes[1,:]
ax3=axes[2,:]
labels=["nCoV","If","H"]
for i in range(3):
    ax=ax1[i]
    ax.plot(time,results[:,i],label=labels[i],linewidth=2)
    ax.set_title(labels[i])
ax1[0].set_ylabel("Cellular Imm. Only")
ax1[0].set_yscale("log")

e_ic = 0
e_ik = 0
e_c = 2
e_k = 0
para=[k_infect,N1,d_If,r_H,d_H,e_ic,e_ik,e_c,e_k]
results=odeint(func,init,time,args=(para,))
for i in range(3):
    ax=ax2[i]
    ax.plot(time,results[:,i],label=labels[i],linewidth=2)
ax2[0].set_ylabel("Humoral Imm. Only")
ax2[0].set_yscale("log")

e_ic = 1
e_ik = 1
e_c = 0
e_k = 0
para=[k_infect,N1,d_If,r_H,d_H,e_ic,e_ik,e_c,e_k]
results=odeint(func,init,time,args=(para,))
for i in range(3):
    ax=ax3[i]
    ax.plot(time,results[:,i],label=labels[i],linewidth=2)
    ax.set_xlabel("time (days)")
ax3[0].set_ylabel("Innate Imm. Only")
ax3[0].set_yscale("log")
fig.subplots_adjust(wspace=0.3)
fig.savefig("Single_Imm_Arm.png",dpi=300)
plt.show()