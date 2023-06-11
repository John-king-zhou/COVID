#calculating & plotting e & Rt time courses
import numpy as np
from E_Calculation import E,E1,E2,E12
from scipy.integrate import odeint
from Equation import func
import warnings
import matplotlib.pyplot as plt
import matplotlib.mathtext as mathtext

warnings.filterwarnings('error')

dt=0.1
time=np.arange(0,50,dt)
for i in range(1,5,1):
    try:
        mean_e = np.loadtxt('dissect_e%i.txt'%i)
    except:
        mean_e=[]
        data=np.loadtxt('Spl_Para%i.txt'%i)
        e_list=[]
        e1_list=[]
        e2_list=[]
        e12_list=[]
        R_list=[]
        s=len(data)
        print('mode %i' % (i))
        for j in range(s):
            ratio = j / s
            rat_str = ['>'] * int(ratio * 50) + ['-'] * (50 - int(ratio * 50))
            rat_str = ''.join(rat_str)
            print('\r' + rat_str + '%.2f %%' % (ratio * 100), end='')
            Para=data[j]
            initial = [0.01, 0, Para[52] / Para[62], 0, 0, Para[53] / Para[65], 0, 0, Para[159], 0, 0, 0, 0, 0, 0,
                       Para[54] / Para[70], 0, Para[160], 0, 0, 0, 0, 0, 0, Para[77] / Para[101], Para[82] / Para[102],
                       Para[84] / Para[103], (Para[88] + Para[54] / Para[70] * Para[90]) / Para[104],
                       Para[91] / Para[105], Para[95] / Para[106], 0, 0]
            results=odeint(func, initial, time, args=(Para,))
            e=E(results,Para)
            e1=E1(results,Para)
            e2=E2(results,Para)
            e12=E12(results,Para)
            e_list.append(e)
            e1_list.append(e1)
            e2_list.append(e2)
            e12_list.append(e12)
        e_list=np.array(e_list)
        e1_list=np.array(e1_list)
        e2_list=np.array(e2_list)
        e12_list=np.array(e12_list)
        mean_e.append(np.mean(e_list, axis=0))#algebraic mean
        mean_e.append(np.mean(e1_list, axis=0))
        mean_e.append(np.mean(e2_list, axis=0))
        mean_e.append(np.mean(e12_list, axis=0))
        mean_e=np.array(mean_e)
        np.savetxt('dissect_e%i.txt'%i,mean_e)

#plotting the time course of e and Rt (figure 2B)


ggcolors=['#2CA02C','#1F77B4','#FF7F0E','#D62728',]

labels=['Mode 1','2','3','4',]
dt=0.1
time=np.arange(0,50,dt)

fig,axes=plt.subplots(1,3,figsize=(6,3))

ids=[1,3,2]
for mode in [1,2,3,4]:
    mean_e=np.loadtxt('dissect_e%i.txt'%mode)[:,0:len(time)]
    for i in range(3):
        ax=axes[i]
        ax.plot(time,mean_e[ids[i],:],color=ggcolors[mode-1],linewidth=4,zorder=20,label=labels[mode-1])
subscript=["i","x","a",]
for i in range(3):
    ax=axes[i]
    ax.set_title("$\epsilon_%s$"%subscript[i])
    ax.set_xlabel('time (days)')
    ax.set_xlim([0,50])
    ax.set_ylim([0,7])
    ax.set_xticks([0,25,50])

axes[1].legend( bbox_to_anchor=(0.5,1.05), markerscale=1.1, fontsize=13,
           loc='lower center', ncol=4, frameon=False, handlelength=1, columnspacing=1.3)
fig.subplots_adjust(left=0.12,right=0.95,bottom=0.18,top=0.8,wspace=0.2)
fig.savefig("dissect_ei_ea.png",dpi=300)
plt.show()