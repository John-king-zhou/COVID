#calculating & plotting e & Rt time courses
import numpy as np
from E_Calculation import E,E1,E2,E12,R0
from scipy.integrate import odeint
from Equation import func
import warnings
import matplotlib.pyplot as plt
import matplotlib.mathtext as mathtext

warnings.filterwarnings('error')


dt=0.1
time=np.arange(0,28.1,dt)
for i in range(1,5,1):
    try:
        mean_e = np.loadtxt('mean%i_e.txt'%i)
        mean_logR = np.loadtxt('mean%i_R.txt'%i)
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
            Rt=R0(results,Para)
            e_list.append(e)
            e1_list.append(e1)
            e2_list.append(e2)
            e12_list.append(e12)
            R_list.append(Rt)
        e_list=np.array(e_list)
        e1_list=np.array(e1_list)
        e2_list=np.array(e2_list)
        e12_list=np.array(e12_list)
        R_list=np.array(R_list)
        innate=e1_list
        adaptive=e2_list+e12_list
        mean_e.append(np.mean(e_list, axis=0))#algebraic mean
        mean_e.append(np.mean(innate, axis=0))
        mean_e.append(np.mean(adaptive, axis=0))
        mean_e.append(10 ** np.mean(np.log10(e_list), axis=0))#geometric mean
        mean_e=np.array(mean_e)
        #R_list[R_list<1e-4]=1e-4
        logR_list=np.log10(R_list)
        mean_logR=np.mean(logR_list,axis=0)
        np.savetxt('mean%i_e.txt'%i,mean_e)
        np.savetxt('mean%i_R.txt'%i,mean_logR)

#plotting the time course of e and Rt (figure 2B)

mathtext.FontConstantsBase.sup1 = 0.5
mathtext.FontConstantsBase.sub1 = 0.2
mathtext.FontConstantsBase.sub2 = 0.3
ggcolors=['#2CA02C','#1F77B4','#FF7F0E','#D62728',]

labels=['Mode 1','2','3','4',]
yticks=[[0,5,10],[0,2,4],[0,4,8]]
dt=0.1
time=np.arange(0,28.1,dt)

fig,ax=plt.subplots(1,1,figsize=(5,3.7))
left, bottom, width, height = [0.58, 0.58, 0.3, 0.3]
ax2 = fig.add_axes([left, bottom, width, height])
ax.grid(True,zorder=-20)
ax2.grid(True,zorder=-20)
time=np.arange(0,28.1,dt)
for mode in [1,2,3,4]:
    mean_e=np.loadtxt('mean%i_e.txt'%mode)
    mean_logR=np.loadtxt('mean%i_R.txt'%mode)
    ax.plot(time,mean_e[3,:],color=ggcolors[mode-1],linewidth=4,zorder=20,label=labels[mode-1])
    ax2.plot(time,mean_logR,color=ggcolors[mode-1],linewidth=2,zorder=20)
xlim=[0,28]
ax2.hlines(0,xlim[0],xlim[1],colors='k',lw=1,zorder=10)
ax2.set_xticks([0,7,14,21,28])
ax2.set_xlim(xlim)
ax2.set_ylabel('log$_{10}$R$_t$',labelpad=0)
ax.hlines(3.6,0,28,colors='k',lw=1.3,zorder=10)
ax.set_ylabel(r'$\epsilon$', fontsize=20, labelpad=10, rotation=0, ha='center', va='center')
ax.set_xlabel('time (days)', fontsize=13, labelpad=1)
ax.set_xlim([0,28])
ax.set_ylim([0,15])
ax.set_yticks([0,6,12])
ax.set_xticks([7,14,21,28])
ax.tick_params(labelsize=13,length=6)

patches = [ plt.plot([],[], color=ggcolors[i],linewidth=2,
            label=labels[i])  for i in [3,0,1,2]]
labels2=[labels[i] for i in [3,0,1,2]]
ax.legend( bbox_to_anchor=(0.5,1.15), markerscale=1.1, fontsize=13,
           loc='upper center', ncol=4, frameon=False, handlelength=1, columnspacing=1.3)
fig.subplots_adjust(left=0.12,right=0.92,bottom=0.13,top=0.9)
fig.savefig('E-t-geomean.png')
fig.savefig('E-t-geomean.svg')
plt.show()