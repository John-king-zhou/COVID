import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.integrate import odeint
from Equation import func
from E_Calculation import *

ggcolors=['#808080','#2CA02C','#1F77B4','#FF7F0E','#D62728','#4DBEEE','#77AC30','#9467BD']

def get_result(Para):
    time0=np.arange(0, 28, 0.1)
    initial = [0.01, 0, Para[52] / Para[62], 0, 0, Para[53] / Para[65], 0, 0, Para[159], 0, 0, 0, 0, 0, 0,
               Para[54] / Para[70], 0, Para[160], 0, 0, 0, 0, 0, 0, Para[77] / Para[101], Para[82] / Para[102],
               Para[84] / Para[103], (Para[88] + Para[54] / Para[70] * Para[90]) / Para[104], Para[91] / Para[105],
               Para[95] / Para[106], 0, 0]
    results = odeint(func, initial, time0, args=(Para,))
    f_antv = f_APC_anv(results,Para)
    f_inf = f_APC_inf(results,Para)
    K_A = Para[128:132]
    results = np.vstack((results.T, f_antv, f_inf)).T
    return results,K_A

fig,axes=plt.subplots(3,4,figsize=(6,4))
fig2,axes2=plt.subplots(2,2,figsize=(6,3))
axes2=axes2.flat
time = np.arange(0, 28, 0.1)
for mode in range(1,5,1):
    ax1=axes[0,mode-1]
    ax2=axes[1,mode-1]
    ax3=axes[2, mode - 1]
    Para=np.loadtxt('Spl_Para%i.txt'%mode)#[0:10,:]
    data1=[]
    data2=[]
    data3=[]
    fCD4=[]
    fCD8=[]
    fB=[]
    for i in range(Para.shape[0]):
        para=Para[i,:]
        results,K_A=get_result(para)
        K_ACD4 = K_A[0]
        K_ACD8 = K_A[1]
        K_AB = K_A[2]
        K_mem = K_A[3]
        f_antv=results[:,32]
        f_inf=results[:,33]
        ratio=f_inf/f_antv
        fCD4.append(f_antv*Hill(f_antv,K_ACD4,2))
        fCD8.append(f_antv*Hill(f_antv,K_ACD8,2))
        fB.append(f_antv*Hill(f_antv,K_AB,2))
        data1.append(f_inf)
        data2.append(f_antv)
        data3.append(ratio)
        axes2[mode-1].plot(time,results[:,4],color=ggcolors[mode])
    mean1=np.mean(data1,axis=0)
    std1=np.std(data1,axis=0)
    mean2=np.mean(data2,axis=0)
    std2=np.std(data2,axis=0)
    mean3=np.mean(data3,axis=0)
    std3=np.std(data3,axis=0)
    ax1.plot(time, mean1, linewidth=2, c=ggcolors[mode])
    ax1.fill_between(time, mean1 - std1, mean1 + std1, facecolor=ggcolors[mode], alpha=0.2)
    ax1.set_xticks([0,7,14,21,28])
    ax1.set_xticklabels([])
    ax1.set_ylim([1,8])
    ax1.set_yticks([1,4,7])
    ax2.plot(time, mean2, linewidth=2, c=ggcolors[mode])
    ax2.fill_between(time, mean2 - std2, mean2 + std2, facecolor=ggcolors[mode], alpha=0.2)
    ax2.set_xticks([0,7,14,21,28])
    ax2.set_xticklabels([])
    ax2.set_ylim([1.1,1.9])
    ax2.set_yticks([1.2,1.6])
    ax3.plot(time, mean3, linewidth=2, c=ggcolors[mode])
    ax3.fill_between(time, mean3 - std3, mean3 + std3, facecolor=ggcolors[mode], alpha=0.2)
    ax3.set_xticks([0,7,14,21,28])
    ax3.set_xlabel('time (days)')
    ax3.set_ylim([1, 6])
    ax3.set_yticks([1,3,5])
    if mode==1:
        ax1.set_ylabel('$f_{inf}^{APC}$')
        ax2.set_ylabel('$f_{AntV}^{APC}$')
        ax3.set_ylabel('$f_{inf}^{APC}/f_{AntV}^{APC}$')
    else:
        ax1.set_yticklabels([])
        ax2.set_yticklabels([])
        ax3.set_yticklabels([])
    print('mode %i finished'%mode)
fig.savefig('APC_response.svg')
plt.show()