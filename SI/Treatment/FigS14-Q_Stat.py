import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import pandas as pd
from Equation import func
from scipy.integrate import odeint
from E_Calculation import E

def Hill(x,k,n):
    return x**n/(k**n+x**n)

def score(result,Para):
    IL6c=2000
    hc=30
    vc=1
    IL6=result[:,26]
    H=result[:,2]
    v=result[:,0]
    vf=v[-1]
    hm=np.min(H)
    IL6m=np.max(IL6)
    Q=(1+1*Hill(hm,hc,1))*(1+2*Hill(IL6c,IL6m,1))*(1+1*Hill(vc,vf,1))
    return Q

ggcolors=['#808080','#2CA02C','#1F77B4','#FF7F0E','#D62728','#4DBEEE','#77AC30','#9467BD']
labels=['Mode 1','Mode 2','Mode 3','Mode 4',]
fig,ax=plt.subplots(1,1,figsize=(7,3))
Time=np.arange(0,50,0.1)
bins=np.arange(1,9.2,0.2)
for i in [1,2,3,4]:
    print('running: mode=',i)
    Paras=np.loadtxt('Spl_Para%i.txt'%i)
    Q_list=[]
    for j in range(len(Paras)):
        ratio = j / Paras.shape[0]
        rat_str = ['>'] * int(ratio * 50) + ['-'] * (50 - int(ratio * 50))
        rat_str = ''.join(rat_str)
        print('\r' + rat_str + '%.2f %%' % (ratio * 100), end='')
        Para=Paras[j,:]
        initial = [0.01, 0, Para[52] / Para[62], 0, 0, Para[53] / Para[65], 0, 0, Para[159], 0, 0, 0, 0, 0, 0,
                   Para[54] / Para[70], 0, Para[160], 0, 0, 0, 0, 0, 0, Para[77] / Para[101], Para[82] / Para[102],
                   Para[84] / Para[103], (Para[88] + Para[54] / Para[70] * Para[90]) / Para[104], Para[91] / Para[105],
                   Para[95] / Para[106], 0, 0]
        result=odeint(func, initial, Time, args=(Para,))
        Q_list.append(score(result,Para))
    print(labels[i-1],'median Q:',np.median(Q_list),'mean Q:',np.mean(Q_list))
    sns.histplot(Q_list,ax=ax,bins=bins,color=ggcolors[i],alpha=0.7)
patches=[mpatches.Patch(color=ggcolors[i], label=labels[i-1]) for i in [1,2,3,4]]
ax.legend(handles=patches,loc='upper right',ncol=4,bbox_to_anchor=(1.03,1.18),
              columnspacing=1,frameon=False,fontsize=11)
ax.set_xlabel('Q',labelpad=2)
fig.subplots_adjust(top=0.85,bottom=0.15,left=0.1,right=0.9)
fig.savefig('Q_stat.png')
fig.savefig('Q_stat.svg')

plt.show()