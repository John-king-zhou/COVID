#Fig 2D
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.integrate import odeint
from Equation import func
from E_Calculation import E
import warnings

warnings.filterwarnings('error')

import matplotlib.mathtext as mathtext
mathtext.FontConstantsBase.sup1 = 0.5
mathtext.FontConstantsBase.sub1 = 0.2
mathtext.FontConstantsBase.sub2 = 0.3

time=np.arange(0,80,0.1)
Types = ['Mode 1', 'Mode 2', 'Mode 3', 'Mode 4']
try:
    dfs=pd.read_csv('v-e.csv')
except:
    dfs=[]
    for mode in range(1,5,1):
        Paras=np.loadtxt('Spl_Para%i.txt'%mode)
        Mem=[]
        Aver_Results=[]
        N=Paras.shape[0]
        points=[]
        print('mode %i' % mode)
        for k in range(Paras.shape[0]):
            ratio=k/N
            rat_str=['>']*int(ratio*50)+['-']*(50-int(ratio*50))
            rat_str=''.join(rat_str)
            print('\r'+rat_str+'%.2f %%' %(ratio*100), end='')
            Para=Paras[k]
            initial = [0.01, 0, Para[52] / Para[62], 0, 0, Para[53] / Para[65], 0, 0, Para[159], 0, 0, 0, 0, 0, 0,
                       Para[54] / Para[70], 0, Para[160], 0, 0, 0, 0, 0, 0, Para[77] / Para[101], Para[82] / Para[102],
                       Para[84] / Para[103], (Para[88] + Para[54] / Para[70] * Para[90]) / Para[104],
                       Para[91] / Para[105],
                       Para[95] / Para[106], 0, 0]
            results=odeint(func, initial, time, args=(Para,))
            v=results[:,0]
            inds=np.arange(0,280,1)
            v[v<1e-2]=1e-2
            logv=np.log10(v)
            e=E(results, Para)
            IL6=results[:, 26]
            points.append([logv[-1],np.max(e)])
        points=np.array(points)
        df=pd.DataFrame({'v_final':points[:,0], 'emax':points[:,1],
                         'Mode':[Types[mode-1] for i in range(points.shape[0])]})
        dfs.append(df)
    dfs=pd.concat(dfs)
    dfs.to_csv('v-e.csv')

dfs['vf']=10**dfs['v_final']
fig,ax=plt.subplots(1,1,figsize=(5,2))

ggcolors=['#2CA02C','#1F77B4','#FF7F0E','#D62728',]
ho= ['Mode 1', 'Mode 2', 'Mode 3', 'Mode 4', ]
markers=['v','o','^','s']
sizes=[80,80,60,80]
for i in range(4):
    inds=(dfs['Mode']==ho[i])
    x=dfs.loc[inds,'emax']
    y=dfs.loc[inds,'vf']
    ax.scatter(x,y,marker=markers[i],facecolor=ggcolors[i],edgecolor='w',lw=0.6,s=sizes[i],alpha=0.9,)
ax.set_xlabel('$\epsilon_{max}$',fontsize=13,labelpad=-2)
ax.set_xscale('log')
ax.set_yscale('log')
xlim=ax.get_xlim()

ylim=ax.get_ylim()
y=np.linspace(ylim[0],ylim[1],100)
x=np.ones(100)*3.6
ax.plot(x,y,c='k',zorder=10,linewidth=1)
ax.annotate(xy=(3.6,0.1),xytext=(1,1),text=r'$\epsilon=\gamma$',color='k',ha='left',
             arrowprops=dict(arrowstyle='->', connectionstyle="arc3", color='k',lw=1),)
ax.set_yticks([0.01,1,100,10000])
ax.set_yticklabels(['0','$10^0$','$10^2$','$10^4$'])
ax.set_ylabel('[nCoV]$_{final}$',labelpad=0)
ax.set_ylim(ylim)
fig.subplots_adjust(bottom=0.2,top=0.9,left=0.18,right=0.9,hspace=0.3)
fig.savefig('v-e.svg')
fig.savefig('v-e.png')
plt.show()
