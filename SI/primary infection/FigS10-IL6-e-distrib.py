#plotting the ensemble-averaged curves of mode 1, 2, 3 and asymptomatic. (figure S3)
#call AveragePlot
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.spatial import ConvexHull
from scipy.integrate import odeint
from Equation import func
from E_Calculation import E,E1,E2,E12
from scipy.stats import pearsonr
import warnings

warnings.filterwarnings('error')

def encircle(x,y, ax, **kw):

    p = np.c_[x,y]

    hull = ConvexHull(p)

    poly = plt.Polygon(p[hull.vertices,:], **kw)

    ax.add_patch(poly)


import matplotlib.mathtext as mathtext
mathtext.FontConstantsBase.sup1 = 0.5
mathtext.FontConstantsBase.sub1 = 0.2
mathtext.FontConstantsBase.sub2 = 0.3

time=np.arange(0,80,0.1)
Types = ['Mode 1', 'Mode 2', 'Mode 3', 'Mode 4']
try:
    dfs=pd.read_csv('IL6-e.csv')
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
            points.append([e[0],e[70],logv[-1],e[-1],np.max(e),np.max(IL6),
                           np.max(logv),np.mean(e[inds])])
        points=np.array(points)
        df=pd.DataFrame({'e0':points[:,0],'e7':points[:,1],'v_final':points[:,2],'e_final':points[:,3],
                         'emax':points[:,4],'IL6':points[:,5],'vmax':points[:,6],'e_aver':points[:,7],
                         'Mode':[Types[mode-1] for i in range(points.shape[0])]})
        dfs.append(df)
    dfs=pd.concat(dfs)
    dfs.to_csv('IL6-e.csv')

dfs['vf']=10**dfs['v_final']


dfs["<e>_mapped_mode"]="Others"
dfs["emax_mapped_mode"]="Others"
e_aver_mapped_mode=np.ones(dfs.shape[0])
e_max_mapped_mode=np.ones(dfs.shape[0])

e_aver_th=0.9
e_max_th=3.6
#classification by average e over the first 4 weeks
inds1=(dfs["IL6"]<1000)&(dfs["e_aver"]>e_aver_th)
inds2=(dfs["IL6"]>1000)&(dfs["IL6"]<2000)&(dfs["e_aver"]>e_aver_th)
inds3=(dfs["IL6"]>2000)&(dfs["e_aver"]>e_aver_th)
inds4=(dfs["IL6"]>2000)&(dfs["e_aver"]<e_aver_th)
dfs.loc[inds1,"<e>_mapped_mode"]="Mode 1"
dfs.loc[inds2,"<e>_mapped_mode"]="Mode 2"
dfs.loc[inds3,"<e>_mapped_mode"]="Mode 3"
dfs.loc[inds4,"<e>_mapped_mode"]="Mode 4"


#classification by maximum e
inds1=(dfs["IL6"]<1000)&(dfs["emax"]>e_max_th)
inds2=(dfs["IL6"]>1000)&(dfs["IL6"]<2000)&(dfs["emax"]>e_max_th)
inds3=(dfs["IL6"]>2000)&(dfs["emax"]>e_max_th)
inds4=(dfs["IL6"]>2000)&(dfs["emax"]<e_max_th)
dfs.loc[inds1,"emax_mapped_mode"]="Mode 1"
dfs.loc[inds2,"emax_mapped_mode"]="Mode 2"
dfs.loc[inds3,"emax_mapped_mode"]="Mode 3"
dfs.loc[inds4,"emax_mapped_mode"]="Mode 4"

fig1,axes=plt.subplots(1,3,figsize=(9,3))
ggcolors=['#2CA02C','#1F77B4','#FF7F0E','#D62728','#808080',]

markers= {'Mode 1':'o','Mode 2':'^','Mode 3':'s','Mode 4':'v','Others':'o',}
sns.scatterplot(data=dfs,x='e_aver',y='IL6',hue='Mode', style='Mode',
                hue_order=['Mode 1', 'Mode 2', 'Mode 3', 'Mode 4', ], markers=markers,
                palette=sns.set_palette(ggcolors),ax=axes[0],s=50,alpha=0.9,legend=False)
sns.scatterplot(data=dfs,x='e_aver',y='IL6',hue='<e>_mapped_mode', style='<e>_mapped_mode',
                hue_order=['Mode 1', 'Mode 2', 'Mode 3', 'Mode 4', 'Others'], markers=markers,
                palette=sns.set_palette(ggcolors),ax=axes[1],s=50,alpha=0.9)
sns.scatterplot(data=dfs,x='emax',y='IL6',hue='emax_mapped_mode', style='emax_mapped_mode',
                hue_order=['Mode 1', 'Mode 2', 'Mode 3', 'Mode 4', 'Others'], markers=markers,
                palette=sns.set_palette(ggcolors),ax=axes[2],s=50,alpha=0.9,legend=False)
axes[0].set_xlabel('<$\epsilon$>',fontsize=12,labelpad=0)
axes[1].set_xlabel('<$\epsilon$>',fontsize=12,labelpad=0)
axes[2].set_xlabel('$\epsilon_{max}$',fontsize=12,labelpad=0)
axes[0].set_ylabel('IL-6$_{max}$(pg/mL)',labelpad=0)
labels= ['Mode 1', '2', '3', '4', 'Others']
marker=['v','o', '^', 's', 'o']
s=[plt.scatter([],[],c=ggcolors[i],label=labels[i],marker=marker[i],s=50) for i in range(0,5,1)]
axes[1].legend(handles=s,frameon=False,ncol=5,loc='lower center',bbox_to_anchor=(0.5,0.98),
          columnspacing=0.8,handlelength=0.4)
axes[0].set_xscale('log')
axes[0].set_yscale('log')
axes[1].set_xscale('log')
axes[1].set_yscale('log')
axes[2].set_xscale('log')
axes[2].set_yscale('log')
ylim=axes[0].get_ylim()
axes[0].text(x=2,y=500,s="symptom-based",ha="right")
axes[1].plot([0.9,0.9],ylim,c='k',linewidth=1)
axes[1].text(x=0.9,y=500,s="<$\epsilon$>=0.9",ha="right")
axes[2].plot([3.6,3.6],ylim,c='k',linewidth=1)
axes[2].text(x=3.6,y=500,s="$\epsilon_{max}$=3.6",ha="right")
for i in range(3):
    if i>0:
        axes[i].set_ylabel("")
    axes[i].set_ylim(ylim)
fig1.subplots_adjust(bottom=0.2)
fig1.savefig("IL6-e_distrib_class.png",dpi=300)
plt.show()