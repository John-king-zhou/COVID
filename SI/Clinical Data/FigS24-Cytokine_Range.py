import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy.spatial import ConvexHull

matplotlib.rcParams["mathtext.default"]='regular'
def encircle(x,y, ax, **kw):
    p = np.c_[x,y]
    hull = ConvexHull(p)
    poly = plt.Polygon(p[hull.vertices,:], **kw)
    ax.add_patch(poly)

ggcolors=['#1F77B4','#FF7F0E','#D62728','#323232']
flierprops = dict(marker='x', markerfacecolor='k', markersize=3,
                  linestyle='none', markeredgecolor='k')
dt=0.1
inds1=[0,1,2,6,7,9,10]
labels = ['IL-2','IL-4',
          'IL-6','IL-10',r'TNF-$\alpha$',
          r'IFN-$\gamma$','type']
CStat=['Mild','Severe','Critical']
unit=['(10$^6$/mL)','(pg/mL)','(pg/mL)','(pg/mL)','(pg/mL)','(pg/mL)','(pg/mL)',r'($\mu$g/mL)']
unit2=['10$^6$/mL','pg/mL','pg/mL','pg/mL','pg/mL','pg/mL','pg/mL',r'$\mu$g/mL']

excl = pd.ExcelFile('Time_Course.xlsx')
print(len(excl.sheet_names))
data=[]
for filename in excl.sheet_names:
    df1 = pd.read_excel(excl, sheet_name=filename, header=0, index_col=None)
    type = df1.loc[0, 'Clinical Condition']
    indices=[10,11,12,13,14,15]
    dt=df1.iloc[1:,indices]
    dt=list(np.max(dt,axis=0))+[type,]
    data.append(dt)
data=pd.DataFrame(data=data,columns=labels)
print(data)
fig,ax=plt.subplots(2,3,figsize=(5,4))
ax=ax.flat
for i in range(len(labels)-1):
    sns.stripplot(data=data, x='type', y=labels[i], ax=ax[i], palette=ggcolors,zorder=-10)
    boxprops=dict(linestyle='-',color='k',facecolor='None',linewidth=1)
    sns.boxplot(data=data,x='type',y=labels[i],ax=ax[i],palette=ggcolors,showfliers=False,boxprops=boxprops,zorder=10)
    ax[i].set_xticks([0,1,2])
    ax[i].set_xticklabels([])
    ax[i].set_xlabel('')
    ax[i].set_ylabel('')
    ax[i].set_title(labels[i],fontsize=11)
    sFormatter1=matplotlib.ticker.ScalarFormatter(useOffset=False, useMathText=True)
    sFormatter1.set_powerlimits((0, 2))
    ax[i].yaxis.set_major_formatter(sFormatter1)
    ax[i].yaxis.offsetText.set_fontsize(9)
    ax[i].yaxis.offsetText.set_position((-0.12, 1))
ax[3].set_xticklabels(CStat,rotation=45)
fig.subplots_adjust(wspace=0.4,hspace=0.6,left=0.13,right=0.95,bottom=0.1,top=0.9)
ax[2].set_yticks(np.arange(0,6000,2000))
ax[3].set_yticks(np.arange(0,60,20))
ax[4].set_yticks(np.arange(0,70,20))
fig.text(x=0.95, y=0.12, s='Unit: pg/mL', ha="right", fontsize=10)
fig.subplots_adjust(bottom=0.2,hspace=0.2)
fig.savefig('ctk_range.png',dpi=300)
plt.show()