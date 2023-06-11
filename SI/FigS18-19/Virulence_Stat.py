import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.patches as mpatches

ggcolors = ['#2CA02C','#1F77B4','#FF7F0E','#D62728']
markers=['o','^','s','d']
labels=['Mode 1','Mode 2','Mode 3','Mode 4']
df=pd.read_csv('Virulence_Sample.csv')
df['lgv']=np.log10(df.loc[:,'vmax'])
df['lggamma']=np.log10(df.loc[:,'gamma'])
fig,ax=plt.subplots(1,1,figsize=(5,3))
gamma_bins=np.arange(0,1.5,0.2)
indices1=df.index[df['mode'] == 1].tolist()
m1,bins=np.histogram(a=df.iloc[indices1,:].lggamma,bins=gamma_bins)
indices2=df.index[df['mode'] == 2].tolist()
m2,bins=np.histogram(a=df.iloc[indices2,:].lggamma,bins=gamma_bins)
indices3=df.index[df['mode'] == 3].tolist()
m3,bins=np.histogram(a=df.iloc[indices3,:].lggamma,bins=gamma_bins)
indices4=df.index[df['mode'] == 4].tolist()
m4,bins=np.histogram(a=df.iloc[indices4,:].lggamma,bins=gamma_bins)
m=[m1,m2,m3,m4]
y=np.zeros(len(m1))
for i in [0,1,2,3]:
    x=m[i]/(m1+m2+m3+m4)
    sns.barplot(x=bins[:-1],y=x,ax=ax,bottom=y,color=ggcolors[i])
    y+=x
patches=[mpatches.Patch(color=ggcolors[i], label=labels[i]) for i in [0,1,2,3]]
ax.legend(handles=patches,loc='upper right',ncol=4,bbox_to_anchor=(1.03,1.15),
              columnspacing=1,frameon=False,fontsize=11)
ax.set_xticklabels(['%1.1f' %i for i in np.arange(0,1.3,0.2)])
ax.set_xlabel(r'$log_{10}\gamma$')
ax.set_ylabel('Frequency')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
fig.subplots_adjust(bottom=0.15)
fig2,ax2=plt.subplots(1,1,figsize=(5,3))

indices=df.loc[:,'gamma']<1e-1
df.loc[indices,'gamma']=1e-1
indices=[indices1,indices2,indices3,indices4]
ho=[1,2,3,4,5]
indices5=df.index[df['mode'] == 5].tolist()
inds=indices4+indices1+indices2+indices3
sns.scatterplot(data=df.iloc[indices5,:],x='lggamma',y='lgv',color='#E0E0E0',ax=ax2,alpha=1,zorder=-10)#
sns.scatterplot(data=df.iloc[inds,:],x='lggamma',y='lgv',hue='mode',palette=ggcolors,ax=ax2,alpha=1)#
ax2.set_xlabel(r'$log_{10}\gamma$')
ax2.set_ylabel('log$_{10}$V$_{max}$')
ax2.grid()
scatters = [plt.scatter([], [], s=60, facecolor=ggcolors[i], label=labels[i]) for i in [0,1,2,3]]
ax2.legend(handles=scatters,loc='upper right',bbox_to_anchor=(1,1.15),ncol=4,frameon=False, handlelength=0.5,columnspacing=0.8)
fig2.subplots_adjust(bottom=0.15)
x=df.loc[:,'lggamma']
y=df.loc[:,'lgv']
indices=y>2
x=x[indices]
y=y[indices]
z1 = np.polyfit(x,y, 1)
p1 = np.poly1d(z1)
print(p1)
x2=np.arange(np.min(x),np.max(x),0.1)
y2=z1[0]*x2+z1[1]
ax2.plot(x2,y2,color='k',linewidth=2)
ax2.text(0,5,s='y=%1.1fx+%1.1f'%(z1[0],z1[1]),ha='center',va='center',color='k',bbox = dict(boxstyle='round', facecolor='white', alpha=0.5))
fig.savefig('gamma_ratio.png')
fig2.savefig('gamma_vmax.png')
fig.savefig('gamma_ratio.svg')
fig2.savefig('gamma_vmax.svg')
plt.show()