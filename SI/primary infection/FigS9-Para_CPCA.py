#performing CPCA on the sampled parameters to get the featured dimensions responsible for the difference
#among three modes and statistics of them. (figure S2)
#
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
import pandas as pd
import xlwt
from ClassPCA import CPCA
from scipy.spatial import ConvexHull
import matplotlib.patches as mpatches


def encircle(x,y, ax, **kw):

    p = np.c_[x,y]

    hull = ConvexHull(p)

    poly = plt.Polygon(p[hull.vertices,:], **kw)

    ax.add_patch(poly)

import matplotlib.mathtext as mathtext
mathtext.FontConstantsBase.sup1 = 0.5
mathtext.FontConstantsBase.sub1 = 0.2
mathtext.FontConstantsBase.sub2 = 0.3
markers=['o', '^', 's', 'D']
bbox = dict(boxstyle="round", fc="w", color='k')
ggcolors = ['#2CA02C','#1F77B4','#FF7F0E','#D62728']

labels=['k$^{APC}_{nCoV}$','k$^{APC}_{If}$','k$_{APC}^{rcr}$','k$^{NK}_{If}$','k$^{NK}_{APC}$','k$^{Neut}_{If}$','k$^{Neut}_{D}$',
        'k$^{Neut}_{Th17}$','k$^{CD4}_{naive}$','k$^{mem}_{CD4}$','k$^{Th1}_{CD4}$','k$^{Th2}_{CD4}$','k$^{Th17}_{CD4}',
        'k$^{Tfh}_{CD4}$','k$^{iTreg}_{CD4}$','k$^{nTreg}_{APC}$','k$^{CD8}_{naive}$','k$^{CD8}_{mem}$','k$^{CTL}_{CD8}$','k^{GC}_{naive}',
        'k$^{GC}_{mem}$','k$_{PB}$','k$^c_1$','k$^c_2$','k$^c_3$','k$^c_4$','k$^k_1$','k$^k_2$','k$^k_3$','k$^k_4$','k$^k_5$',
        'K$_{ACD4}$','K$_{ACD8}$','K$_{AB}$','K$_{mem}$','naive B','CD4+T$_N$','CD8+T$_N$']
labels2=['Mode 1','Mode 2','Mode 3','Mode 4',]
para=[]
length=[0,]
indices=np.hstack((np.arange(0,25,1),np.arange(27,36,1),np.array([128,129,130,131,150,159,160])))
indices=np.delete(indices,[10,20,24])
for i in [1,2,3,4]:
    a=np.loadtxt('Spl_Para%i.txt'%(i))[:,indices]
    for j in range(len(a)):
        para.append(a[j])
    length.append(length[i-1]+len(a))
para=np.array(para)
para=np.log10(para/(np.loadtxt('Para.txt')[indices]))
c=0.5
weight=[c,c,c,c]
data,eigu=CPCA(para,4,length,weight,2,std=1)
print(eigu.shape)
data=-data
eigu=-eigu
modes=np.hstack([[np.ones(length[i+1]-length[i])*(i+1)] for i in range(4)])
data=np.vstack((data.T,modes)).T
df=pd.DataFrame(data,columns=['CPCA1','CPCA2','Mode'])
for i in range(data.shape[0]):
    df.iloc[i, 2] = 'Mode %i' % (df.iloc[i, 2])

g=sns.jointplot(data=df,x='CPCA1',y='CPCA2',hue='Mode',palette=sns.set_palette(ggcolors),kind='scatter')
ax=g.ax_joint
for i in [1,2,3,4]:
    encircle(data[length[i-1]:length[i],0],data[length[i-1]:length[i],1],ax, ec="k", fc=ggcolors[i-1],alpha=0.3, linewidth=0, zorder=-10)
inds0=[]
for i in range(len(labels)):
    eigu[i,:]=eigu[i,:]*4
    ax.annotate('', xytext=(eigu[i,0], eigu[i,1]), xy=(0, 0),
                        arrowprops=dict(arrowstyle='<-', connectionstyle="arc3", color='k',lw=1),c='k')
    if abs(eigu[i,0])>1 or abs(eigu[i,1])>1:
        inds0.append(i)
        ha='right' if eigu[i,0]<0 else 'left'
        va='top' if eigu[i,0]>0 else 'bottom'
        ax.annotate(labels[i], xy=(eigu[i, 0], eigu[i, 1]), xytext=(eigu[i, 0] * 1.1, eigu[i, 1] * 1.1),
                    c='k', bbox=bbox, ha=ha, va=va)
scatters=[plt.scatter([],[],color=ggcolors[mode-1],alpha=1,label=labels2[mode-1]) for mode in [1,2,3,4]]
ax.legend(handles=scatters,loc='lower right',bbox_to_anchor=(1.03,0),handlelength=0.8,columnspacing=0.5,fontsize=12,
          framealpha=0,markerscale=1.3)
xlim=ax.get_xlim()
ylim=ax.get_ylim()
ax.hlines(0,xlim[0],xlim[1],lw=1,zorder=-20,colors='k')
ax.vlines(0,ylim[0],ylim[1],lw=1,zorder=-20,colors='k')
ax.set_xlim(xlim)
ax.set_ylim([-3.6,2.8])
ax.set_xlabel('CPCA$_1$',fontsize=13)
ax.set_ylabel('CPCA$_2$',fontsize=13)
fig=plt.gcf()
fig.set_size_inches(5.5,5.5)
fig.savefig('Para_CPCA.svg')
fig.savefig('Para_CPCA.png')

n0=para.shape[1]
data2,eigv2,eigu2=CPCA(para,4,length,weight,n0,std=1,eigvalue=True)
workbook = xlwt.Workbook()
worksheet = workbook.add_sheet('Classification')
for i in range(n0):
    worksheet.write(0,i,'lambda%i'%i)
    worksheet.write(1,i,eigv2[i])
    worksheet.write(2,i,'eigenvector%i'%i)
    for j in range(3,3+n0):
        worksheet.write(j,i,eigu2[j-3,i])
workbook.save('Para_CPCA.xls')#output eigenvectors and eigenvalues from CPCA

def nlargest(a,N):
    return np.argsort(a)[::-1][:N]
flierprops = dict(marker='x', markerfacecolor='k', markersize=3,
                  linestyle='none', markeredgecolor='k')

standard=np.loadtxt('Para.txt')
inds=[indices[k] for k in inds0]
fig,ax=plt.subplots(nrows=1,ncols=1,figsize=(6.5,3.5))
ax.tick_params(axis='x',labelsize=11)
ax.tick_params(axis='y',labelsize=12)
dfs=[]
for j in range(4):
    Para=np.loadtxt('Spl_Para%i.txt'%(j+1))/standard
    pj=[Para[:, k] for k in inds]
    for k in range(len(pj)):
        p=pj[k]
        df=pd.DataFrame({'Mode': ['Mode %i'%(j+1) for m in range(len(p))], 'parameter': [labels[inds0[k]] for m in range(len(p))],
                         'k': np.log10(p)})
        dfs.append(df)
data1=pd.concat(dfs)
sns.boxplot(x='parameter', y='k', data=data1, hue="Mode", palette=sns.set_palette(ggcolors),
            width=0.6, linewidth=1.0, ax=ax, flierprops=flierprops)
patches=[mpatches.Patch(color=ggcolors[mode-1], label=labels2[mode-1]) for mode in [1,2,3,4]]
ax.legend(handles=patches,fontsize=12,bbox_to_anchor=(1.02,1.02),loc='upper right',ncol=4,columnspacing=0.5,framealpha=0)
ax.set_xlabel('')
ax.set_ylabel('')
fig.suptitle('Log Parameter Fold Change')
ax.set_ylim([-0.9,0.9])
ax.set_yticks([-0.8,-0.4,0,0.4,0.8])
fig.subplots_adjust(wspace=0.35, hspace=0.2,left=0.1,right=0.95)
fig.savefig('CPCAParaStat.svg')
fig.savefig('CPCAParaStat.png')

plt.show()