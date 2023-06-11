#Using CPCA method to get the principal variables of mild/moderate, severe and critical patients
#Need Classification1.xls (generated by Classification_1.py)
#figure S6A,B
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from ClassPCA import CPCA
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
labels = ['WBC$_{max}$','Neut$_{max}$','Lymph$_{min}$',
          'Mono$_{max}$','IL-2$_{max}$','IL-4$_{max}$',
          'IL-6$_{max}$','IL-10$_{max}$',r'TNF-$\alpha$$_{max}$',
          r'IFN-$\gamma$$_{max}$','$\epsilon$*$_{max}$']
labels2 = ['WBC$_{max}$','Neut$_{max}$','Lymph$_{min}$',
          'Mono$_{max}$','IL-2$_{max}$','IL-4$_{max}$',
          'IL-6$_{max}$','IL-10$_{max}$',r'TNF-$\alpha$$_{max}$',
          r'IFN-$\gamma$$_{max}$','$\epsilon$*$_{max}$','Type']
labels3 = ['WBC','Neut','Lymph',
          'Mono','IL-2','IL-4',
          'IL-6','IL-10',r'TNF-$\alpha$$',
          r'IFN-$\gamma$','$\epsilon$*']
CStat=['Mild/Moderate','Severe','Critical']
unit=['(10$^6$/mL)','(pg/mL)','(pg/mL)','(pg/mL)','(pg/mL)','(pg/mL)','(pg/mL)',r'($\mu$g/mL)']
unit2=['10$^6$/mL','pg/mL','pg/mL','pg/mL','pg/mL','pg/mL','pg/mL',r'$\mu$g/mL']
data1=[]
length=[0,]
data0=[]
excl = pd.ExcelFile('Time_Course.xlsx')
print(len(excl.sheet_names))
for filename in excl.sheet_names:
    df1 = pd.read_excel(excl, sheet_name=filename, header=0, index_col=None)
    type = df1.loc[0, 'Clinical Condition']
    if np.isnan(type):
        continue
    data=df1.iloc[1:,3:]
    data = np.array(data).T
    data[[2,4,6],:]*=100
    data_arr=[]
    zero_array=0
    for j in range(len(labels)+3):
        if j in [2,4,6]:
            continue
        x=data[j,:]
        x=x[~np.isnan(x)]
        if len(x)==0:
            if j>=4:
                print(labels[j-3],filename,'zero array')
            else:
                print(labels[j],filename,'zero array')
            zero_array=1
            continue
        if j==3:
            data_arr.append(np.min(x))
        else:
            data_arr.append(np.max(x))
    if zero_array:
        print(filename,'continue')
        continue
    data_arr.append(type)
    data0.append(data_arr)
data0=np.array(data0)
type_arr=data0[:,-1]
data1=data0[np.argsort(type_arr),:]
length=[0,]
type_arr=data1[:,-1]
for i in range(3):
    n=len(type_arr[type_arr==i+1])
    length.append(length[i]+n)
print(length)
c=0.15
weight=[c,c,c]
data1=data1
data,eigu1=CPCA(data1[:,:-1],3,length,weight,2,std=1)
modes=np.hstack([[np.ones(length[i+1]-length[i])*(i+1)] for i in range(3)])
data=np.vstack((data.T,data1[:,-1])).T
df=pd.DataFrame(data,columns=['CPCA1','CPCA2','Type'])
df3=pd.DataFrame(data1,columns=labels2)
for i in range(data.shape[0]):
    df.iloc[i,2]=CStat[int(df.iloc[i,2])-1]
    df3.iloc[i,-1]=str(int(df3.iloc[i,-1]))
# g=sns.jointplot(data=df,x='CPCA1',y='CPCA2',hue='Type',palette=sns.set_palette(ggcolors),
#                 kind='scatter',marginal_kws={'kde_kws': {'common_norm': False}})
g = sns.JointGrid()
sns.scatterplot(data=df,x='CPCA1',y='CPCA2',hue='Type',palette=sns.set_palette(ggcolors),ax=g.ax_joint)
# sns.histplot(data=dfs,x='e',hue='type',palette=sns.set_palette(ggcolors),ax=g.ax_marg_x,
#              bins=np.linspace(0,0.5,20),multiple='dodge',linewidth=0,stat='density',common_norm=False)
sns.kdeplot(data=df,x='CPCA1',hue='Type',palette=sns.set_palette(ggcolors),ax=g.ax_marg_x,
             linewidth=1,common_norm=False)
# sns.histplot(data=dfs,y='logIL6',hue='type',palette=sns.set_palette(ggcolors),ax=g.ax_marg_y,
#              bins=np.linspace(0.2,4.2,20),multiple='dodge',linewidth=0,stat='density',common_norm=False)
sns.kdeplot(data=df,y='CPCA2',hue='Type',palette=sns.set_palette(ggcolors),ax=g.ax_marg_y,
             linewidth=1,common_norm=False)

for i in [1,2,3]:
    encircle(data[length[i-1]:length[i],0],data[length[i-1]:length[i],1],g.ax_joint, ec="k", fc=ggcolors[i-1],alpha=0.2, linewidth=0)
bbox = dict(boxstyle="round", fc="w", color=ggcolors[3])
ax=g.ax_joint
eigu=eigu1
inds1=[]
for i in range(len(labels3)):
    eigu[i,:]=eigu[i,:]*3.5
    ax.annotate('', xytext=(eigu[i,0], eigu[i,1]), xy=(0, 0),
                        arrowprops=dict(arrowstyle='<-', connectionstyle="arc3", color='k',lw=1),c='k')
    if np.linalg.norm(eigu[i,:])>0.9:
        inds1.append(i)
        if i==0:
            ax.annotate(labels3[i], xy=(eigu[i, 0], eigu[i, 1]), xytext=(eigu[i, 0] * 1.2-0.5,
                                                                         eigu[i, 1] * 1.2),
                        c='k', bbox=bbox, ha='center', va='center')
        elif i==1:
            ax.annotate(labels3[i], xy=(eigu[i, 0], eigu[i, 1]), xytext=(eigu[i, 0] * 1.2-1,
                                                                         eigu[i, 1] * 1.2),
                        c='k', bbox=bbox, ha='center', va='center')
        elif i==6:
            ax.annotate(labels3[i], xy=(eigu[i, 0], eigu[i, 1]), xytext=(eigu[i, 0] * 1.2-1, eigu[i, 1] * 1.2),
                        c='k', bbox=bbox, ha='center',va='center')
        else:
            ax.annotate(labels3[i], xy=(eigu[i, 0], eigu[i, 1]), xytext=(eigu[i, 0] * 1.2, eigu[i, 1] * 1.2),
                        c='k', bbox=bbox, ha='center',va='center')
print(inds1)
xlim=list(ax.get_xlim())
xlim[1]=7.5
xlim=[-8,4]
ylim=ax.get_ylim()
ax.hlines(0,xlim[0],xlim[1],lw=1,zorder=-10,colors='k')
ax.vlines(0,ylim[0],ylim[1],lw=1,zorder=-10,colors='k')
ax.set_xlim(xlim)
ax.set_ylim(ylim)
ax.legend(loc='lower right',bbox_to_anchor=(1.05,-0.03),frameon=False,title=None,handlelength=0.8)#
g.ax_marg_x.get_legend().remove()
g.ax_marg_y.get_legend().remove()
g.ax_joint.set_xlabel('CPCA$_1$')
g.ax_joint.set_ylabel('CPCA$_2$')
fig=plt.gcf()
fig.set_size_inches(5,5)
fig.subplots_adjust(bottom=0.1,left=0.1)
fig.savefig('Clinical_BM.svg',format='svg')

fig2,axes2=plt.subplots(nrows=2,ncols=4,figsize=(6,4))
axes2=axes2.flat
axes2[-1].set_visible(False)
for k in range(len(inds1)):
    sns.boxplot(x='Type',y=labels[inds1[k]],data=df3,palette=sns.set_palette(ggcolors), width=0.6, linewidth=1.0,
                ax=axes2[k], flierprops=flierprops)
    yy=max(df3[labels[inds1[k]]])
    axes2[k].set_ylabel('')
    axes2[k].set_xticks([])
    axes2[k].tick_params(labelsize=11)
    axes2[k].set_xlabel(labels[inds1[k]],fontsize=11)
    sFormatter1=matplotlib.ticker.ScalarFormatter(useOffset=False, useMathText=True)
    sFormatter1.set_powerlimits((0, 2))
    axes2[k].yaxis.set_major_formatter(sFormatter1)
    ylim=axes2[k].get_ylim()
    axes2[k].set_ylim([ylim[0],1.15*(ylim[1]-ylim[0])+ylim[0]])
    axes2[k].yaxis.offsetText.set_fontsize(9)
    axes2[k].yaxis.offsetText.set_position((-0.12, 1))
fig2.subplots_adjust(top=0.85,bottom=0.1,left=0.1,right=0.9,wspace=0.4,hspace=0.4)
fig2.text(x=0.6, y=0.88, s='Cell: 10$^6$/mL\t Cytokine: pg/mL', fontsize=10)
fig2.savefig('Clinical_BMStat.svg')
fig2.savefig('Clinical_BMStat.png')
plt.show()