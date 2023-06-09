#scatter map of 95 patients IL-6 and e*
#figure 5B
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import seaborn as sns
import matplotlib.mathtext as mathtext
from scipy.spatial import ConvexHull

mathtext.FontConstantsBase.sup1 = 0.5
mathtext.FontConstantsBase.sub1 = 0.2
mathtext.FontConstantsBase.sub2 = 0.3

ggcolors=['#1F77B4','#FF7F0E','#D62728','#808080']
markers=['o','s','^']
flierprops = dict(marker='_', markerfacecolor='#77AC30', markersize=6,
                  linestyle='none', markeredgecolor='#77AC30')
boxprops=dict(linestyle='-', color='k', facecolor='None', linewidth=1)

counts=[20,27,134,15]
Vars = ['WBC', 'Ne', 'Ly', 'Mo', 'Ne%', 'Ly%', 'Mo%', 'IL-2', 'IL-4','IL-6','IL-10', 'TNF', 'IFN','e']
CStat=['Mild/Moderate','Severe','Critical','Unlabeled']

def encircle(x,y, ax, **kw):

    p = np.c_[x,y]

    hull = ConvexHull(p)

    poly = plt.Polygon(p[hull.vertices,:], **kw)

    ax.add_patch(poly)

count=0

excl = pd.ExcelFile('Data.xlsx')
print(excl.sheet_names)
selected=[1,3,5,2,4,6,9,10,13]
indices=[np.arange(1,91,1),np.arange(92,191,1),np.hstack((np.arange(192,204,1),np.arange(205,217,1)))]
indivs=[[] for j in range(3)]

for i in range(len(excl.sheet_names)):
    df=pd.read_excel(excl, sheet_name=excl.sheet_names[i],header=0,index_col=None)
    for j in range(3):
        df1 = df.iloc[0:25, indices[j]]
        df1=np.array(df1,dtype=float)
        indivs[j].append(df1)
        if i in [1,2,3]:
            percentage=df1/(indivs[j][0])*100
            percentage[percentage>100]=100
            indivs[j].append(percentage)
for j in range(3):
    indivs[j].append((indivs[j][2]+indivs[j][6])*(indivs[j][4]+indivs[j][6])/1e4)
zorders=[-10,10,-10]

multiple=[[] for j in range(3)]
for j in [0,1,2]:
    data1=indivs[j][13]
    data2=indivs[j][9]
    for k in range(data1.shape[1]):
        if np.sum(1-np.isnan(data1[:,k]))>2 and np.sum(1-np.isnan(data2[:,k]))>2:
            multiple[j].append(k)
print(len(multiple[0]),len(multiple[1]),len(multiple[2]))
dfs=[]
for i in range(3):
    IL6=indivs[i][9]
    IL6=IL6[0:25,np.array(multiple[i])]
    e=indivs[i][13]
    e=e[0:25,np.array(multiple[i])]
    e=np.nanmean(e,axis=0)
    print(i,'average e',np.mean(e))
    IL6=np.nanmax(IL6,axis=0)
    indices=(np.isnan(e)+np.isnan(IL6))
    print(indices)
    e=e[~indices]
    IL6=IL6[~indices]
    df=pd.DataFrame({'e':e,'IL6':IL6,'type':[CStat[i] for j in range(len(e))]})
    dfs.append(df)
dfs=pd.concat(dfs)
dfs['logIL6']=np.log10(dfs.IL6)

g = sns.JointGrid()
sns.scatterplot(data=dfs,x='e',y='logIL6',hue='type',palette=sns.set_palette(ggcolors),ax=g.ax_joint)
sns.kdeplot(data=dfs,x='e',hue='type',palette=sns.set_palette(ggcolors),ax=g.ax_marg_x,
             linewidth=1,common_norm=False)
sns.kdeplot(data=dfs,y='logIL6',hue='type',palette=sns.set_palette(ggcolors),ax=g.ax_marg_y,
             linewidth=1,common_norm=False)

ax=g.ax_joint
for i in range(3):
    data=dfs[dfs.type==CStat[i]]
    encircle(np.array(data.e),np.array(data.logIL6), ax, ec="k", fc=ggcolors[i], alpha=0.4, zorder=-10, linewidth=0)
ax.set_xlabel('$<\epsilon*>$',labelpad=0)
ax.set_ylabel('$log_{10}IL-6_{max}$(pg/mL)',labelpad=0)
ax.get_legend().remove()
g.ax_marg_x.get_legend().remove()
g.ax_marg_y.get_legend().remove()
g.ax_marg_y.set_ylim([0.2,4.2])
z1 = np.polyfit(dfs.e,dfs.logIL6, 1)
p1 = np.poly1d(z1)
x2=np.arange(0,0.52,0.01)
y2=z1[0]*x2+z1[1]
r=pearsonr(dfs.e,np.log10(dfs.IL6))
print(p1,'pearson=',r[0])
ax.plot(x2,y2,color='k',linewidth=2,zorder=10)
ax.text(x=0.4,y=1,s='p=%1.2f'%r[0])
ax.set_ylim([0.2,4.2])
ax.set_xlim([-0.05,0.55])
fig=plt.gcf()
fig.set_size_inches(3.5,3.5)
fig.subplots_adjust(bottom=0.15,top=0.95,left=0.15,right=0.95)
fig.savefig('Scattermap_IL6_e_multiple.svg')
plt.show()