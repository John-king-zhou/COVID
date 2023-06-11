#generate the aligned while averaged curves for mild/moderate, severe and critical patients
#figure S8
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

ggcolors=['#1F77B4','#FF7F0E','#D62728']
import matplotlib.mathtext as mathtext
mathtext.FontConstantsBase.sup1 = 0.2
mathtext.FontConstantsBase.sub1 = 0.2
mathtext.FontConstantsBase.sub2 = 0.3

Types=['Mild/Moderate','Severe','Critical']

bt_name = ['WBC', 'Neutrophil', 'Neutrophil%', 'Lymphocyte', 'Lymphocyte%', 'Monocyte', 'Monocyte%', 'IL-2', 'IL-4',
           'IL-6', 'IL-10', r'TNF-$\alpha$', r'IFN-$\gamma$','$\epsilon$*']
labels=['$\epsilon*$=(M%+N%)(M%+L%)','$\epsilon^{\#}$=N%'+r'$\times$'+'L%',
        '$E*$=(M+N)(M+L)']

CStat=['Mild/Moderate','Severe','Critical','Unlabeled']

unit = ['(10$^9$/mL)','10$^9$/mL','10$^9$/mL','10$^9$/mL','','','','pg/mL','pg/mL','pg/mL','pg/mL','pg/mL','pg/mL','pg/mL','']
w=1
excl = pd.ExcelFile('Data.xlsx')
fig,axes=plt.subplots(nrows=1,ncols=3,figsize=(6,2.5))
ax=axes.flat
selected=[1,3,5,2,4,6,9,10,13,14,15]
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
    indivs[j].append((indivs[j][2])*(indivs[j][4])/1e4)
    indivs[j].append((indivs[j][1]+indivs[j][5])*(indivs[j][3]+indivs[j][5]))
zorders=[-10,10,-10]

multiple=[[] for j in range(3)]
for j in [0,1,2]:
    data1=indivs[j][13]
    data2=indivs[j][9]
    for k in range(data1.shape[1]):
        if np.sum(1-np.isnan(data1[:,k]))>2 and np.sum(1-np.isnan(data2[:,k]))>2:
            multiple[j].append(k)
print(len(multiple[0]),len(multiple[1]),len(multiple[2]))

for j in [0,1,2]:
    for i in range(3):
        index=selected[i+8]
        data0=indivs[j][index]
        data=data0[0:25,np.array(multiple[j])]
        mean=np.nanmean(data,axis=1)
        indices=np.isnan(mean)
        time=np.arange(0,25,1)
        t=time[~indices]
        y=mean[~indices]
        if i!=2:
            y[y>0.3]=0.3
        ax[i].plot(t,y,c=ggcolors[j],zorder=zorders[j])
        ax[i].scatter(t,y,c=ggcolors[j], label=Types[j], s=16, zorder=zorders[j])
        ax[i].set_title(labels[i],fontsize=12)
        ax[i].tick_params(labelsize=12)
        ax[i].set_xticks([0,10,20])
        if i!=2:
            ax[i].set_ylim(0, 0.32)
            ax[i].set_yticks([0, 0.1, 0.2, 0.3])
            ax[i].set_yticklabels([0, 0.1, 0.2, '$\geq 0.3$'])

handles, labels = ax[0].get_legend_handles_labels()
fig.text(x=0.5,y=0.04,s='Days post admission',ha='center',va='center')
fig.legend(handles,labels,loc='lower right',bbox_to_anchor=(0.92,0.8),ncol=3,
           fontsize=12,markerscale=1.5,framealpha=0,handlelength=0.6,columnspacing=0.5)
fig.subplots_adjust(hspace=0.6,wspace=0.3,left=0.1,right=0.95,top=0.75,bottom=0.2)
fig.savefig('e_star.svg')
fig.savefig('e_star.png')


dfs=[]
for i in range(3):
    e=indivs[i][13]
    e2=indivs[i][14]
    E=indivs[i][15]
    e=e[0:25,np.array(multiple[i])]
    e=np.nanmean(e,axis=0)
    e2=e2[0:25,np.array(multiple[i])]
    e2=np.nanmean(e2,axis=0)
    E=E[0:25,np.array(multiple[i])]
    E=np.nanmean(E,axis=0)
    print(i,'average e',np.mean(e))
    print(i,'average e#',np.mean(e2))
    print(i,'average E',np.mean(E))
    e=e[~np.isnan(e)]
    e2=e2[~np.isnan(e2)]
    E=E[~np.isnan(E)]
    df=pd.DataFrame({'e1':e,'e2':e2,'E':E,'type':[CStat[i] for j in range(len(e))]})
    dfs.append(df)
dfs=pd.concat(dfs)
print(dfs.shape)
fig2,axes2=plt.subplots(3,1,figsize=(6,4))
xlabel=['e1','e2','E']
for i in range(3):
    sns.kdeplot(data=dfs,x=xlabel[i],hue='type',common_norm=False,palette=ggcolors,ax=axes2[i],legend=False)
ax1,ax2,ax3=axes2
ax1.set_xticks([0,0.2,0.4,0.6])
ax1.set_xlim(0,0.6)
ax1.set_yticks([])
ax1.set_xlabel('$\epsilon*$',labelpad=0,fontsize=14)
ax2.set_xticks([0,0.1,0.2,0.3,0.4])
ax2.set_xlim(0,0.4)
ax2.set_yticks([])
ax2.set_xlabel('$\epsilon^{\#}}$',labelpad=0,fontsize=14)
ax3.set_xticks([0,10,20,30,40])
ax3.set_xlim(0,40)
ax3.set_yticks([])
ax3.set_xlabel('$E*$',labelpad=0)
fig2.subplots_adjust(hspace=0.45,top=0.95)
fig2.savefig('e_star_distribution.svg')
plt.show()

