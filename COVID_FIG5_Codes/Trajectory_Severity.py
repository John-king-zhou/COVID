#generate averaged curves for mild/moderate, severe and critical patients
#figure 5
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ggcolors=['#1F77B4','#FF7F0E','#D62728','#4DBEEE','#77AC30','#9467BD']
import matplotlib.mathtext as mathtext
mathtext.FontConstantsBase.sup1 = 0.5
mathtext.FontConstantsBase.sub1 = 0.2
mathtext.FontConstantsBase.sub2 = 0.3

Types=['Mild/Moderate','Severe','Critical']

bt_name = ['WBC', 'Neutrophil', 'Neutrophil%', 'Lymphocyte', 'Lymphocyte%', 'Monocyte', 'Monocyte%', 'IL-2', 'IL-4',
           'IL-6', 'IL-10', r'TNF-$\alpha$', r'IFN-$\gamma$','$\epsilon$*']
correspondence=[0,1,2,3,7,8,9,10,11,12]
unit = ['(10$^9$/mL)','10$^9$/mL','10$^9$/mL','10$^9$/mL','','','','pg/mL','pg/mL','pg/mL','pg/mL','pg/mL','pg/mL','pg/mL','']
w=1
excl = pd.ExcelFile('Data.xlsx')
fig,axes=plt.subplots(nrows=3,ncols=3,figsize=(6,4))
ax=axes.flat
selected=[1,3,5,2,4,6,9,10,13]
indices=[np.arange(1,91,1),np.arange(92,191,1),np.hstack((np.arange(192,204,1),np.arange(205,217,1)))]
indivs=[[] for j in range(3)]

for i in range(len(excl.sheet_names)):
    df=pd.read_excel(excl, sheet_name=excl.sheet_names[i],header=0,index_col=None)
    for j in range(3):
        #df1=df.iloc[0:20,indices[j]]
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
for j in [0,1,2]:
    for i in range(len(selected)):
        data0=indivs[j][selected[i]]
        data=data0[0:25,np.array(multiple[j])]
        mean=np.nanmean(data,axis=1)
        indices=np.isnan(mean)
        time=np.arange(0,25,1)
        t=time[~indices]
        y=mean[~indices]
        if i==2:
            y[y>1]=1
        if i==6:
            y[y>1000]=1000
        if i==8:
            y[y>0.3]=0.3
        ax[i].plot(t,y,c=ggcolors[j],zorder=zorders[j])
        ax[i].scatter(t,y,c=ggcolors[j], label=Types[j], s=16, zorder=zorders[j])
        ax[i].set_title(bt_name[selected[i]],fontsize=12)
        ax[i].tick_params(labelsize=10)
        ax[i].set_xticks([0,10,20])
fig.subplots_adjust(wspace=0.3,hspace=0.4)

ax[0].set_yticks([0,10,20])
ax[1].set_yticks([0,1,2])
ax[2].set_yticks([0.5,1])
ax[2].set_yticklabels([0.5,'$\geq 1$'])
ax[3].set_yticks([50,70,90])
ax[4].set_yticks([0,20,40])
ax[5].set_yticks([0,5,10,15])
ax[6].set_yticks([0,500,1000])
ax[6].set_yticklabels([0,500,'$\geq 1000$'])
ax[8].set_ylim(0,0.32)
ax[8].set_yticks([0, 0.1, 0.2, 0.3])
ax[8].set_yticklabels([0, 0.1, 0.2, '$\geq 0.3$'])
handles, labels = ax[0].get_legend_handles_labels()
fig.text(x=0.9, y=0.92, ha='right',va='bottom', s='Cell: 10$^6$/mL\tCytokine: pg/mL', fontsize=12)
fig.text(x=0.5,y=0.04,s='Days post admission',ha='center',va='center')
fig.legend(handles,labels,loc='lower right',bbox_to_anchor=(0.92,0.84),ncol=3,
           fontsize=12,markerscale=1.5,framealpha=0,handlelength=0.6,columnspacing=0.5)
fig.subplots_adjust(hspace=0.6,wspace=0.3,top=0.82)
fig.savefig('Clinical W%i_new_multiple.svg'%w)
fig.savefig('Clinical W%i_new_multiple.png'%w)
plt.show()

