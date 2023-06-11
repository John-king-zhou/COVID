# Comparison between mode 4.1 and 4.2 (figure S3)

import numpy as np
import matplotlib.pyplot as plt
import warnings
from AveragePlot import averPlot_Supply

warnings.filterwarnings('error')

ggcolors=['#808080','#1F77B4','#FF7F0E','#D62728','#2CA02C','#4DBEEE','#77AC30','#9467BD']
import matplotlib.mathtext as mathtext
# mathtext.FontConstantsBase.sup1 = 0.5
# mathtext.FontConstantsBase.sub1 = 0.2
# mathtext.FontConstantsBase.sub2 = 0.3

# loda data, simulated in '_S_classification_mode4.py'
time = np.arange(0, 50, 0.1)
Traj, Aver_Results = np.load('Data_M4_S_comparison.npy', allow_pickle=True)


## -----------------------------------------------------------------------------------
## classifications by S_CD8 Changes (S_CD8(t=end)-SCD8(0) > 0, enough)----------------

# --------- typical variables ---------------------------------------------------------
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(9, 2.5))
axes = axes.flat
Traj_sort=[[],[]]
count = [0,0]
for j in range(len(Traj[-2])):
    cd4_tmp = np.array(Traj[-2][j])
    cd8_tmp = np.array(Traj[-1][j])

    y_tmp = np.array([np.log10(cd8_tmp),np.log10(cd4_tmp)])
    if cd8_tmp[-1]-cd8_tmp[0]>0:
        color_tmp = 'g'
        Traj_sort[0].append(y_tmp)
    else:
        color_tmp = 'r'
        Traj_sort[1].append(y_tmp)
        
    for i in range(len(axes)):
        axes[i].plot(time, y_tmp[i], alpha=0.1, lw=0.8, c=color_tmp)
        
    
# mean traj by classification
Traj_sort[0] = np.array(Traj_sort[0])
Traj_sort[1] = np.array(Traj_sort[1])
for j in range(Traj_sort[0].shape[1]):
    axes[j].plot(time, np.mean(Traj_sort[0][:,j], axis=0),
                 'g--', label=r'$\Delta J_{S,CD8} > 0$')
    axes[j].plot(time, np.mean(Traj_sort[1][:,j], axis=0),
                 'r--', label=r'$\Delta J_{S,CD8} \leq 0$')

# mean traj of all
# axes[0].plot(time, np.log10(Aver_Results[:,-2]), alpha=1, c='k', ls='--')
# axes[1].plot(time, np.log10(Aver_Results[:,-1]), alpha=1, c='k', ls='--')

axes[0].legend(loc='best', ncol=1, columnspacing=2,
               )

axes[0].set_ylabel(r'$log_{10}\left(J_{S,CD8}\right)$')
axes[1].set_ylabel(r'$log_{10}\left(J_{S,CD4}\right)$')

axes[0].set_xlabel('Days')
axes[1].set_xlabel('Days')
plt.subplots_adjust(top=0.8,
                    bottom=0.2,
                    left=0.15,
                    right=0.95,
                    hspace=0.2,
                    wspace=0.5)

# fig.savefig('S_CD8_CD4_Illustrate.png',dpi=200)
fig.savefig('S1B_S_CD8_CD4_Illustrate.svg',dpi=200)


# ---------------- All variables -------------------------------------------------
Traj_good=[[] for j in range(40)]
Traj_bad=[[] for j in range(40)]
for j in range(len(Traj[0])):
    cd8_tmp = np.array(Traj[-1][j])
    if cd8_tmp[-1]-cd8_tmp[0]>0:
        for p in range(len(Traj)):
            # log S_cd8 and S_cd4
            if p in [38,39]:
                Traj_good[p].append(np.log10(Traj[p][j]))
            else:
                Traj_good[p].append(Traj[p][j])
    else:
        for p in range(len(Traj)):
            if p in [38,39]:
                Traj_bad[p].append(np.log10(Traj[p][j]))
            else:
                Traj_bad[p].append(Traj[p][j])

Aver_good=[]
Aver_bad=[]
Std_good=[]
Std_bad=[]
for j in range(len(Traj)):
    X=np.array(Traj_good[j])
    Aver_good.append(np.mean(X, axis=0))
    Std_good.append(np.std(X, axis=0))

    X=np.array(Traj_bad[j])
    Aver_bad.append(np.mean(X, axis=0))
    Std_bad.append(np.std(X, axis=0))
Aver_good=np.array(Aver_good).T
Aver_bad=np.array(Aver_bad).T
Std_good=np.array(Std_good).T
Std_bad=np.array(Std_bad).T
# plot
fig, axes = plt.subplots(nrows=7, ncols=6, figsize=(10, 8))
axes = axes.flat

averPlot_Supply(axes,time,Aver_good,'g',4.1,Std_good)
averPlot_Supply(axes,time,Aver_bad,'r',4.2,Std_good)

handles,labels=axes[0].get_legend_handles_labels()
fig.legend(handles,labels,fontsize=11,ncol=2,bbox_to_anchor=(0.91,1.005),loc='upper right',frameon=False,columnspacing=0.5)
fig.text(x=0.1,y=0.94,s='Cell Density: 10$^6$/mL\t Ig:$\mu$g/mL\nCytokine: pg/mL',fontsize=11)
fig.subplots_adjust(top=0.9,bottom=0.05,left=0.1,right=0.9,wspace=0.4, hspace=0.6)

# fig.savefig('Mode4-1-and-2.png',dpi=200)
fig.savefig('S1C_Mode4-1-and-2.svg',dpi=200)

plt.show()

