#Type_Specific: function for plotting the illustration of mode 1, 2 and 3 with mean+std (figure 2A)
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

xmajorLocator = MultipleLocator(10)

#axes=plotting subplots; time=integration time point array;
# Mean=list containing the trajectory of the 30 variables; #color=str of rgb or name for a color;
# n=int mode-1
def Type_Specific(ax,time,Mean,Std,color,n):
    labels=['Mode 1','Mode 2','Mode 3','Mode 4']
    ulim=[52,5,3.5,0.9,3.2,3000]
    llim=[0,-2.1,-0.1,-2,0,-100]
    ticks=[[0,25,50],[-2,0,2,4],[0,1,2,3],[-2,-1,0,1],[0,1,2,3],[0,1000,2000,3000]]
    if n==0:
        ax[0].set_title('Epithelial Cell',fontsize=13)
        ax[1].set_title('log$_{10}$ Viral Load',fontsize=13)
        ax[2].set_title('Ag Presentation',fontsize=13)
        ax[3].set_title('log$_{10}$ CD8$^+$T',fontsize=13)
        ax[4].set_title('log$_{10}$ Abs',fontsize=13)
        ax[5].set_title('IL-6',fontsize=15)
    indices=[[2,],[0,],[36,],[32,],[33,],[26,]]
    sFormatter1 = matplotlib.ticker.ScalarFormatter(useOffset=False, useMathText=True)
    sFormatter1.set_powerlimits((-1, 2))
    ax[5].yaxis.set_major_formatter(sFormatter1)
    for i in range(6):
        mean=np.sum(Mean[indices[i],:],axis=0)
        std=np.sum(Std[indices[i],:],axis=0)
        ax[i].plot(time, mean, linewidth=3, c=color)
        ax[i].fill_between(time, mean-std, mean+std, facecolor=color, alpha=0.2)
        ax[i].set_ylim(llim[i],ulim[i])
        yticks=ticks[i]
        ax[i].set_yticks(yticks)
        ax[i].set_xlim([0, 28])
        ax[i].set_xticks([0, 7, 14, 21, 28])
        ax[i].set_xticklabels([])
        ax[i].tick_params(top=False, bottom=True, left=False, right=False)
        ax[i].grid(True, zorder=-20, lw=0.3)

    ax[0].set_ylabel(labels[n], fontsize=13)

    if n==3:
        for j in range(6):
            ax[j].set_xticks([0,7,14,21,28])
            ax[j].set_xticklabels([0,1,2,3,4])