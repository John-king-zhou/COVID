#averPlot: function for plotting the assemble mean for all 24 variables and R0, e1, e2, e12 and 1/s,
#to be called in Mean_Comparison.py (generate figure S3)
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
# mode=int mode
def averPlot(axes,time,Mean,color,mode):
    Vars=['log$_{10}$(nCoV)', 'If', 'H', 'D', 'APC$^l$', 'APC$^u$', 'NK', 'Neut', 'CD4+T$_N$', 'CD4+T$_A$',
          'Th1', 'Th2', 'Th17', 'Tfh', 'Treg$^a$', 'Treg$^r$', 'CD4+T$_M$', 'CD8+T$_N$', 'CD8+T$_A$',
          'CTL', 'CD8+T$_M$', 'GC B', 'Plasma B', 'B$_M$' ,'IL-2', 'IL-4', 'IL-6', r'IL-10/TGF-$\beta$',
          r'TNF-$\alpha$', r'IFN-$\gamma$', 'Ab', 'A','R', r'$\epsilon$', r'$\epsilon_i$', r'$\epsilon_a$',
          r'$\epsilon_k$', r'$\epsilon_c$', 'CD4+T','CD8+T']
    axes[0].plot(time,np.zeros(len(time)),c='k',lw=1)
    for i in range(Mean.shape[1]):
        ax=axes[i]
        ax.tick_params(labelsize=11)
        sFormatter1=matplotlib.ticker.ScalarFormatter(useOffset=False, useMathText=True)
        sFormatter1.set_powerlimits((-1, 2))
        ax.yaxis.set_major_formatter(sFormatter1)
        ax.yaxis.offsetText.set_fontsize(8)
        ax.yaxis.offsetText.set_position((-0.12,1))
        x=Mean[:,i]
        ax.plot(time, x, linewidth=2, c=color, label='Mode %i'%(mode))
        ax.xaxis.set_major_locator(xmajorLocator)
        ax.set_title(Vars[i],fontsize=11,pad=3)
        ax.set_xticks([0,7,14,21,28,35])
        if i==0:
            ax.set_ylim([-3.2,5])
        if i in range(34,40):
            ax.set_xlabel('time(days)',fontsize=11,labelpad=-0.5)
    axes[40].set_visible(False)
    axes[41].set_visible(False)
    # axes[34].set_visible(False)
    # axes[35].set_visible(False)

#axes=plotting subplots; time=integration time point array;
# Mean=list containing the trajectory of the 30 variables; #color=str of rgb or name for a color;
# n=int mode-1
def Type_Specific(ax,time,Mean,Std,color,n):
    labels=['Mode 1','Mode 2','Mode 3','Mode 4']
    ulim=[52,5,3.5,0.9,3.2,3000]
    llim=[0,-2.1,-0.1,-2,0,-100]
    if n==0:
        ax[0].set_title('Epithelial Cell',fontsize=13)
        ax[1].set_title('log$_{10}$ Viral Load',fontsize=13)
        ax[2].set_title('Ag Presentation',fontsize=13)
        ax[3].set_title('log$_{10}$ CD8$^+$T',fontsize=13)
        ax[4].set_title('log$_{10}$ Abs',fontsize=13)
        ax[5].set_title('IL-6',fontsize=15)
    indices=[[2,],[0,],[36,],[32,],[33,],[26,]]
    for i in range(6):
        mean=np.sum(Mean[indices[i],:],axis=0)
        std=np.sum(Std[indices[i],:],axis=0)
        ax[i].plot(time, mean, linewidth=3, c=color)
        ax[i].fill_between(time, mean-std, mean+std, facecolor=color, alpha=0.2)
        ax[i].set_ylim(llim[i],ulim[i])
        yticks=np.arange(1,4,1)*(ulim[i]-llim[i])/4+llim[i]
        ax[i].set_yticks(yticks)
        ax[i].set_yticklabels([])
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
            ax[j].tick_params(labelsize=14)

def Type_Specific2(ax,time,Mean,Std,color,n):
    labels=['Mode 1','Mode 2','Mode 3','Mode 4']
    ulim=[52,5,2.8,4.5,1500,3300]
    llim=[0,-2.1,-0.1,0,0,-100]
    if n==0:
        ax[0].set_title('Epithelial Cell',fontsize=13)
        ax[1].set_title('log$_{10}$ Viral Load',fontsize=13)
        ax[2].set_title('Activated APC',fontsize=13)
        ax[3].set_title('CD8$^+$T',fontsize=13)
        ax[4].set_title('Abs',fontsize=13)
        ax[5].set_title('IL-6',fontsize=15)
    indices=[[2,],[0,],[4,],[32,],[33,],[26,]]
    for i in range(6):
        mean=np.sum(Mean[indices[i],:],axis=0)
        std=np.sum(Std[indices[i],:],axis=0)
        ax[i].plot(time, mean, linewidth=3, c=color)
        ax[i].fill_between(time, mean-std, mean+std, facecolor=color, alpha=0.2)
        ax[i].set_ylim(llim[i],ulim[i])
        yticks=np.arange(1,4,1)*(ulim[i]-llim[i])/4+llim[i]
        ax[i].set_yticks(yticks)
        ax[i].set_yticklabels([])
        ax[i].set_xlim([0, 28])
        ax[i].set_xticks([0, 7, 14, 21, 28])
        ax[i].set_xticklabels([])
    ax[0].set_ylabel(labels[n], fontsize=13)

    if n==3:
        for j in range(6):
            ax[j].set_xticks([0,7,14,21,28])
            ax[j].set_xticklabels([0,1,2,3,4])
            ax[j].tick_params(labelsize=14)
