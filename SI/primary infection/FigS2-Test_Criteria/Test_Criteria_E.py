import numpy as np
from E_Calculation import E,R0
from scipy.integrate import odeint
from Equation import func
import warnings
import matplotlib.pyplot as plt
import matplotlib.mathtext as mathtext

warnings.filterwarnings('error')
mathtext.FontConstantsBase.sup1 = 0.5
mathtext.FontConstantsBase.sub1 = 0.2
mathtext.FontConstantsBase.sub2 = 0.3
ggcolors=['#2CA02C','#1F77B4','#FF7F0E','#D62728',]

labels=['Mode 1','2','3','4',]
yticks=[[0,5,10],[0,2,4],[0,4,8]]

fig,ax1=plt.subplots(1,2,figsize=(8,3.7))
left, bottom, width, height = [0.3, 0.62, 0.15, 0.25]
ax21 = fig.add_axes([left, bottom, width, height])
left, bottom, width, height = [0.6, 0.62, 0.15, 0.25]
ax22 = fig.add_axes([left, bottom, width, height])
ax2=[ax21,ax22]
titles=["a=800, b=1600", "a=1200, b=2400"]
dt=0.1
time=np.arange(0,28.1,dt)
for ii in range(2):
    for i in range(1, 5, 1):
        data = np.loadtxt('./Test_Criteria/Spl_Para%i_%i.txt' % (ii,i))#[0:50,:]
        e_list = []
        R_list = []
        s = len(data)
        print('mode %i' % (i))
        for j in range(s):
            ratio = j / s
            rat_str = ['>'] * int(ratio * 50) + ['-'] * (50 - int(ratio * 50))
            rat_str = ''.join(rat_str)
            print('\r' + rat_str + '%.2f %%' % (ratio * 100), end='')
            Para = data[j]
            initial = [0.01, 0, Para[52] / Para[62], 0, 0, Para[53] / Para[65], 0, 0, Para[159], 0, 0, 0, 0, 0, 0,
                       Para[54] / Para[70], 0, Para[160], 0, 0, 0, 0, 0, 0, Para[77] / Para[101],
                       Para[82] / Para[102],
                       Para[84] / Para[103], (Para[88] + Para[54] / Para[70] * Para[90]) / Para[104],
                       Para[91] / Para[105], Para[95] / Para[106], 0, 0]
            results = odeint(func, initial, time, args=(Para,))
            e = E(results, Para)
            Rt=R0(results,Para)
            e_list.append(e)
            R_list.append(Rt)
        e_list = np.array(e_list)
        mean_e=10 ** np.mean(np.log10(e_list), axis=0)  # geometric mean
        mean_e = np.array(mean_e)
        logR_list = np.log10(R_list)
        mean_logR = np.mean(logR_list, axis=0)


        ax1[ii].grid(True,zorder=-20)
        ax2[ii].grid(True,zorder=-20)
        ax1[ii].plot(time,mean_e,color=ggcolors[i-1],linewidth=4,zorder=20,label=labels[i-1])
        ax2[ii].plot(time,mean_logR,color=ggcolors[i-1],linewidth=2,zorder=20)
    xlim=[0,28]
    ax2[ii].hlines(0,xlim[0],xlim[1],colors='k',lw=1,zorder=10)
    ax2[ii].set_xticks([0,7,14,21,28])
    ax2[ii].set_xlim(xlim)
    ax2[ii].set_ylabel('log$_{10}$R$_t$',labelpad=0)
    ax1[ii].hlines(3.6,0,28,colors='k',lw=1.3,zorder=10)
    ax1[ii].set_ylabel(r'$\epsilon$', fontsize=20, labelpad=10, rotation=0, ha='center', va='center')
    ax1[ii].set_xlabel('time (days)', fontsize=13, labelpad=1)
    ax1[ii].set_xlim([0,28])
    ax1[ii].set_ylim([0,15])
    ax1[ii].set_yticks([0,6,12])
    ax1[ii].set_xticks([7,14,21,28])
    ax1[ii].tick_params(labelsize=13,length=6)
    ax1[ii].set_title(titles[ii])
fig.subplots_adjust(left=0.12,right=0.92,bottom=0.13,top=0.9)
fig.savefig('./Test_Criteria/E-t-test_criteria.png',dpi=300)
plt.show()