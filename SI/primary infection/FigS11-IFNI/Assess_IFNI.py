#assess the impact of IFN-I
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import multiprocessing
import seaborn as sns
from scipy.integrate import odeint
from Equation import func_IFNI
from E_Calculation_IFNI import E, Gamma
from mpl_toolkits import mplot3d
import matplotlib.patches as mpatches
import copy

ggcolors=['#808080','#1F77B4','#FF7F0E','#D62728','#4DBEEE','#77AC30','#9467BD']
import matplotlib.mathtext as mathtext
mathtext.FontConstantsBase.sup1 = 0.4

def Hill(x,k,n):
    return x**n/(k**n+x**n)

def get_results(Para0,IFNI):
    alpha_IFNI=1
    h_APC_IFNI=2
    h_NK_IFNI=2
    K_IFNI=14532
    Para=np.hstack((Para0, alpha_IFNI, h_APC_IFNI, h_NK_IFNI, K_IFNI, IFNI))
    dt=0.1
    Time=np.arange(0,35,dt)

    initial = [0.01, 0, Para[52] / Para[62], 0, 0, Para[53] / Para[65], 0, 0, Para[159], 0, 0, 0, 0, 0, 0,
               Para[54] / Para[70], 0, Para[160], 0, 0, 0, 0, 0, 0, Para[77] / Para[101], Para[82] / Para[102],
               Para[84] / Para[103], (Para[88] + Para[54] / Para[70] * Para[90]) / Para[104],
               Para[91] / Para[105], Para[95] / Para[106], 0, 0]
    results=odeint(func_IFNI, initial, Time, args=(Para,))
    e=E(results,Para)
    gamma=Gamma(results,Para)
    results=np.vstack((results.T,gamma,e,)).T
    return results

def get_results_merge(args):
    return get_results(*args)

def score(vf,hm,IL6m):
    IL6c=2000
    hc=30
    vc=1
    S=(1+1*Hill(hm,hc,1))*(1+2*Hill(IL6c,IL6m,1))*(1+1*Hill(vc,vf,1))
    return S

ls=['-','--','-.']
labels=["log$_{10}$nCoV","IL-6","$\gamma [H]/[H]_0$","$\epsilon$"]
bounds=[[-2,5],[0,3000],[0,4],[0,15],]
time = np.arange(0, 35, 0.1)
ggcolors=['#808080','#2CA02C','#1F77B4','#FF7F0E','#D62728','#4DBEEE','#77AC30','#9467BD']
sFormatter1 = matplotlib.ticker.ScalarFormatter(useOffset=False, useMathText=True)
sFormatter1.set_powerlimits((-1, 2))

if __name__=='__main__':
    cores=multiprocessing.cpu_count()
    pool=multiprocessing.Pool(processes=int(cores)-2)
    IFNI_list=[0,300,1000]
    fig,axes=plt.subplots(4,4,figsize=(8,6))
    for mode in range(1,5,1):
        print("Mode = %i"%mode)
        ax=axes[mode-1,:]
        Paras=np.loadtxt('Spl_Para%i.txt'%mode)#[0:5,:]
        for i in range(3):
            IFN_I=IFNI_list[i]
            args=[[Paras[j],IFN_I] for j in range(Paras.shape[0])]
            LOGV=[]
            IL6=[]
            GAMMA=[]
            EP=[]
            for y in pool.imap(get_results_merge, args):
                v=y[:,0]
                v[v<1e-6]=1e-6
                LOGV.append(np.log10(v))
                IL6.append(y[:,26])
                GAMMA.append(y[:,-2])
                EP.append(y[:,-1])
            X=[LOGV,IL6,GAMMA,EP]
            for j in range(4):
                ax[j].plot(time,np.mean(X[j],axis=0),c=ggcolors[mode],ls=ls[i],label="%i"%IFN_I)
                ax[j].set_xlim([0,35])
                ax[j].set_xticks(np.arange(0,36,7))
                if mode!=4:
                    ax[j].set_xticklabels([])
                ax[j].set_ylim(bounds[j])
                ax[j].grid()
        ax[0].set_ylabel("Mode %i"%mode)
    pool.close()
    for i in range(4):
        axes[0,i].set_title(labels[i])
        axes[3,i].set_xlabel("time (days)")
        axes[i,1].yaxis.set_major_formatter(sFormatter1)
    axes[0,0].legend(fontsize=8,title="IFN-I (pg/mL)",title_fontsize=8)
    fig.tight_layout()
    fig.savefig("Assess_IFNI.png",dpi=300)
    plt.show()
