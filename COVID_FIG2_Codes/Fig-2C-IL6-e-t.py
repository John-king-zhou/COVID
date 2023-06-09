import numpy as np
import matplotlib.pyplot as plt
import multiprocessing
from scipy.integrate import odeint
import matplotlib
from Equation import func as func
from E_Calculation import E
import warnings

warnings.filterwarnings('error')

ggcolors=['#808080','#2CA02C','#1F77B4','#FF7F0E','#D62728','#4DBEEE','#77AC30','#9467BD']
markers=['o','^','s','v']
import matplotlib.mathtext as mathtext
mathtext.FontConstantsBase.sup1 = 0.5
mathtext.FontConstantsBase.sub1 = 0.2
mathtext.FontConstantsBase.sub2 = 0.3

def get_result(Para):
    time0=np.arange(0, 50, 0.1)
    initial = [0.01, 0, Para[52] / Para[62], 0, 0, Para[53] / Para[65], 0, 0, Para[159], 0, 0, 0, 0, 0, 0,
               Para[54] / Para[70], 0, Para[160], 0, 0, 0, 0, 0, 0, Para[77] / Para[101], Para[82] / Para[102],
               Para[84] / Para[103], (Para[88] + Para[54] / Para[70] * Para[90]) / Para[104], Para[91] / Para[105],
               Para[95] / Para[106], 0, 0]
    results = odeint(func, initial, time0, args=(Para,))
    e = E(results, Para)
    IL6 = results[:,26]
    return e,IL6

if __name__=='__main__':
    fig, ax = plt.subplots(1,1,figsize=(5,2))
    cores = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=int(cores-2))
    alphas=[0.03,0.05,0.05,0.03]
    for mode in [1,2,3,4]:
        print('mode=',mode)
        try:
            IL6=np.loadtxt('IL6_%i.txt'%mode)
            e=np.loadtxt('E_%i.txt'%mode)
        except:
            Paras=np.loadtxt('Spl_Para%i.txt'%mode)
            IL6=[]
            e=[]
            N=Paras.shape[0]
            count=0
            for y in pool.imap(get_result,Paras):
                count+=1
                ratio=count/N
                rat_str=['>']*int(ratio*50)+['-']*(50-int(ratio*50))
                rat_str=''.join(rat_str)
                print('\r'+rat_str+'%.2f %%' %(ratio*100), end='')
                e.append(y[0])
                IL6.append(y[1])
            e=np.array(e)
            IL6=np.array(IL6)
            np.savetxt('IL6_%i.txt'%mode,IL6)
            np.savetxt('E_%i.txt'%mode,e)
        n=IL6.shape[0]
        for i in range(n):
            ax.plot(e[i,:],IL6[i,:],color=ggcolors[mode],alpha=alphas[mode-1],linewidth=1,zorder=-10)
        mean_e=10 ** np.mean(np.log10(e), axis=0)
        mean_IL6=np.mean(IL6,axis=0)
        ax.plot(mean_e,mean_IL6,color=ggcolors[mode],linewidth=2,zorder=-mode)
        T=len(mean_e)
        T=[80,140,170,90]
        t=T[mode-1]
        loc=np.array([mean_e[t],mean_IL6[t]])
        velocity=np.array([mean_e[t+1],mean_IL6[t+1]])-loc
        ax.quiver(loc[0],loc[1],velocity[0],velocity[1],color=ggcolors[mode],angles='xy',headwidth=5)
    sFormatter1=matplotlib.ticker.ScalarFormatter(useOffset=False, useMathText=True)
    sFormatter1.set_powerlimits((-1, 1))
    ax.yaxis.set_major_formatter(sFormatter1)
    ax.set_xlim([0,10])
    ax.set_xticks([5,10])
    ax.set_yticks([0,1000,2000,3000])
    ax.set_ylim([0,3000])
    ax.set_xlabel(r'$\epsilon$',fontsize=13,labelpad=-3)
    ax.set_ylabel('[IL-6] (pg/mL)')
    fig.subplots_adjust(top=0.9,bottom=0.18,left=0.18,right=0.9,wspace=0.4, hspace=0.6)
    # fig.savefig('IL6-e-t.png')
    fig.savefig('IL6-e-t.svg')
    plt.show()