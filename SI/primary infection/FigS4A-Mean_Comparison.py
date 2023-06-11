#plotting the ensemble-averaged curves of mode 1, 2, 3 and asymptomatic. (figure S3)
#call AveragePlot
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing
from scipy.integrate import odeint
from AveragePlot import averPlot
from Equation import func as func
from E_Calculation import E,E1,E2,E12,R0,E_kill,E_clear
import warnings

warnings.filterwarnings('error')

ggcolors=['#808080','#2CA02C','#1F77B4','#FF7F0E','#D62728','#4DBEEE','#77AC30','#9467BD']
import matplotlib.mathtext as mathtext
mathtext.FontConstantsBase.sup1 = 0.5
mathtext.FontConstantsBase.sub1 = 0.2
mathtext.FontConstantsBase.sub2 = 0.3

def get_result(Para):
    time0=np.arange(0, 35, 0.1)
    initial = [0.01, 0, Para[52] / Para[62], 0, 0, Para[53] / Para[65], 0, 0, Para[159], 0, 0, 0, 0, 0, 0,
               Para[54] / Para[70], 0, Para[160], 0, 0, 0, 0, 0, 0, Para[77] / Para[101], Para[82] / Para[102],
               Para[84] / Para[103], (Para[88] + Para[54] / Para[70] * Para[90]) / Para[104], Para[91] / Para[105],
               Para[95] / Para[106], 0, 0]
    results = odeint(func, initial, time0, args=(Para,))
    r0 = R0(results, Para)
    e = E(results, Para)
    e1 = E1(results, Para)
    e2 = E2(results, Para)
    e12 = E12(results, Para)
    e_k = E_kill(results, Para)
    e_c = E_clear(results, Para)
    CD4T = np.sum(results[:, 8:17], axis=1)
    CD8T = np.sum(results[:, 17:21], axis=1)
    CD4T[CD4T < 1e-4] = 1e-4
    CD8T[CD8T < 1e-4] = 1e-4
    results = np.vstack((results.T, np.log10(r0), e, e1, e2+e12, e_k, e_c, CD4T, CD8T)).T
    return results

if __name__=='__main__':
    time0 = np.arange(0, 80, 0.1)
    time = np.arange(0, 35, 0.1)
    fig, axes = plt.subplots(nrows=7, ncols=6, figsize=(10, 8))
    axes = axes.flat
    cores = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=int(cores-2))
    for mode in [1,2,3,4]:
        print('mode=',mode)
        Paras = np.loadtxt('Spl_Para%i.txt' % mode)
        print(Paras.shape[0])
        try:
            result=np.loadtxt('Spl_mean%i.txt'%mode)
        except:
            Paras=np.loadtxt('Spl_Para%i.txt'%mode)
            Traj=[[] for j in range(40)]
            Aver_Results=[]
            N=Paras.shape[0]
            count=0
            for results in pool.imap(get_result,Paras):
                count+=1
                ratio=count/N
                rat_str=['>']*int(ratio*50)+['-']*(50-int(ratio*50))
                rat_str=''.join(rat_str)
                print('\r'+rat_str+'%.2f %%' %(ratio*100), end='')
                for j in range(results.shape[1]):
                    if j==0:
                        v = results[:, j]
                        v[v < 1e-6] = 1e-6
                        Traj[j].append(np.log10(v))
                    else:
                        Traj[j].append(results[:, j])
            for j in range(results.shape[1]):
                X=np.array(Traj[j])
                Aver_Results.append(np.mean(X, axis=0))
            Aver_Results=np.array(Aver_Results).T
            np.savetxt('Spl_mean%i.txt'%mode, Aver_Results)
            result=Aver_Results
        print(result.shape,time.shape)
        averPlot(axes,time,result,ggcolors[mode],mode)
    handles,labels=axes[0].get_legend_handles_labels()
    fig.legend(handles,labels,fontsize=11,ncol=2,bbox_to_anchor=(0.91,1.005),loc='upper right',frameon=False,columnspacing=0.5)
    fig.text(x=0.1,y=0.94,s='Cell Density: 10$^6$/mL\t Ig:$\mu$g/mL\nCytokine: pg/mL',fontsize=11)
    fig.subplots_adjust(top=0.9,bottom=0.05,left=0.1,right=0.9,wspace=0.4, hspace=0.6)
    fig.savefig('Spl_mean_curve.svg')
    fig.savefig('Spl_mean_curve.png')
    plt.show()