#plotting the ensemble-averaged curves of mode 1, 2, 3 and asymptomatic. (figure S3)
#call AveragePlot
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing
from scipy.integrate import odeint
from Equation import func as func
from E_Calculation import E,E1,E2,E12,R0,E_kill,E_clear,S_CD8,S_CD4
import warnings

warnings.filterwarnings('error')

ggcolors=['#808080','#1F77B4','#FF7F0E','#D62728','#2CA02C','#4DBEEE','#77AC30','#9467BD']
import matplotlib.mathtext as mathtext
# mathtext.FontConstantsBase.sup1 = 0.5
# mathtext.FontConstantsBase.sub1 = 0.2
# mathtext.FontConstantsBase.sub2 = 0.3

def get_result(Para):
    time0 = np.arange(0, 50, 0.1)
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
    s_cd4 = S_CD4(results, Para)
    s_cd8 = S_CD8(results, Para)
    results = np.vstack((results.T, r0, e, e1, e2+e12, e_k, e_c, s_cd4, s_cd8)).T
    return results

if __name__=='__main__':
    mode = 4
    print('mode=',mode)
    time = np.arange(0, 50, 0.1)
    try:
        Traj, Aver_Results = np.load('Data_M%i_S_comparison.npy' %mode, allow_pickle=True)
    except:
        cores = multiprocessing.cpu_count()
        pool = multiprocessing.Pool(processes=int(cores-2))

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
                elif j==32:
                    # R-->lg(R)
                    Traj[j].append(np.log10(results[:, j]))
                else:
                    Traj[j].append(results[:, j])
        for j in range(results.shape[1]):
            X=np.array(Traj[j])
            Aver_Results.append(np.mean(X, axis=0))
        Aver_Results=np.array(Aver_Results).T
        pool.close()

        np.save('Data_M%i_S_comparison.npy' %mode, np.array([Traj, Aver_Results], dtype=object))
            
    ## plot S CD4 and CD8----------------------------------------------------------------

    # fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10, 4))
    # axes = axes.flat
    # for j in range(len(Traj[-2])):
    #     cd4_tmp = Traj[-2][j]
    #     cd8_tmp = Traj[-1][j]
    #     axes[0].plot(time, np.log10(cd4_tmp), alpha=0.2, lw=0.5, c='b')
    #     axes[1].plot(time, np.log10(cd8_tmp), alpha=0.2, lw=0.5, c='g')
    # axes[0].plot(time, np.log10(Aver_Results[:,-2]), alpha=1, c='k', ls='--')
    # axes[1].plot(time, np.log10(Aver_Results[:,-1]), alpha=1, c='k', ls='--')
    # axes[0].set_title('lg(S_CD4)')
    # axes[1].set_title('lg(S_CD8)')
    # axes[0].set_xlabel('Days')
    # axes[1].set_xlabel('Days')
    # plt.tight_layout()

    ## -----------------------------------------------------------------------------------
    ## classifications by S_CD8 Changes (S_CD8(t=end)-SCD8(0) > 0, enough)----------------

    # --------- 4 typical variables ---------------------------------------------------------
    fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(10, 9))
    axes = axes.flat
    Traj_sort=[[],[]]
    for j in range(len(Traj[-2])):
        lgv = np.array(Traj[0][j])
        lgRt = np.array(Traj[32][j])
        cd4_tmp = np.array(Traj[-2][j])
        cd8_tmp = np.array(Traj[-1][j])
        ek_tmp = np.array(Traj[-4][j])
        ec_tmp = np.array(Traj[-3][j])

        y_tmp = np.array([lgv, lgRt, np.log10(cd4_tmp),np.log10(cd8_tmp),ek_tmp,ec_tmp])
        if cd8_tmp[-1]-cd8_tmp[0]>0:
            color_tmp = 'g'
            Traj_sort[0].append(y_tmp)
        else:
            color_tmp = 'r'
            Traj_sort[1].append(y_tmp)

        for i in range(len(axes)):
            axes[i].plot(time, y_tmp[i], alpha=0.2, lw=0.5, c=color_tmp)
        
    # mean traj by classification
    Traj_sort[0] = np.array(Traj_sort[0])
    Traj_sort[1] = np.array(Traj_sort[1])
    for j in range(Traj_sort[0].shape[1]):
        axes[j].plot(time, np.mean(Traj_sort[0][:,j], axis=0), 'g--')
        axes[j].plot(time, np.mean(Traj_sort[1][:,j], axis=0), 'r--')

    # mean traj of all
    axes[0].plot(time, Aver_Results[:,0], alpha=1, c='k', ls='--')
    axes[1].plot(time, Aver_Results[:,32], alpha=1, c='k', ls='--')
    axes[2].plot(time, np.log10(Aver_Results[:,-2]), alpha=1, c='k', ls='--')
    axes[3].plot(time, np.log10(Aver_Results[:,-1]), alpha=1, c='k', ls='--')
    axes[4].plot(time, Aver_Results[:,-4], alpha=1, c='k', ls='--')
    axes[5].plot(time, Aver_Results[:,-3], alpha=1, c='k', ls='--')

    axes[0].set_title('lg(nCoV)')
    axes[1].set_title('lg(Rt)')
    axes[2].set_title('lg(S_CD4)')
    axes[3].set_title('lg(S_CD8)')
    axes[4].set_title(r'$\epsilon_k$')
    axes[5].set_title(r'$\epsilon_c$')

    axes[4].set_xlabel('Days')
    axes[5].set_xlabel('Days')
    plt.tight_layout()

    # ---------------- All variables -------------------------------------------------
    # Traj_good=[[] for j in range(40)]
    # Traj_bad=[[] for j in range(40)]
    # for j in range(len(Traj[0])):
    #     cd8_tmp = np.array(Traj[-1][j])
    #     if cd8_tmp[-1]-cd8_tmp[0]>0:
    #         for p in range(len(Traj)):
    #             Traj_good[p].append(Traj[p][j])
    #     else:
    #         for p in range(len(Traj)):
    #             Traj_bad[p].append(Traj[p][j])

    # Aver_good=[]
    # Aver_bad=[]
    # Std_good=[]
    # Std_bad=[]
    # for j in range(len(Traj)):
    #     X=np.array(Traj_good[j])
    #     Aver_good.append(np.mean(X, axis=0))
    #     Std_good.append(np.std(X, axis=0))

    #     X=np.array(Traj_bad[j])
    #     Aver_bad.append(np.mean(X, axis=0))
    #     Std_bad.append(np.std(X, axis=0))
    # Aver_good=np.array(Aver_good).T
    # Aver_bad=np.array(Aver_bad).T
    # Std_good=np.array(Std_good).T
    # Std_bad=np.array(Std_bad).T
    # # plot
    # fig, axes = plt.subplots(nrows=7, ncols=6, figsize=(10, 8))
    # axes = axes.flat

    # averPlot(axes,time,Aver_good,'g',4.1,Std_good)
    # averPlot(axes,time,Aver_bad,'r',4.2,Std_good)

    # handles,labels=axes[0].get_legend_handles_labels()
    # fig.legend(handles,labels,fontsize=11,ncol=2,bbox_to_anchor=(0.91,1.005),loc='upper right',frameon=False,columnspacing=0.5)
    # fig.text(x=0.1,y=0.94,s='Cell Density: 10$^6$/mL\t Ig:$\mu$g/mL\nCytokine: pg/mL',fontsize=11)
    # fig.subplots_adjust(top=0.9,bottom=0.05,left=0.1,right=0.9,wspace=0.4, hspace=0.6)

    plt.show()
    
