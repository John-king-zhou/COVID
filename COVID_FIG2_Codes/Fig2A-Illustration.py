#plotting figure 1A, call for AverPlot
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from AveragePlot import *
from Equation import func
from E_Calculation import *
import warnings

warnings.filterwarnings('error')

plt.rcParams["xtick.major.size"] = 5
plt.rcParams["ytick.major.size"] = 5

ggcolors=['#808080','#2CA02C','#1F77B4','#FF7F0E','#D62728','#4DBEEE','#77AC30','#9467BD']

Para1=np.loadtxt('Spl_Para1.txt')
Para2=np.loadtxt('Spl_Para2.txt')
Para3=np.loadtxt('Spl_Para3.txt')
Para4=np.loadtxt('Spl_Para4.txt')
Paras=[Para1,Para2,Para3,Para4]
dt=0.1
time=np.arange(0,50,dt)

fig,axes=plt.subplots(nrows=4,ncols=6,figsize=(9.5,6.2))

for i in [0,1,2,3]:
    print(i)
    Typei=Paras[i]
    Traj=[[] for j in range(38)]
    try:
        Aver_Results=np.loadtxt('Mean%i.txt'%(i+1))
        Std_Results=np.loadtxt('Std%i.txt'%(i+1))
    except:
        Aver_Results=[]
        Std_Results=[]
        for k in range(len(Typei)):
            ratio=k/len(Typei)
            rat_str=['>']*int(ratio*50)+['-']*(50-int(ratio*50))
            rat_str=''.join(rat_str)
            print('\r'+rat_str+'%.2f %%' %(ratio*100), end='')
            Para=Typei[k]
            initial = [0.01, 0, Para[52] / Para[62], 0, 0, Para[53] / Para[65], 0, 0, Para[159], 0, 0, 0, 0, 0, 0,
                       Para[54] / Para[70], 0, Para[160], 0, 0, 0, 0, 0, 0, Para[77] / Para[101], Para[82] / Para[102],
                       Para[84] / Para[103], (Para[88] + Para[54] / Para[70] * Para[90]) / Para[104],
                       Para[91] / Para[105], Para[95] / Para[106], 0, 0]
            results=odeint(func,initial,time,args=(Para,))
            v=results[:,0]
            for j in range(32):
                if j==0:
                    try:
                        Traj[j].append(np.log10(results[:,j]))
                    except:
                        v=results[:,j]
                        v[v<10**(-4)]=10**(-4)
                        Traj[j].append(np.log10(v))
                        continue
                else:
                    Traj[j].append(results[:,j])
            CD8T = np.sum(results[:, 17:21], axis=1)
            Traj[32].append(np.log10(CD8T))
            #Traj[32].append(CD8T)
            Ig = results[:, 30]
            Ig[Ig < 1] = 1
            Traj[33].append(np.log10(Ig))
            #Traj[33].append(Ig)
            CD4T = np.sum(results[:, 8:14], axis=1) + results[:, 16]
            Traj[34].append(np.log10(CD4T))
            #Traj[34].append(CD4T)
            B = np.sum(results[:, 21:23], axis=1)+Para[150]
            Traj[35].append(np.log10(B))
            Traj[36].append(f_APC_anv(results,Para)*results[:,4])
            Traj[37].append(f_APC_inf(results, Para) * results[:, 4])
        for j in range(len(Traj)):
            A=np.array(Traj[j])
            Aver_Results.append(np.mean(A,axis=0))
            Std_Results.append(np.std(A,axis=0))
        Aver_Results=np.array(Aver_Results)
        Std_Results=np.array(Std_Results)
        np.savetxt('Mean%i.txt'%(i+1),Aver_Results)
        np.savetxt('Std%i.txt'%(i+1),Std_Results)
    Type_Specific(axes[i,:],time,Aver_Results,Std_Results,color=ggcolors[i+1],n=i)
fig.subplots_adjust(bottom=0.1,top=0.9,left=0.07,right=0.97,wspace=0.3)
fig.text(0.5,0.02,s='time (weeks)',fontsize=16,horizontalalignment='center')
fig.savefig('Illustration.svg')
fig.savefig('Illustration.png',dpi=300)
plt.show()