#plot the ensemble-averaged trajectory in the space of v-IL6-e
#run in parallel for efficiency
#figure S12 D~F
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import multiprocessing
import seaborn as sns
from scipy.integrate import odeint
from Equation import func_med
from E_Calculation import E,E1,E2,E12
from mpl_toolkits import mplot3d
import matplotlib.patches as mpatches
import copy

ggcolors=['#808080','#1F77B4','#FF7F0E','#D62728','#4DBEEE','#77AC30','#9467BD']
import matplotlib.mathtext as mathtext
mathtext.FontConstantsBase.sup1 = 0.4

def Hill(x,k,n):
    return x**n/(k**n+x**n)

def get_results(Para0,period1,period2,period3,period4,mode=0):
    alpha=1
    epsilon_IFNI=1
    h_APC_IFNI=0
    h_NK_IFNI=0
    k_c_5=0.006
    Ig_ex=0
    beta_GC=1
    d_lymph_GC=0
    Para=np.hstack((Para0, alpha, epsilon_IFNI, h_APC_IFNI, h_NK_IFNI, k_c_5, Ig_ex, beta_GC, d_lymph_GC))
    funcc=func_med
    dt=0.1
    Time=np.arange(0,80,dt)

    initial = [0.01, 0, Para[52] / Para[62], 0, 0, Para[53] / Para[65], 0, 0, Para[159], 0, 0, 0, 0, 0, 0,
               Para[54] / Para[70], 0, Para[160], 0, 0, 0, 0, 0, 0, Para[77] / Para[101], Para[82] / Para[102],
               Para[84] / Para[103], (Para[88] + Para[54] / Para[70] * Para[90]) / Para[104],
               Para[91] / Para[105], Para[95] / Para[106], 0, 0]
    results0=odeint(funcc, initial, Time, args=(Para,))
    T=[0,7,14,80]
    times=[np.arange(T[i],T[i+1],dt) for i in range(3)]
    Traj=[]
    for i in range(3):
        if period1[0]<=T[i] and period1[1]>=T[i+1]:
            Para[161] = 0.8
        else:
            Para[161] = 1
        if period2[0]<=T[i] and period2[1]>=T[i+1]:
            Para[162] = 0.95
            Para[163] = 0.1
            Para[164] = 0.1
        else:
            Para[162] = 1
            Para[163] = 0
            Para[164] = 0
        if period3[0]<=T[i] and period3[1]>=T[i+1]:
            Para[166] = 400
        else:
            Para[166] = 0
        if period4[0]<=T[i] and period4[1]>=T[i+1]:
            Para[167] = 0.6
            Para[168] = 0.1
        else:
            Para[167]=1
            Para[168]=0
        if i==0:
            initial_i=initial*1
        else:
            initial_i=Traj[i-1][-1,:]
        results_i=odeint(funcc,initial_i,times[i],args=(Para,))
        Traj.append(results_i)
    results=np.vstack(Traj)
    e0=E(results0,Para)
    ei0=E1(results0,Para)
    ea0=e0-ei0
    e=E(results,Para,period2=period2,period3=period3)
    ei=E1(results,Para,period2=period2,period3=period3)
    ea=e-ei#E2(results,Para,period2=period2,period3=period3)+E12(results,Para,period2=period2,period3=period3)
    results=np.vstack((results.T,e,ei,ea)).T
    results0=np.vstack((results0.T,e0,ei0,ea0)).T
    return results,results0

def get_results_merge(args):
    return get_results(*args)

def score(vf,hm,IL6m):
    IL6c=2000
    hc=30
    vc=1
    S=(1+1*Hill(hm,hc,1))*(1+2*Hill(IL6c,IL6m,1))*(1+1*Hill(vc,vf,1))
    return S

dt = 0.1
if __name__=='__main__':
    line_labels=['No treatment','Treated']
    line_styles=['--','-']
    ds=[(5,0.6),(5,0)]
    fig,ax=plt.subplots(1,1)
    lines=[plt.plot([], [], linewidth=2, color='k', linestyle=line_styles[i], dashes=ds[i],label=line_labels[i])[0] for i in range(2)]
    scatter_labels=['p.i. 7','p.i. 14','healthy',r'[nCoV]$\neq$0']
    markers=['o','s','^','X']
    colors=['k','k','k','k']
    fc=['w','w','w','w']
    scatters=[plt.scatter([], [], marker=markers[i], s=60, edgecolor=colors[i], facecolor=fc[i], label=scatter_labels[i]) for i in range(4)]
    fig.legend(handles=lines,labels=line_labels,loc='lower left',bbox_to_anchor=(0.2,0.4),ncol=1,frameon=False,handlelength=1.5,columnspacing=0.3)
    fig.legend(handles=scatters,labels=scatter_labels,loc='lower left',bbox_to_anchor=(0.2,0.6),ncol=2,frameon=False,handlelength=1.5,columnspacing=0.3)
    fig.savefig('MedTraj_Legend.svg',format='svg')

    TP1=[(0,100),(0,0),(0,0),(0,0)]
    TP2=[(0,100),(0,7),(0,0),(7,14)]
    TP3=[(0,100),(0,7),(14,100),(7,14)]
    TreatmentPlans=[TP1,TP2,TP3]
    cores=multiprocessing.cpu_count()
    pool=multiprocessing.Pool(processes=int(cores)-4)
    try:
        Scrs=np.loadtxt('Scrs.txt')
    except:
        Scrs=[]
    for mode in range(2,5,1):
        TP=TreatmentPlans[mode-2]
        fig=plt.figure(figsize=(3.7, 2.7))
        ax=plt.axes(projection='3d', proj_type='ortho')
        fig2=plt.figure(figsize=(3.7, 2.7))
        ax2=plt.axes(projection='3d', proj_type='ortho')
        print(mode)
        print(TP)
        try:
            results0=np.loadtxt('Non_MedTraj%i.txt'%(mode))
            results=np.loadtxt('MedTraj%i.txt'%(mode))
        except:
            Paras=np.loadtxt('Spl_Para%i.txt'%mode)
            Sc=[]
            v0=[]
            e0=[]
            IL60=[]
            ei0=[]
            ea0=[]
            v=[]
            e=[]
            IL6=[]
            ei=[]
            ea=[]
            args=[(Paras[i],TP[0],TP[1],TP[2],TP[3],mode) for i in range(Paras.shape[0])]
            count=0
            for y in pool.imap(get_results_merge,args):
                results=y[0]
                results0=y[1]
                scr=score(results[-1,0],min(results[:,2]),max(results[:,26]))
                scr0=score(results0[-1,0],min(results0[:,2]), max(results0[:,26]))
                Sc.append((scr-scr0)/scr0)
                v0.append(results0[:,0])
                v.append(results[:,0])
                e0.append(results0[:,32])
                e.append(results[:,32])
                ei0.append(results0[:,33])
                ei.append(results[:,33])
                ea0.append(results0[:,34])
                ea.append(results[:,34])
                IL60.append(results0[:,26])
                IL6.append(results[:,26])
                count+=1
            v0=np.array(v0)
            e0=np.array(e0)
            ei0=np.array(ei0)
            ea0=np.array(ea0)
            IL60=np.array(IL60)
            v=np.array(v)
            e=np.array(e)
            ei=np.array(ei)
            ea=np.array(ea)
            IL6=np.array(IL6)
            Sc = np.array(Sc)
            Q1 = np.percentile(Sc, 25, interpolation='midpoint')
            Q3 = np.percentile(Sc, 75, interpolation='midpoint')
            IQR = Q3 - Q1
            indices = (Sc > Q3 + 1.5 * IQR) + (Sc < Q1 - 1.5 * IQR)
            print('mode=%i'%mode,'delta Q/Q=%1.4f'%np.mean(Sc))#[~indices]
            Scrs.append(np.mean(Sc))
            # v0=np.mean(v0[~indices,:],axis=0)
            # e0=np.mean(e0[~indices,:],axis=0)
            # IL60=np.mean(IL60[~indices,:],axis=0)
            # v=np.mean(v[~indices,:],axis=0)
            # e=np.mean(e[~indices,:],axis=0)
            # IL6=np.mean(IL6[~indices,:],axis=0)
            v0=np.mean(v0,axis=0)
            e0=np.mean(e0,axis=0)
            ei0=np.mean(ei0,axis=0)
            ea0=np.mean(ea0,axis=0)
            IL60=np.mean(IL60,axis=0)
            v=np.mean(v,axis=0)
            e=np.mean(e,axis=0)
            ei=np.mean(ei,axis=0)
            ea=np.mean(ea,axis=0)
            IL6=np.mean(IL6,axis=0)
            results0=np.array([v0,e0,IL60,ei0,ea0]).T
            results=np.array([v,e,IL6,ei,ea]).T
            np.savetxt('Non_MedTraj%i.txt'%(mode),results0)
            np.savetxt('MedTraj%i.txt'%(mode),results)
        v0 = results0[:, 0]
        e0 = results0[:, 1]
        IL60 = results0[:, 2]
        ei0 = results0[:, 3]
        ea0 = results0[:, 4]
        v = results[:, 0]
        e = results[:, 1]
        IL6 = results[:, 2]
        ei = results[:, 3]
        ea = results[:, 4]
        v0[v0>3000]=3000
        v[v>3000]=3000
        if mode==4:
            e[e>4]=4
            ea[ea>2]=2
        '''
        IL-6-e-v
        '''
        ax.plot(IL60,e0,np.zeros(len(v0)),linewidth=1,linestyle='-',c=ggcolors[mode-1],alpha=0.5)
        ax.plot(IL6,e,np.zeros(len(v)),linewidth=1,c=ggcolors[mode-1],alpha=0.5)
        ax.plot(IL60,e0,v0,linewidth=2,linestyle='-',dashes=(5,0.6),c=ggcolors[mode-1],label='No treatment')
        ax.scatter(IL60[int(7/dt)],e0[int(7/dt)],v0[int(7/dt)],s=30,edgecolor='k',facecolor='w',marker='o',zorder=-1)
        ax.scatter(IL60[int(14/dt)],e0[int(14/dt)],v0[int(14/dt)],s=30,edgecolor='k',facecolor='w',marker='s',zorder=-1)

        ax.plot(IL6,e,v,linewidth=2,c=ggcolors[mode-1],label='Treated')
        ax.scatter(IL6[0],e[0],v[0],marker='^',s=40,edgecolor='k',facecolor='w')
        ax.scatter(IL6[int(7/dt)],e[int(7/dt)],v[int(7/dt)],s=30,edgecolor='k',facecolor='w',marker='o')
        ax.scatter(IL6[int(14/dt)],e[int(14/dt)],v[int(14/dt)],s=30,edgecolor='k',facecolor='w',marker='s')
        '''
        ei-ea-v
        '''
        ax2.plot(ei0,ea0,np.zeros(len(v0)),linewidth=1,linestyle='-',c=ggcolors[mode-1],alpha=0.5)
        ax2.plot(ei,ea,np.zeros(len(v)),linewidth=1,c=ggcolors[mode-1],alpha=0.5)
        ax2.plot(ei0,ea0,v0,linewidth=2,linestyle='-',dashes=(5,0.6),c=ggcolors[mode-1],label='No treatment')
        ax2.scatter(ei0[int(7/dt)],ea0[int(7/dt)],v0[int(7/dt)],s=30,edgecolor='k',facecolor='w',marker='o',zorder=-1)
        ax2.scatter(ei0[int(14/dt)],ea0[int(14/dt)],v0[int(14/dt)],s=30,edgecolor='k',facecolor='w',marker='s',zorder=-1)

        ax2.plot(ei,ea,v,linewidth=2,c=ggcolors[mode-1],label='Treated')
        ax2.scatter(ei[0],ea[0],v[0],marker='^',s=40,edgecolor='k',facecolor='w')
        ax2.scatter(ei[int(7/dt)],ea[int(7/dt)],v[int(7/dt)],s=30,edgecolor='k',facecolor='w',marker='o')
        ax2.scatter(ei[int(14/dt)],ea[int(14/dt)],v[int(14/dt)],s=30,edgecolor='k',facecolor='w',marker='s')
        '''
        adjusting main text figure
        '''
        if mode!=4:
            ax.scatter(IL6[-1],e[-1],v[-1],marker='^',s=40,edgecolor='k',facecolor='w',zorder=100)
            ax2.scatter(ei[-1],ea[-1],v[-1],marker='^',s=40,edgecolor='k',facecolor='w',zorder=100)
        ax.tick_params(labelsize=11,pad=-2)
        ax2.tick_params(labelsize=11,pad=-2)
        if mode==4:
            ax.plot([IL60[-1], IL60[-1]], [e0[-1], e0[-1]], [v0[-1], 0], color='k', linewidth=1, alpha=1)
            ax.scatter(IL60[-1],e0[-1],v0[-1],marker='X',s=30,edgecolor='k',facecolor='w',label=r'[nCoV]$\neq$0',zorder=10)
            ax2.plot([ei0[-1], ei0[-1]], [ea0[-1], ea0[-1]], [v0[-1], 0], color='k', linewidth=1, alpha=1)
            ax2.scatter(ei0[-1],ea0[-1],v0[-1],marker='X',s=30,edgecolor='k',facecolor='w',label=r'[nCoV]$\neq$0',zorder=10)
        if mode==2:
            ax.set_xticks([0,500,1000])
            ax.set_yticks([0,3,6,9])
            ax.set_zticks([0,1000,2000])
            ax.set_zlim([0,2000])
            ax.view_init(azim=156, elev=29)
            ax2.set_zticks([0,1000,2000])
            ax2.set_zlim([0,2000])
            ax2.view_init(azim=156, elev=29)
        if mode==3:
            ax.set_xticks([0,1000,2000])
            ax.set_yticks([0,3,6])
            ax.set_zticks([0,1000,2000,3000])
            ax.set_zticklabels([0, 1, 2, r'$\geq 3$'])
            ax.set_zlim([0,3100])
            ax.view_init(azim=151, elev=27)
            ax2.set_zticks([0,1000,2000,3000])
            ax2.set_zticklabels([0, 1, 2, r'$\geq 3$'])
            ax2.view_init(azim=151, elev=27)
        if mode==4:
            ax.set_xticks([0,1000,2000])
            ax.set_yticks([0,1,2])
            ax.set_zticks([0, 1000, 2000, 3000])
            ax.set_zlim([0,3100])
            ax.set_yticklabels([0,1,2])
            ax.set_zticklabels([0, 1, 2, r'$\geq 2$'])
            ax.view_init(azim=168,elev=21)
            ax2.set_xticks([0,0.5])
            ax2.set_yticks([0,1,2])
            ax2.set_zticks([0, 1000, 2000, 3000])
            ax2.set_zlim([0,3100])
            ax2.set_zticklabels([0, 1, 2, r'$\geq 3$'])
            ax2.view_init(azim=144,elev=26)
        sFormatter1=matplotlib.ticker.ScalarFormatter(useOffset=False, useMathText=True)
        sFormatter1.set_powerlimits((-1, 2))
        ax.xaxis.set_major_formatter(sFormatter1)
        ax.xaxis.offsetText.set_visible(False)
        ax.xaxis.offsetText.set_position((1, -0.2))
        if mode==2:
            ax.zaxis.set_major_formatter(sFormatter1)
            ax.zaxis.offsetText.set_visible(False)
            ax.zaxis.offsetText.set_position((0.5, 0.5, 0.5))
            ax2.zaxis.set_major_formatter(sFormatter1)
            ax2.zaxis.offsetText.set_visible(False)
            ax2.zaxis.offsetText.set_position((0.5, 0.5, 0.5))
        ax.text2D(0.8,0.85,'Mode %i'%(mode),transform=ax.transAxes,fontsize=11)
        ax.text2D(0.8,0.78,r'$\Delta Q/Q=%1.2f$'%Scrs[mode-2],transform=ax.transAxes,fontsize=9)#[~indices]
        ax2.text2D(0.8,0.85,'Mode %i'%(mode),transform=ax2.transAxes,fontsize=11)
        ax2.text2D(0.8,0.78,r'$\Delta Q/Q=%1.2f$'%Scrs[mode-2],transform=ax2.transAxes,fontsize=9)#[~indices]
        ax.set_xlabel(r'[IL-6](10$^3$pg/mL)',labelpad=-1,fontsize=12)
        ax.set_ylabel(r'$\epsilon$',labelpad=-1,fontsize=13)
        ax2.set_xlabel(r'$\epsilon_i$',labelpad=2,fontsize=13)
        ax2.set_ylabel(r'$\epsilon_a$',labelpad=2,fontsize=13)
        if mode==2:
            ax.set_zlabel(r'[nCoV]($10^8/mL$)', labelpad=-2, fontsize=12)
            ax2.set_zlabel(r'[nCoV]($10^8/mL$)', labelpad=-2, fontsize=12)
        else:
            ax.set_zlabel(r'[nCoV]($10^9/mL$)',labelpad=-2,fontsize=12)
            ax2.set_zlabel(r'[nCoV]($10^9/mL$)',labelpad=-2,fontsize=12)
        fig.subplots_adjust(left=0.05,right=0.85,bottom=0.1,top=1)
        fig2.subplots_adjust(left=0.05,right=0.85,bottom=0.1,top=1)
        fig.savefig('MedTraj%i.svg'%(mode))
        fig.savefig('MedTraj%i.png'%(mode))
        fig2.savefig('MedTraj%i2.svg'%(mode))
        fig2.savefig('MedTraj%i2.png'%(mode))
    np.savetxt('Scrs.txt',Scrs)
    plt.show()
    pool.close()