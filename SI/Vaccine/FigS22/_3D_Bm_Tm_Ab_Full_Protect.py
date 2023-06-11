import numpy as np
import matplotlib.pyplot as plt
import multiprocessing

from Equation import func
from scipy.integrate import odeint
from E_Calculation import E,E1,E2,E12
import os
from scipy import interpolate
from PIL import Image

set2colors = ['#fc8d62', '#66c2a5', '#a6d854',
              '#e78ac3', '#ffd92f', '#e5c494', '#b3b3b3']
ggcolors=['#808080','#1F77B4','#FF7F0E','#D62728','#4DBEEE','#77AC30','#9467BD']
import matplotlib.mathtext as mathtext
mathtext.FontConstantsBase.sup1 = 0.4
time0 = np.arange(0, 50, 0.1)

ls=['-','--']
colorseq = [set2colors[1], set2colors[0], set2colors[3]]
path = os.path.split(os.path.abspath(__file__))[0]

def fig_3d_gif(fig, ax, angle_1st, name):
    file = []
    # ax.legend(loc='lower left', bbox_to_anchor=[0.9, 0.4],
    #           frameon=False, labelspacing=0.5, ncol=1)
    for angle in range(0, 360, 3):
        ax.view_init(angle_1st, angle)
        file.append(path+'/3d_tmp/3d_vis_'+str(angle)+'.png')
        fig.savefig(file[-1], dpi=150)

    frame = []
    for f in file:
        new_frame = Image.open(f)
        frame.append(new_frame)

    # Save into a GIF file that loops forever
    frame[0].save(path+'/3d_vis'+name+'.gif', format='GIF',
                  append_images=frame[1:],
                  save_all=True,
                  duration=80, loop=0)

def Bm_defined_Ab(Para,Bm):
    initial = [0.01, 0, Para[52] / Para[62], 0, 0, Para[53]/Para[65], 0, 0, Para[159]/2, 0, 0, 0, 0, 0, 0,
               Para[54]/Para[70], 0, Para[160]/2, 0, 0, 0, 0, 0, 0, Para[77] /
               Para[101], Para[82] / Para[102],
               Para[84]/Para[103], (Para[88]+Para[54]/Para[70]*Para[90]) /
               Para[104], Para[91]/Para[105], Para[95]/Para[106],
                0, 0]
    return (1 + Para[51] * initial[25] / (initial[25] + Para[116])) * Bm * Para[100] / Para[107]

def get_time_seq(Para, v0, CD8Tm, CD4Tm, Bm, Ab, A):

    if A != 0:
        initial = [0.01, 0, Para[52] / Para[62], 0, 0, Para[53]/Para[65], 0, 0, Para[159]/2, 0, 0, 0, 0, 0, 0,
                   Para[54]/Para[70], 0, Para[160]/2, 0, 0, 0, 0, 0, 0, Para[77] /
                   Para[101], Para[82] / Para[102],
                   Para[84]/Para[103], (Para[88]+Para[54]/Para[70]*Para[90]) /
                   Para[104], Para[91]/Para[105], Para[95]/Para[106],
                   0, 0]
        
        initial[31] = A
        initial[30] = Ab
        initial[23] = Bm
        initial[20] = CD8Tm
        initial[16] = CD4Tm
        initial[0] = v0
    elif A == 0:
        initial = [0.01, 0, Para[52] / Para[62], 0, 0, Para[53]/Para[65], 0, 0, Para[159], 0, 0, 0, 0, 0, 0,
                   Para[54]/Para[70], 0, Para[160], 0, 0, 0, 0, 0, 0, Para[77] /
                   Para[101], Para[82] / Para[102],
                   Para[84]/Para[103], (Para[88]+Para[54]/Para[70]*Para[90]) /
                   Para[104], Para[91]/Para[105], Para[95]/Para[106],
                   0, 0]

    results = odeint(func, initial, time0, args=(Para,))
    e = E(results, Para)
    e1 = E1(results, Para)
    e1 = E1(results, Para)
    e2 = E2(results, Para)
    e12 = E12(results, Para)
    results = np.vstack((results.T, e, e1, e2, e12)).T
    return results

def find_Tm_with_fixed_bm_varied_Ab(Para,Bm,Ab_seq,v0,CD4Tm,A):

    Tm_seq=np.zeros(Ab_seq.shape)
    dTm_min=0.02
    for Ab,i in zip(Ab_seq,range(len(Ab_seq))):
        Tm_tmp=[]
        max_v=[]
        dTm=0+dTm_min
        if i!=0:
            Tm_tmp.append(Tm_seq[i-1])
        else:
            Tm_tmp.append(Tm_seq[i]+.1)
        sol = get_time_seq(Para, v0, Tm_tmp[-1], CD4Tm, Bm, Ab, A)
        max_v.append(max(sol[:,0]))

        Tm_tmp.append(Tm_tmp[-1])
        
        # if i==0:
        #     while len(Tm_tmp) < 100:
        #         sol = get_time_seq(Para, v0, Tm_tmp[-1], CD4Tm, Bm, Ab, A)
        #         max_v.append(max(sol[:,0]))

        #         # Condition for output
        #         if  max_v[-1] - v0 < 10 and max_v[-1] > v0:
        #             break
                
        #         # change learning rate
        #         if max_v[-2]>v0 and max_v[-2]-max_v[-1]>10: # speed up when max_v-vo_0<10
        #             dTm = 0.5*abs(max_v[-1]-v0)/(abs(max_v[-1]-max_v[-2])/abs(Tm_tmp[-1]-Tm_tmp[-2]))
        #             dTm = max([dTm, dTm_min])
        #         elif max_v[-2]==v0:
        #             dTm = 0 + dTm_min
        #         dTm = min([dTm,Tm_tmp[-1]/10+dTm_min])
        #         # update new Tm
        #         if max_v[-2]>v0 and max_v[-1]>v0:
        #             Tm_tmp.append(max(0, Tm_tmp[-1]+dTm))
        #         else:
        #             Tm_tmp.append(max(0, Tm_tmp[-1]/2))
        # else:
        if i==0:
            dTm=0.05

            for j in range(100):
                sol = get_time_seq(Para, v0, Tm_tmp[-1], CD4Tm, Bm, Ab, A)
                max_v.append(max(sol[:,0]))
                if  (max_v[-1]==v0 and max_v[-2]>v0) \
                    or (max_v[-1]>v0 and max_v[-2]==v0):
                    break
                if max_v[-1]==v0:
                    Tm_tmp.append(Tm_tmp[-1]-dTm)
                elif max_v[-1]>v0:
                    Tm_tmp.append(Tm_tmp[-1]+dTm)
            dTm=0+dTm_min

        for j in range(100):
            sol = get_time_seq(Para, v0, Tm_tmp[-1], CD4Tm, Bm, Ab, A)
            max_v.append(max(sol[:,0]))

            if  (max_v[-1]==v0 and max_v[-2]>v0) \
                or (max_v[-1]>v0 and max_v[-2]==v0):
                break
            if max_v[-1]==v0:
                Tm_tmp.append(Tm_tmp[-1]-dTm)
            elif max_v[-1]>v0:
                Tm_tmp.append(Tm_tmp[-1]+dTm)

        Tm_seq[i] = Tm_tmp[-1]

    return Tm_seq


def get_results_merge(args):
    return find_Tm_with_fixed_bm_varied_Ab(*args)


def prt_proc(ratio):
    tot_len=80
    rat_str = ['>'] * int(ratio * tot_len) + ['-'] * (tot_len - int(ratio * tot_len))
    position=int(tot_len*np.random.rand())
    str_par = '|{:.2f}%|'.format(ratio*100)
    str_par = list(str_par)
    rat_str[position:position+len(str_par)]=str_par
    rat_str = ''.join(rat_str)
    print('\r' + rat_str, end='')


if __name__=='__main__':
    # select Mode
    mode = 3
    Parameter=np.loadtxt(path+"\GeoMean_para_Mode%i.txt" %mode)

    '''
    Initial conditions
    '''
    A=1
    v0=0.01
    cd4tm=0.02

    # for full protection
    Bmseq=np.linspace(0,0.5,20)
    # CD8Tmseq=np.linspace(0,1,20)
    Abseq=np.linspace(100,1500,15)

    try:
        dt = np.load(path+'\\3D_AbBmTm_protect_m%d.npy' %mode, allow_pickle=True).item()
        Bmseq, Abseq = dt['bm'],dt['Ab']
    except:
        cores=multiprocessing.cpu_count()
        pool=multiprocessing.Pool(processes=int(cores)-2)
        args=[]
        dt={'Para':Parameter,'bm':Bmseq,'Ab':Abseq,'Tm':[]}

        for bmi in Bmseq:
            args.append([Parameter,bmi,Abseq,v0,cd4tm,A])

        count=0
        for results in pool.imap(get_results_merge, args):
            count += 1
            ratio = count / len(args)
            prt_proc(ratio)

            dt['Tm'].append(results)

        np.save(path+'\\3D_AbBmTm_protect_m%d.npy' %mode, dt)
        pool.close()


    # ------------ 3D plot of Bm, Ab, Tm ----------------
    XX,YY = np.meshgrid(Bmseq,Abseq)
    ZZ_Tm = dt['Tm'][0]
    for seq in dt['Tm'][1:]:
        ZZ_Tm =np.vstack((ZZ_Tm,seq))
    ZZ_Tm = ZZ_Tm.T
    
    fig = plt.figure()
    ax3d = fig.add_subplot(projection='3d')

    ZZ_plot=np.copy(ZZ_Tm)
    ZZ_plot[Bm_defined_Ab(Parameter, XX)>YY]=np.NaN
    # ax3d.plot_wireframe(XX,YY,ZZ_plot,rstride=1,cstride=1,linewidth=.5)
    ax3d.plot_surface(XX, YY, ZZ_plot, alpha=0.3, color='green')

    # for j in [0,4,9]:
    #     ax3d.plot(XX[j,:],YY[j,:],ZZ_Tm[j,:],
    #               lw=2,
    #               label='%i'%Abseq[j])

    # Bm induced Ab
    Bm_axis = np.arange(min(Bmseq),max(Bmseq),0.01)
    Ab_by_Bm = Bm_defined_Ab(Parameter, Bm_axis)
    Bm_tmp = Bm_axis[Ab_by_Bm<max(Abseq)]
    Ab_by_Bm = Ab_by_Bm[Ab_by_Bm<max(Abseq)]
    Bm_tmp = Bm_tmp[Ab_by_Bm>min(Abseq)]
    Ab_by_Bm = Ab_by_Bm[Ab_by_Bm>min(Abseq)]

    f_zz = interpolate.interp2d(XX.flat,YY.flat,ZZ_Tm, kind='linear')
    zz_new = f_zz(Bm_tmp,Ab_by_Bm)
    zz_line = np.array([zz_new[i,i] for i in range(len(zz_new))])

    ax3d.plot(Bm_tmp,Ab_by_Bm,zz_line,
              lw=2.5, alpha=0.8, ls='--',
              c='red',label=r'$Steady State$')

    ax3d.legend(title=r'$[Ab]_{0}$',
                loc='lower left', 
                bbox_to_anchor=[0.9, 0.4],
                frameon=True)
    ax3d.set_xlabel(r'$[Bm]_{0}$')
    ax3d.set_ylabel(r'$[Ab]_{0}$')
    ax3d.set_zlabel(r'$[CD8+Tm]_{0}$')
    ax3d.set_title('Full Proteciton Surface')
    ax3d.view_init(30, 49)

    # fig_3d_gif(fig, ax3d, 30, 'Protection_Surface')

    plt.show()
