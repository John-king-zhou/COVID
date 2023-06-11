import numpy as np
import multiprocessing
from Equation import func
from scipy.integrate import odeint
from E_Calculation import E_kill,E_clear

def initial_from_param(Para, Bm=0, T8m=0, T4m=0, A0=0, is_postvax=1):
    if is_postvax:
        initial = [0.01, 0, Para[52] / Para[62], 0, 0, Para[53]/Para[65], 0, 0, Para[159], 0, 0, 0, 0, 0, 0,
           Para[54]/Para[70], 0, Para[160], 0, 0, 0, 0, 0, 0, Para[77] / Para[101], Para[82] / Para[102],
           Para[84]/Para[103], (Para[88]+Para[54]/Para[70]*Para[90])/Para[104], Para[91]/Para[105], Para[95]/Para[106],
           0, 0]
        Ig0 = (1 + Para[51] * initial[25] / (initial[25] + Para[116])) * Bm * Para[100] / Para[107]
        initial[31] = A0
        initial[30] = Ig0
        initial[23] = Bm
        initial[20] = T8m
        initial[16] = T4m
    else:
        initial = [0.01, 0, Para[52] / Para[62], 0, 0, Para[53]/Para[65], 0, 0, Para[159], 0, 0, 0, 0, 0, 0,
           Para[54]/Para[70], 0, Para[160], 0, 0, 0, 0, 0, 0, Para[77] / Para[101], Para[82] / Para[102],
           Para[84]/Para[103], (Para[88]+Para[54]/Para[70]*Para[90])/Para[104], Para[91]/Para[105], Para[95]/Para[106],
           0, 0]
    return initial

def get_Bm(Ab,Para):
    initial = [0.01, 0, Para[52] / Para[62], 0, 0, Para[53] / Para[65], 0, 0, Para[159] / 2, 0, 0, 0, 0, 0, 0,
               Para[54] / Para[70], 0, Para[160] / 2, 0, 0, 0, 0, 0, 0, Para[77] /
               Para[101], Para[82] / Para[102],
               Para[84] / Para[103], (Para[88] + Para[54] / Para[70] * Para[90]) /
               Para[104], Para[91] / Para[105], Para[95] / Para[106],
               0, 0]
    Bm = Ab / (1 + Para[51] * initial[25] / (initial[25] + Para[116])) * Para[107] / Para[100]
    return Bm

def get_results(Para,v,CD4Tm,CD8Tm,Ab,type):
    Time = np.arange(0, 20, 0.1)
    A=1
    initial = [v, 0, Para[52] / Para[62], 0, 0, Para[53] / Para[65], 0, 0, 0, 0, 0, 0, 0, 0, 0,
               Para[54] / Para[70], 0, 0, 0, 0, 0, 0, 0, 0, Para[77] / Para[101], Para[82] / Para[102],
               Para[84] / Para[103], (Para[88] + Para[54] / Para[70] * Para[90]) / Para[104], Para[91] / Para[105],
               Para[95] / Para[106], 0, 0]
    if type=='A':
        Bm=0
    elif type=='B':
        Bm=get_Bm(Ab,Para)
    initial[31] = A
    initial[30] = Ab
    initial[23] = Bm
    initial[20] = CD8Tm
    initial[16] = CD4Tm
    results = odeint(func, initial, Time, args=(Para,))
    ek=E_kill(results,Para)
    ec=E_clear(results,Para)
    results = np.vstack((results.T, ec, ek)).T
    return results

def find_edge_line(threshold, map):
    x_indices = []
    y_indices = []
    for i in range(len(map[:, 0])):
        i_seq = np.array(map[i, :])-threshold
        for j in range(1, len(i_seq)):
            if i_seq[j]*i_seq[j-1] < 0:
                x_indices.append(j)
                y_indices.append(i)
                break
    return x_indices, y_indices

def find_edge(type, Para, v0, CD4Tm, T8m_seq, Ab_seq, mode):
    # find simulation boaders
    Ab_bd_indices=[]
    ek=[]
    ec=[]
    N=len(T8m_seq)
    for i in reversed(range(N)):
        start=0
        CD8Tm=T8m_seq[i]
        if len(Ab_bd_indices)!=0:
            start=Ab_bd_indices[-1]
        Ab_bd=len(Ab_seq)-1
        for j in range(start,N):
            Ab=Ab_seq[j]
            results=get_results(Para,v0,CD4Tm,CD8Tm,Ab,type=type)
            vt=results[:,0]
            if np.max(vt)<=vt[0]:
                Ab_bd=j
                ec.append(results[0,32])
                ek.append(results[0,33])
                break
        Ab_bd_indices.append(Ab_bd)
    return Ab_seq[Ab_bd_indices],np.array(ec),np.array(ek),mode

def find_edge_merge(args):
    return find_edge(*args)

if __name__=='__main__':
    for type in ['A','B']:
        CD4Tm=0.02
        N2=50
        T8m_seq=np.linspace(0,1,N2)
        Ab_seq=np.linspace(0,2000,N2)
        # Boader
        # simlation boader
        Ab_bd = []
        BD_ek = []
        BD_ec = []

        mode_save =[]

        cores = multiprocessing.cpu_count()
        pool = multiprocessing.Pool(processes=int(cores/3))
        for mode in range(1,5,1):
            print('mode %i started'%mode)
            Paras=np.loadtxt('Spl_Para%i.txt'%mode)
            args=[]
            for j in range(Paras.shape[0]):
                args.append([type,Paras[j,:],0.01,CD4Tm,T8m_seq,Ab_seq,mode])
            iterations=len(args)
            count = 0
            # find_edge_merge(args[0]) # testing
            for Ab,ec,ek,mode in pool.imap(find_edge_merge,args):
                Ab_bd.append(Ab)
                BD_ec.append(ec)
                BD_ek.append(ek)
                mode_save.append(mode)

                ratio = count / iterations
                rat_str = ['>'] * int(ratio * 50) + ['-'] * (50 - int(ratio * 50))
                rat_str = ''.join(rat_str)
                print('\r' + rat_str + '%.2f %%' % (ratio * 100), end='')
                count += 1

        dt_save={'bd_Bm':Ab_bd, 'bd_ec':BD_ec,'bd_ek':BD_ek,'mode':mode_save,'CD8Tm':T8m_seq}
        np.save('3types_border_data%s.npy'%type, dt_save)
        pool.close()