import numpy as np
import multiprocessing
from scipy.integrate import odeint
from Latin_Hypercube import LHSample
from Equation import func
from Type_Characterization2 import Type_Characterization, WithinPhys
from E_Calculation import E
import time
import pandas as pd


# warnings.filterwarnings('error')
def Hill(x, k, n):
    if k == 0:
        return 1
    else:
        return x ** n / (k ** n + x ** n)


def score(results):
    vf=results[-1,0]
    hm=np.min(results[:,2])
    IL6m=np.max(results[:,26])
    Hss=results[0,2]

    IL6c = 2000
    hc = 30 / 50 * Hss
    vc = 1
    S = (1 + 1 * Hill(hm, hc, 1)) * (1 + 2 * Hill(IL6c, IL6m, 1)) * (1 + 1 * Hill(vc, vf, 1))
    return S


def get_mode(Para):
    Time = np.arange(0, 50, 0.1)
    initial = [0.01, 0, Para[52] / Para[62], 0, 0, Para[53] / Para[65], 0, 0, Para[159], 0, 0, 0, 0, 0, 0,
               Para[54] / Para[70], 0, Para[160], 0, 0, 0, 0, 0, 0, Para[77] / Para[101], Para[82] / Para[102],
               Para[84] / Para[103], (Para[88] + Para[54] / Para[70] * Para[90]) / Para[104], Para[91] / Para[105],
               Para[95] / Para[106], 0, 0]
    result = odeint(func, initial, Time, args=(Para,))
    e = E(result,Para)
    result=np.vstack((result.T,e)).T
    if WithinPhys(result):
        mode = Type_Characterization(result)
        return mode, Para, result
    else:
        return 100, Para, result

sample_indices1 = np.arange(0,25,1)
sample_indices1 = np.delete(sample_indices1,[8,9,11,12,13,14,15,17,18,19,21,22,23,10,20,24])
sample_indices2 = np.hstack((np.array([8,9,11,12,13,14,15,17,18,19,21,22,23]),np.arange(36,52,1)))
sample_indices3 = [150,159,160]+list(np.arange(128,132,1))+list(np.arange(27,36,1))
sample_indices4 = np.array([25,61,151,52])
para = np.loadtxt('Para.txt')
n = len(para)
LogParaBound = []
for i in range(len(para)):
    a = np.log10(para[i])
    if i in sample_indices1:
        LogParaBound.append([a - np.log10(5), a + np.log10(5)])
    elif i in sample_indices2:
        LogParaBound.append([a - np.log10(2), a + np.log10(2)])
    elif i in sample_indices3:
        LogParaBound.append([a - np.log10(5), a + np.log10(5)])
    elif i in sample_indices4:
        LogParaBound.append([a - np.log10(5), a + np.log10(5)])
    else:
        LogParaBound.append([a, a])
LogParaBound = np.array(LogParaBound)

if __name__ == '__main__':
    cores = multiprocessing.cpu_count()
    pool = multiprocessing.Pool()
    N = 1 * 10 ** 4
    time0 = time.perf_counter()
    LogPara_matrix = np.array(LHSample(len(LogParaBound), LogParaBound, N))
    Para_matrix = np.power(10, LogPara_matrix)
    Para_mode = np.zeros(5)
    #data = dict.fromkeys(['k_infect', 'N1', 'H0', 'd_If', 'mode', 'vmax', 'IL6max', 'Q'], [])
    data = {'k_infect':[],'N1':[],'H0':[],'d_If':[],'gamma':[],'mode':[],'vmax':[],'IL6max':[],'Q':[],'emax':[],'H':[],'e':[]}
    num=0
    for y, para, result in pool.imap(get_mode, Para_matrix):
        Para_mode[0] += 1
        data['k_infect'] += [para[25], ]
        data['N1'] += [para[151], ]
        data['H0'] += [para[52] / para[62], ]
        data['d_If'] += [para[61], ]
        data['gamma'] += [para[25] * para[151] * para[61] * para[52] / para[62], ]
        v = result[:, 0]
        e = result[:, -2]
        data['vmax'] += [np.max(v), ]
        data['IL6max'] += [np.max(result[:, 26]), ]
        data['Q'] += [score(result), ]
        data['emax'] += [np.max(e), ]
        data['H'] += [result[np.argmax(v),2]/result[0,2],]
        data['e'] += [e[np.argmax(v)],]
        if 0<y<=4:
            data['mode'] += [y,]
        else:
            data['mode'] += [5,]
        if y <= 4:
            Para_mode[y] += 1
        num += 1

        # print
        ratio = num / N
        rat_str = ['>'] * int(ratio * 50) + ['-'] * (50 - int(ratio * 50))
        rat_str = ''.join(rat_str)
        print('\r' + rat_str + '%.2f %%' % (ratio * 100), end='')
        if num % 1000 == 0:
            print('processes completed:', int(num), 'Within_Phys_Range:', Para_mode[0],
                  'Mode 1:%i/%i' % (Para_mode[1], Para_mode[0]),
                  'Mode 2:%i/%i' % (Para_mode[2], Para_mode[0]),
                  'Mode 3:%i/%i' % (Para_mode[3], Para_mode[0]),
                  'Asymptomatic:%i/%i' % (Para_mode[4], Para_mode[0]),
                  'time used:%1.1f s' % (time.perf_counter() - time0))
    pool.close()
    df=pd.DataFrame(data)
    df.to_csv('Virulence_Sample.csv')
    del (data)

