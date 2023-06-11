import numpy as np
import multiprocessing
from scipy.integrate import odeint
from Latin_Hypercube import LHSample
from Equation import func
from Type_Characterization2 import *
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

def change(Para,v):
    IFNI = [0,1]
    Innate = [0,1,5,6,7]
    # Innate = [0,1]
    AP =  [8,18,21]
    if v=='SARS':
        Para[:,Innate]*=10
        Para[:,AP]*=0.5
        return Para
    elif v=='IAV':
        Para[:, 52]*=0.7
        Para[:,IFNI]*=1.5
        return Para
    else:
        return Para


def get_mode(Para):
    Time = np.arange(0, 80, 0.1)
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



if __name__ == '__main__':
    cores = multiprocessing.cpu_count()
    pool = multiprocessing.Pool()
    time0 = time.perf_counter()
    Para_mode = np.zeros(5)
    data = {'virus':[],'mode':[],'vmax':[],'IL6max':[],'emax':[],'Q':[]}
    N=10**3
    for vtype in ['SARS-CoV-2', 'SARS', 'IAV']:
        # LogParaBound = np.loadtxt('logparabound.txt')
        # LogPara_matrix = np.array(LHSample(len(LogParaBound), LogParaBound, N))
        # Para0 = np.power(10, LogPara_matrix)
        ts=[[] for i in range(3)]
        Para0 = np.vstack([np.loadtxt('Spl_Para%i.txt'%i) for i in range(1,5,1)])
        N = Para0.shape[0]
        print(N)
        Para = change(Para0, vtype)
        Para_mode = np.zeros(5)
        num = 0
        for y, para, result in pool.imap(get_mode, Para):
            if y <= 4:
                Para_mode[0] += 1
                if y != 0:
                    Para_mode[y] += 1
                    data['virus'] += [vtype, ]
                    data['mode'] += [y, ]
                    v = result[:, 0]
                    IL6 = result[:, 26]
                    e = result[:, -1]
                    data['vmax'] += [np.max(v), ]
                    data['IL6max'] += [np.max(IL6), ]
                    data['emax'] += [np.max(e), ]
                    data['Q'] += [score(result), ]
                    ts[0].append(v)
                    ts[1].append(IL6)
                    ts[2].append(e)
            num += 1
            # print
            ratio = num / N
            rat_str = ['>'] * int(ratio * 50) + ['-'] * (50 - int(ratio * 50))
            rat_str = ''.join(rat_str)
            print('\r' + rat_str + '%.2f %%' % (ratio * 100), end='')
            if num % 1000 == 0:
                print('virus type: %s, processes completed:' % vtype, int(num), 'Within_Phys_Range:', Para_mode[0],
                      'Mode 1:%i/%i' % (Para_mode[1], Para_mode[0]),
                      'Mode 2:%i/%i' % (Para_mode[2], Para_mode[0]),
                      'Mode 3:%i/%i' % (Para_mode[3], Para_mode[0]),
                      'Mode 4:%i/%i' % (Para_mode[4], Para_mode[0]),
                      'time used:%1.1f s' % (time.perf_counter() - time0))
        ts=np.array(ts)
        np.save('%s_ts.npy'%vtype,ts)
    df = pd.DataFrame(data)
    df.to_csv('Virus_Sample.csv')
    del (data)
    pool.close()

