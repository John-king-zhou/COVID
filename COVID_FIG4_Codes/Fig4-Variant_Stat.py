#Fig 4B, C, D
import numpy as np
import multiprocessing
from scipy.integrate import odeint
from Equation import func
from Type_Characterization2 import Type_Characterization, WithinPhys
from E_Calculation import E
import time
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

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
    if v=='alpha':
        Para[:,25]*=1.3
        Para[:,52]*=1.3
        return Para
    elif v=='delta1':
        Para[:,25]*=1.6
        Para[:,52]*=1.5
        Para[:,151]*=1.2
        return Para
    else:
        return Para

def get_mode(Para):
    Time = np.arange(0, 80, 0.1)
    initial = [Para[-1], 0, Para[52] / Para[62], 0, 0, Para[53] / Para[65], 0, 0, Para[159], 0, 0, 0, 0, 0, 0,
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

if __name__=='__main__':
    try:
        df=pd.read_csv('Variant_Sample2.csv')
    except:
        cores = multiprocessing.cpu_count()
        pool = multiprocessing.Pool()
        time0 = time.perf_counter()
        data = {'gamma': [], 'mode': [], 'vmax': [], 'IL6max': [], 'Q': [], 'emax': [], 'v4': [], 'vtype': []}
        vt = []
        for vtype in ['SARS-CoV-2', 'alpha', 'delta1']:
            Para0 = np.vstack([np.loadtxt('Spl_Para%i.txt' % mode) for mode in range(1, 5, 1)])
            N=Para0.shape[0]
            Para=change(Para0,vtype)
            Para_matrix = np.hstack((Para, np.ones((N, 1)) * 0.01))
            Para_mode = np.zeros(5)
            num = 0
            for y, para, result in pool.imap(get_mode, Para_matrix):
                if y <= 4:
                    Para_mode[0] += 1
                    if y != 0:
                        Para_mode[y] += 1
                    data['gamma'] += [para[25] * para[151] * para[61] * para[52] / para[62], ]
                    data['mode'] += [y, ]
                    v = result[:, 0]
                    e = result[:, -1]
                    data['vmax'] += [np.max(v), ]
                    data['IL6max'] += [np.max(result[:, 26]), ]
                    data['Q'] += [score(result), ]
                    data['emax'] += [np.max(e), ]
                    data['v4'] += [v[40], ]
                    data['vtype'] += [vtype, ]
                    vt.append(v)
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
                          'Asymptomatic:%i/%i' % (Para_mode[4], Para_mode[0]),
                          'time used:%1.1f s' % (time.perf_counter() - time0))
        df = pd.DataFrame(data)
        vt = np.array(vt)
        np.savetxt('Variant_vt2.txt', vt)
        df.to_csv('Variant_Sample2.csv')
        del (data)
        pool.close()

    ggcolors = ['#2CA02C', '#1F77B4', '#FF7F0E', '#D62728', ]
    set2colors = ['#66c2a5', '#fc8d62', '#e78ac3', '#ffd92f', '#e5c494', '#b3b3b3']
    markers = ['o', '^', 's', 'd']
    flierprops = dict(marker='x', markerfacecolor='k', markersize=3,
                      linestyle='none', markeredgecolor='k')
    labels = ['Mode 1', '2', '3', '4']
    vtype = ['SARS-CoV-2', 'alpha', 'delta1']
    xlabels = ['nCoV', 'Alpha\nB.1.1.7', 'Delta\nB.1.617.2']
    xlabel_simple = ['WT', 'Alpha', 'Delta']
    xlabel_greek = ['19', r'$\alpha$', r'$\delta$']

    df = pd.read_csv('Variant_Sample2.csv')
    df['lgv4'] = np.log10(df.loc[:, 'v4'])
    df['lgvmax'] = np.log10(df.loc[:, 'vmax'])
    # mode 1/2/3/4 ratio
    fig1, ax1 = plt.subplots(1, 1, figsize=(3, 1.5))
    m = np.zeros((4, len(vtype)))
    for v in range(len(vtype)):
        vt = vtype[v]
        inds = df.index[df['vtype'] == vt].tolist()
        print(vt,len(inds))
        df_v = df.iloc[inds, :]
        for mode in [1, 2, 3, 4]:
            indices = df_v.index[df_v['mode'] == mode].tolist()
            m[mode-1, v] = len(indices)
    y = np.zeros(len(vtype))
    for i in range(4):
        x = m[i] / np.sum(m, axis=0)
        print(labels[i], 'nCoV:', x[0], 'Alpha:', x[1], 'Delta:', x[2])
        ax1.barh(np.arange(0, 3, 1), x + y, color=ggcolors[i], zorder=-10 * i)
        y += x
    ax1.set_xticks([0, 0.5, 1])
    ax1.set_xlim([0, 1])
    ax1.set_ylabel('')
    ax1.set_yticks([0, 1, 2])
    ax1.set_yticklabels(xlabel_simple, fontsize=10)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    fig1.subplots_adjust(bottom=0.18, right=0.95, left=0.15)
    fig1.savefig('Variant_Ratio2.svg')
    fig1.savefig('Variant_Ratio2.png')
    fig0, ax0 = plt.subplots(1, 1)

    fig2, ax2 = plt.subplots(1, 1, figsize=(3, 1))
    mean_lgv4=[]
    std_lgv4=[]
    for i in range(3):
        vt = vtype[i]
        inds = df.index[df['vtype'] == vt].tolist()
        print(vt,len(inds))
        df_v = df.iloc[inds, :]
        mean = np.mean(df_v['lgv4'])
        std = np.std(df_v['lgv4'])
        mean_lgv4.append(mean)
        std_lgv4.append(std)
    sns.barplot(data=df, y='vtype', x='lgv4', ax=ax2, ci='sd',palette=set2colors,edgecolor='k', errcolor='k', errwidth=1, capsize=0.3, linewidth=1, orient='h')

    y = [1, 1.8, 1]
    gamma = [3.6, 6.1, 10.4]
    for i in range(3):
        c=set2colors[i]
        if i==2:
            c='k'
        ax2.text(y=i, x=y[i], s='$\gamma$=%1.1f' % gamma[i], color=c, ha='left', va='center')
    ax2.set_xlabel('$log_{10}nCoV$')
    ax2.set_xticks([1, 2])
    ax2.set_yticklabels(xlabel_simple, fontsize=10)
    ax2.set_ylabel('')
    fig2.subplots_adjust(bottom=0.2, right=0.95, left=0.15)
    fig2.savefig('Variant_VLoad_Day4.svg')
    fig2.savefig('Variant_VLoad_Day4.png')

    fig6, ax6 = plt.subplots(1, 1, figsize=(3, 3))
    time = np.arange(0, 30, 0.1)
    vt = np.loadtxt('Variant_vt2.txt')[:, 0:len(time)]

    t_inds = np.arange(70, 290, 70)
    y_add=[0.1,0.4,0.5]
    x_add=[2,0,0]
    for v in range(len(vtype)):
        indices = df.index[df['vtype'] == vtype[v]].tolist()
        y = vt[indices, :]
        y[y < 1e-6] = 1e-6
        y = np.mean(np.log10(y), axis=0)
        ax6.plot(time, y, color=set2colors[v], linewidth=2)
        ax6.scatter(time[t_inds], y[t_inds], color=set2colors[v], marker=markers[v])
        ax6.text(x=time[np.argmax(y)]+x_add[v], y=np.max(y) + y_add[v], s=xlabel_simple[v], color=set2colors[v])
    ax6.annotate(xy=(4, -3.2), xytext=(4.2, -1.5), text='COVID test', color='k', ha='left',
                 arrowprops=dict(arrowstyle='->', connectionstyle="arc3", color='k', lw=1), )
    ax6.vlines(4, -4, 4, linewidth=0.5, linestyles='--', colors='k')
    ax6.hlines(0, 30, 0, linewidth=1, colors='k')
    ax6.spines['top'].set_visible(False)
    ax6.spines['right'].set_visible(False)
    ax6.set_xlim([0, 31])
    ax6.set_xticks([0, 10, 20, 30])
    ax6.set_ylim([-3.2, 4.5])
    ax6.set_yticks([-3, 0, 3])
    ax6.set_ylabel('log$_{10}$nCoV', labelpad=-5)
    ax6.set_xlabel('time (days)', labelpad=0)
    fig6.subplots_adjust(left=0.15, bottom=0.15)
    fig6.savefig('Variant_VDynamics.svg')
    fig6.savefig('Variant_VDynamics.png')

    plt.show()