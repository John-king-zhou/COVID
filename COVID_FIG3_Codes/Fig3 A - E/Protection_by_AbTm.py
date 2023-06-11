#figure 3C
import matplotlib.mathtext as mathtext
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing
from Equation import func
from scipy.integrate import odeint
from E_Calculation import E, E1, E2, E12
from scipy.spatial import ConvexHull

# Settings of global variables and colors
set2colors = ['#fc8d62', '#66c2a5', '#a6d854',
              '#e78ac3', '#ffd92f', '#e5c494', '#b3b3b3']
ggcolors = ['#808080', '#1F77B4', '#FF7F0E', '#D62728',
            '#2CA02C', '#4DBEEE', '#77AC30', '#9467BD']
mathtext.FontConstantsBase.sup1 = 0.4
mode_name = ['Mode 1', 'Mode 2', 'Mode 3', 'Mode 4']
time0 = np.arange(0, 50, 0.1)


def get_results(Para, v0, CD8Tm, Ab, A, CD4Tm=0.02):

    if A != 0:
        initial = [0.01, 0, Para[52] / Para[62], 0, 0, Para[53]/Para[65], 0, 0, Para[159]/2, 0, 0, 0, 0, 0, 0,
                   Para[54]/Para[70], 0, Para[160]/2, 0, 0, 0, 0, 0, 0, Para[77] /
                   Para[101], Para[82] / Para[102],
                   Para[84]/Para[103], (Para[88]+Para[54]/Para[70]*Para[90]) /
                   Para[104], Para[91]/Para[105], Para[95]/Para[106],
                   0, 0]
        initial[31] = A
        initial[30] = Ab
        initial[23] = 0
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
    results = np.vstack((results.T, e, e1, e2+e12)).T
    return results

def get_results_merge(args):
    return get_results(*args)

def classification(vv):
    if max(vv) <= vv[0]:
        return 0
    else:
        return 1


def encircle(x, y, ax, **kw):
    p = np.c_[x, y]
    hull = ConvexHull(p)

    # find left border
    a = np.copy(p[hull.vertices, :])
    left_p = np.argmin(a[:, 0])
    down_p = np.argmin(a[:, 1])
    downleft_seq = np.array([atmp if atmp[0] <= a[down_p, 0] and atmp[1] <= a[left_p, 1] else [np.NaN, np.NaN]
                             for atmp in a])
    downleft_seq = downleft_seq[~np.isnan(downleft_seq[:, 0]), :]
    downleft_seq = downleft_seq[np.argsort(downleft_seq, axis=0)[:, 0]]
    return downleft_seq


if __name__ == '__main__':
    np.random.seed(2020)
    # select Mode and load parameter set
    mode = 3
    Parameter = np.loadtxt('GeoMean_para_Mode%i.txt' % mode)

# --------- Plot the primary infection trajectory -----------------------------
    fig, ax = plt.subplots(1, 4, figsize=(6, 2))
    colorseq = ['gray', set2colors[1], set2colors[0], ]

    v0 = 0.01
    CD8Tm = 0
    Ab = 0
    A = 0
    args = (Parameter, v0, CD8Tm, Ab, A)
    count = 0
    results = get_results(*args)
    v = results[:, 0]
    v[v < 1e-4] = 1e-4
    lgv = np.log10(v)
    ax[0].plot(time0, lgv, c=colorseq[count],
               lw=2, ls='--', zorder=100)
    ax[1].plot(time0, results[:, 26], c=colorseq[count],
               lw=2, ls='--', zorder=100)
    ax[2].plot(time0, results[:, 33], c=colorseq[count],
               lw=2, ls='--', zorder=100)
    ax[3].plot(time0, results[:, 34], c=colorseq[count],
               lw=2, ls='--', zorder=100)

# ----------- Sampling Trajectories with different initial memory cells ---------------
    try:
        data = np.load('Asingle_M%i_timecourse_spl.npy' %
                       mode, allow_pickle=True).item()
        X = data['timeseq']
        initial = np.array(data['ini'])
    except:
        print('Sampling initial Bmem, CD8+Tmem!')
        A = 1
        M = 1000
        T_Mem = np.random.uniform(0, 0.5, M)
        Ab_Mem = np.random.uniform(0,1500, M)

        cores = multiprocessing.cpu_count()
        pool = multiprocessing.Pool(processes=int(cores)-2)
        X = [[], [], [], [], []]
        initial = []
        args = []
        for j in range(M):
            args.append([Parameter, v0, T_Mem[j], Ab_Mem[j], A, 0.2])
        count = 0

        for results in pool.imap(get_results_merge, args):
            count += 1
            ratio = count / M
            rat_str = ['>'] * int(ratio * 50) + ['-'] * (50 - int(ratio * 50))
            rat_str = ''.join(rat_str)
            print('\r' + rat_str + '%.2f %%' % (ratio * 100), end='')

            v = results[:, 0]
            v[v < 1e-4] = 1e-4
            lgv = np.log10(v)
            e1 = results[:, 33]
            e2 = results[:, 34]
            IL6 = results[:, 26]
            ee = results[:, 32]

            X[0].append(lgv)
            X[1].append(IL6)
            X[2].append(ee)
            X[3].append(e1)
            X[4].append(e2)

            initial.append(results[0, :])

        # save sampling data
        X = [np.array(x) for x in X]
        pool.close()
        initial = np.array(initial)
        data = {'timeseq': X, 'ini': initial}
        np.save('Asingle_M%i_timecourse_spl.npy' % mode, data)

# ---------------- Plot the sampling trajectories ---------------------------------------
    y_pred = np.array([classification(X[0][z_i])
                       for z_i in range(X[0].shape[0])])
    Y1 = [XX[y_pred == 0] for XX in X]
    Y2 = [XX[y_pred == 1] for XX in X]
    Y = [Y1, Y2]
    type_lab = ['Protected', 'Infected']
    for i in range(2):
        vv = Y[i][0]
        il6 = Y[i][1]
        ee = Y[i][2]
        ee1 = Y[i][3]
        ee2 = Y[i][4]
        for j, seq in zip(range(6), [vv, il6, ee1, ee2]):
            ax[j].plot(time0, np.mean(seq, axis=0), color=colorseq[i + 1], lw=3)
            ax[j].fill_between(time0,
                               -np.std(seq, axis=0) + np.mean(seq, axis=0),
                               np.std(seq, axis=0) + np.mean(seq, axis=0),
                               color=colorseq[i + 1], alpha=0.5)
    for a in ax:
        a.set_xlim([0, 25])
        a.set_xticks([0, 10, 20])
    ax[0].set_ylim([-3.99, 4])
    ax[0].set_yticks([-3, 0, 3])
    ax[0].set_title(r'$log_{10}nCoV$', fontsize=12)
    ax[1].set_title('IL-6', fontsize=12)
    import matplotlib

    sFormatter1 = matplotlib.ticker.ScalarFormatter(useOffset=False, useMathText=True)
    sFormatter1.set_powerlimits((0, 2))
    ax[1].yaxis.set_major_formatter(sFormatter1)
    ax[2].set_title(r'$\epsilon_i$', fontsize=13)
    ax[2].set_yticks([0, 1])
    ax[3].set_title(r'$\epsilon_a$', fontsize=13)
    ax[3].set_yticks([0, 4, 8])

    ax[0].plot([], [], c=colorseq[0],lw=2, label='Primary Infection',ls='--')
    ax[0].plot([], [], color=colorseq[1], label=type_lab[0], lw=3)
    ax[0].plot([], [], color=colorseq[2], label=type_lab[1], lw=3)
    fig.legend(handlelength=2.5,frameon=False,ncol=3,bbox_to_anchor=(0.5,1.05),loc='upper center')
    fig.subplots_adjust(left=0.05, bottom=0.2, hspace=0.35, top=0.8, right=0.95)
    fig.text(x=0.5, y=0.05, s='time after infection t$^*$ (days)', ha='center', va='center')

# -------------------- Protection border in the Bm, CD8+Tm space ----------------------------------
    figbt, axbt = plt.subplots(figsize=(2.2, 2.2))
    Abseq = np.copy(initial[:, 30])
    bmseq = np.copy(initial[:, 23])
    cd8tmseq = np.copy(initial[:, 20])
    mkr = ['o', 'X']
    fillstl = ['full', 'none']
    fc = [set2colors[1], 'w', ]
    ec = ['w', set2colors[0], ]
    lw = [0.6,1]
    count = [0, 0, 0]
    for i in range(len(bmseq)):
        axbt.scatter(Abseq[i], cd8tmseq[i], s=40, marker=mkr[y_pred[i]],edgecolor=ec[y_pred[i]],
                  facecolor=fc[y_pred[i]], lw=lw[y_pred[i]])
    for pred in [1,0]:
        axbt.scatter([],[],s=40, marker=mkr[pred],edgecolor=ec[pred],
                  facecolor=fc[pred],label=type_lab[pred], lw=lw[pred])
    # find the borders of protected and infected.
    encir_type = list(set(y_pred))
    for i in encir_type:
        leftborder = encircle(Abseq[y_pred == i], cd8tmseq[y_pred == i],
                              axbt, ec="k", fc=colorseq[i+1], alpha=0.3, linewidth=0, zorder=-10)
        if i == 0:
            leftborder = leftborder[:-1]
            
            leftborder = np.vstack((np.array([leftborder[0,0],np.max(cd8tmseq)]),
                                leftborder,
                                np.array([leftborder[-1,0],np.min(cd8tmseq)])))
            axbt.plot(leftborder[:, 0], leftborder[:, 1],
                      lw=3, ls='-', color='k')
            axbt.fill_betweenx(leftborder[:, 1], leftborder[:, 0], leftborder[:, 0]-1000,
                               color=colorseq[i+2], alpha=0.3, zorder=-0.3)
            axbt.fill_betweenx(leftborder[:, 1], leftborder[:, 0], leftborder[:, 0]+1000,
                               color=colorseq[i+1], alpha=0.3, zorder=-0.3)

    axbt.legend(frameon=False, handletextpad=0, 
                ncol=3, columnspacing=0.1,
                labelspacing=0.1, loc='upper center',
                bbox_to_anchor = (0.38,1.23))
    axbt.set_xlim([0, 1300])
    axbt.set_ylim([0, np.max(cd8tmseq)])
    axbt.set_yticks([0.2, 0.4])
    axbt.set_xticks([0, 500, 1000])
    axbt.set_xlabel('Ab (t$^*$)')
    axbt.set_ylabel('CD8+T$_M$ (t$^*$)')
    bbox = dict(boxstyle="square", fc="w", color='k')
    axbt.text(x=1200,y=0.4,s='exp',bbox=bbox,ha='right')
    plt.tick_params(zorder=100)
    figbt.subplots_adjust(bottom=0.25,left=0.25)

    # Saving figures
    figbt.savefig("AFig_BmTm_GeoMode%i.png" % mode, dpi=300)
    fig.savefig("AFig_type_M%i.png" % mode, dpi=300)
    
    figbt.savefig("AFig_BmTm_GeoMode%i.svg" % mode, dpi=300)
    fig.savefig("AFig_type_M%i.svg" % mode, dpi=300)

    plt.show()
