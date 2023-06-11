import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
from Equation import func
from Read import Data
from numba import jit
import multiprocessing
from tqdm import tqdm
from Type_Characterization2 import PeakNum
import os

# import warnings
# warnings.filterwarnings('ignore')

@jit
def get_lgVL(Para):
    '''
    To get Viral Load (VL) time courses
    unit of VL in the model: 10^6/ml
    '''

    time0 = np.arange(0, 50, 0.1)

    initial = [1e-2, 0, Para[52] / Para[62], 0, 0, Para[53]/Para[65], 0, 0, 
               Para[159], 0, 0, 0, 0, 0, 0,
               Para[54]/Para[70], 0, Para[160], 0, 0, 0, 0, 0, 0, Para[77] /
               Para[101], Para[82] / Para[102],
               Para[84]/Para[103], (Para[88]+Para[54]/Para[70]*Para[90]) /
               Para[104], Para[91]/Para[105], Para[95]/Para[106],
               0, 0]

    results = odeint(func, initial, time0, args=(Para,))

    v = results[:, 0]
    v[v < 1e-10] = 1e-10
    v[v > 1e20] = 1e20
    lgVL = np.log10(v)
    # transform unit to '1/ml' and nasal VL is 1/10 of VL in lung
    lgVL = lgVL + 5

    return lgVL, time0


@jit
def residual(para, fit_data):
    '''
    to coumpute residual for log10(VL) of each time points
    '''

    # prediction
    pred_lgVL, pred_time = get_lgVL(para)
    delta_t = pred_time[1]-pred_time[0]

    # read fit_data
    lgVL_exp = np.copy(fit_data['lgVL'])
    LOD = np.copy(fit_data['LOD'])
    time_exp = np.copy(fit_data['time'])
    is_day0 = np.copy(fit_data['is_day0'])
    time_exp_to_num = np.array(time_exp/delta_t,
                               dtype=np.int32)

    # residual

    # set exp. and simulated lgVL minimum = lg(LOD)
    lgVL_exp[lgVL_exp <= LOD] = LOD
    pred_lgVL[pred_lgVL <= LOD] = LOD

    # find minimum errors if day 0 not known
    if is_day0:
        error_by_points = (lgVL_exp-pred_lgVL[time_exp_to_num])**2
        t_init = 0
    else:
        dy2 = [np.sum((lgVL_exp-pred_lgVL[time_exp_to_num+i])**2)
               for i in range(int(10/delta_t))]
        t_init = np.argmin(dy2)*delta_t
        error_by_points = (
            lgVL_exp-pred_lgVL[time_exp_to_num+np.argmin(dy2)])**2

    # delete > 2 peaks lgVL
    if PeakNum(pred_lgVL, LOD)>1:
        error_by_points = error_by_points+10

    return error_by_points, t_init


@jit
def mean_squared_error(param, dataset):
    mse_seq = []
    for dt in dataset:
        res_tmp, _=residual(param, dt)
        mse_seq.append(res_tmp.mean())
    return mse_seq


def mean_squared_error_multiproc(x):
    return mean_squared_error(x[0], x[1])


def plot_results(data):
    fig, ax = plt.subplots(5, 6, figsize=(12, 8))
    ax = ax.flat

    for j in range(len(data)):
        dt_tmp = data[j]
        res, t_init = residual(dt_tmp['best_param'], dt_tmp)
        pred_lgVL, pred_t = get_lgVL(dt_tmp['best_param'])

        num_below_LOD = dt_tmp['lgVL'] <= dt_tmp['LOD']

        ax[j].plot(pred_t - t_init, pred_lgVL-5, c='red',)
        ax[j].plot(dt_tmp['time'][~num_below_LOD], 
                   dt_tmp['lgVL'][~num_below_LOD]-5, 
                   marker='o', markersize=4, c='k', lw=0)
        ax[j].plot(dt_tmp['time'][num_below_LOD], 
                   dt_tmp['lgVL'][num_below_LOD]-5, 
                   marker='o', markersize=4, c='k', lw=0,
                   markerfacecolor='none')
        ax[j].hlines(dt_tmp['LOD']-5, -100, 100,
                     color='gray', ls='--', lw=0.5,
                     zorder=-1)
        ax[j].text(x = 15, y = 3, 
                   s = 'RMSE=%.2f' % (np.sqrt(res.mean())))
        ax[j].set_title(dt_tmp['id'])
        ax[j].set_xlim([-1, 40])
        ax[j].set_ylim([-5, 5])
        ax[j].set_xticks([0, 10, 20, 30])

    ax[26].set_xlabel('Days post symptom onset', fontsize=18,
                      x=1.05)
    ax[12].set_ylabel(r'$log_{10}nCoV$', fontsize=18)
    fig.tight_layout()
    fig.savefig('Fig_Fit_VL.svg', dpi=200)



if __name__ == '__main__':

    try:
        fit_by_patients = np.load('VL_fit.npy', allow_pickle=True)
        print('Fit data loaded!')
    except:
        print('Start Computation!')
        # Load viral load data
        data = Data(virus='SARS-CoV-2', sheet_name='Kim2021')

        fit_by_patients = []
        for i in range(len(data.data['id'])):
            fit_by_patients.append(data.single_patient_data(i))
            fit_by_patients[-1].update({'best_param': [], 'mse': 1000, 'is_day0': False})

        # read sampling parameters
        input_args = []
        for mode in [1, 2, 3]:
            param_set=np.loadtxt('Spl_Para%s.txt' % mode)
            for i, p in enumerate(param_set):
                input_args.append([p, fit_by_patients])

        # computing errors
        mse_save = []
        pool = multiprocessing.Pool(processes=os.cpu_count()-2)
        for y, i in zip(pool.imap(mean_squared_error_multiproc, input_args),
                        tqdm(range(len(input_args)))):
            mse_save.append(y)
        mse_save = np.array(mse_save)
        pool.close()

        # sort best fit
        for i in range(len(data.data['id'])):
            best_num = np.argmin(mse_save[:, i])
            fit_by_patients[i]['best_param']=input_args[best_num][0]
            fit_by_patients[i]['mse']=mse_save[best_num, i]
        
        np.save('VL_fit.npy', fit_by_patients)

    
    plot_results(fit_by_patients)