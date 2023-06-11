import numpy as np
import multiprocessing as mp
from sklearn.utils import resample

from _functions import *


# -------- Bootstrapping for prediction uncertainty -----------------------
class boostrap_method(object):
    def __init__(self, data_original, eth_fit) -> None:
        self.data = data_original
        self.eth = eth_fit
        pass

    def single_sample(self, j):
        dt_tmp = dict.copy(self.data)
        for i in range(len(dt_tmp['name'])):
            xtmp = np.copy(dt_tmp['normAb'][i])
            ytmp = np.copy(dt_tmp['normT'][i])
            xtmp = xtmp[~np.isnan(xtmp)]
            ytmp = ytmp[~np.isnan(ytmp)]
            dt_tmp['normAb'][i] = resample(xtmp, replace=True)
            dt_tmp['normT'][i] = resample(ytmp, replace=True)
        mu_tmp, sigma_tmp, dt_tmp = fit_distribution(dt_tmp, 
                                                     is_figure=False,
                                                     is_disp=False)
        
        return [efficacy(mu_tmp[i], sigma_tmp[i], self.eth) for i in range(len(mu_tmp))]

    def boostrap(self, N_boostrap = 10):
        eff_boostrap = []
        args = [j for j in range(N_boostrap)]
        p = mp.Pool(processes=mp.cpu_count()-2)
        for y, i in zip(p.imap(self.single_sample, args),
                        args):
            print('\r Boostrapping: [%.1f%%]' %(100*i/N_boostrap),
                  flush=True, end='')
            eff_boostrap.append(y)
        p.close()
        # for i in args:
        #     eff_boostrap.append(single_sample(i))
        #     print(i)
        return np.array(eff_boostrap)