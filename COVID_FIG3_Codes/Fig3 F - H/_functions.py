import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.stats.mstats import gmean
from scipy.optimize import minimize
import seaborn as sns
from scipy.integrate import quad
from scipy.optimize import curve_fit

set2colors = ['#fc8d62', '#66c2a5', '#a6d854',
              '#e78ac3', '#ffd92f', '#e5c494', '#b3b3b3']

# ---------------- Fitting Distribution ------------------
def negative_log_likelihood_func(theta, x, limit):
    mu, sigma = theta
    lkhfunc = 1
    for val in x:
        if val < limit:
            lkhfunc = lkhfunc*norm.cdf(val, loc=mu, scale=sigma)
        else:
            lkhfunc = lkhfunc*norm.pdf(val, loc=mu, scale=sigma)
    return -np.log10(lkhfunc)


def fit_distribution(dt, c1=0, c2=0, is_figure=True, is_disp=True):
    '''
    Fit cellular (T) and humoral (antibody) data with a logarithmic two-dimensional gaussian distribution
    by minimizing negative log likelihood function.

    Parameters
    ----------
    dt: dict
        individual-level data including normalized T cells, normalized antibody titers,
        and detection limits of T cell and antibody data.
    c1, c2: optionals
        basic levels of T cells and antibody titers.
    is_figure: bool, optional
        if True, plot the histogram of raw data and the fitted distribution.

    Returns
    -------
    mu: list
        the mean values.
    sigma: list
        the standard deviations.
    '''
    limit = dt['limit']

    if is_figure:
        figdis, axdis = plt.subplots(len(limit), 2)
        axdis[0, 0].set_title('Normalized Antibody titre')
        axdis[0, 1].set_title("Normalized T response")

    mu = []
    sigma = []
    for i in range(len(dt['name'])):
        x = dt['normAb'][i]+c2
        x = x[~np.isnan(x)]
        xlimit = limit[i][0]
        x[x <= 0] = xlimit
        logx = np.log10(x)
        resx = minimize(negative_log_likelihood_func, [[np.mean(logx), 1]], method='Nelder-Mead',
                        args=(logx, np.log10(xlimit)),
                        options={'disp': is_disp})

        y = dt['normT'][i]+c1
        y = y[~np.isnan(y)]
        ylimit = limit[i][1]
        y[y <= 0] = ylimit
        logy = np.log10(y)
        resy = minimize(negative_log_likelihood_func, [[0, 0.5]], method='Nelder-Mead',
                        args=(logy, np.log10(ylimit)),
                        options={'disp': is_disp})

        mu.append([resx.x[0], resy.x[0]])
        sigma.append([resx.x[1], resy.x[1]])
        dt['logAb'].append(logx)
        dt['logT'].append(logy)

        if is_figure:
            axdis[i, 0].hist(logx, density=True)
            axdis[i, 0].plot(np.linspace(-2, 4, 100),
                            norm.pdf(np.linspace(-2, 4, 100), mu[i][0], sigma[i][0]),
                            c=set2colors[i])
            axdis[i, 1].hist(logy, density=True)
            axdis[i, 1].plot(np.linspace(-2, 4, 100),
                            norm.pdf(np.linspace(-2, 4, 100), mu[i][1], sigma[i][1]),
                            c=set2colors[i])
            axdis[i, 0].set_ylabel(dt['name'][i])

        if is_disp and len(logx) == len(logy):
            print(dt['name'][i]+' Pearson=%.2f' %
                np.corrcoef(np.array([logx, logy]))[0, 1])
    
    return mu, sigma, dt
    

# --------- Fitting Efficacy -----------------------------------------------------

def f_eff(y, muy, sy, mux, sx, eth):
    return norm.pdf(y, loc=muy, scale=sy)*(1-norm.cdf(eth-y, loc=mux, scale=sx))


def efficacy(mu, sigma, eth):
    mux, muy = mu
    sx, sy = sigma
    return quad(f_eff, -np.inf, np.inf, args=(muy, sy, mux, sx, eth))[0]


def fit_eff(xdata, eth):
    mux = xdata[:, 0]
    muy = xdata[:, 1]
    sx = xdata[:, 2]
    sy = xdata[:, 3]
    E = np.array([quad(f_eff, -np.inf, np.inf, args=(muy[j], sy[j], mux[j], sx[j], eth))[0]
                  for j in range(len(mux))])
    return E


def root_mean_square_error(y, yp):
    return np.sqrt(np.mean((y-yp)**2))


def fit_e_threshold(VE, mu, sigma, e0=0.1):
    '''
    Fit the immune efficacy threshold by immune distribution of T cells and antibody titers and vaccine efficacy.

    VE = Normalized Size[log(T)+log(Ab) > immune efficacy threshold]

    Parameters
    ----------
    VE: array_like
        vaccie efficacy reported in clinical trials.
    mu, sigma: list
        logarithmic mean values and standard deviations of vaccine-induced T cells and antibody titers.

    Returns
    -------
    eth_fit: array
        estimated immune efficacy threshold
    '''

    # fitting data
    x_for_fit = np.array(
        [np.array([mu[i][0], mu[i][1], sigma[i][0], sigma[i][1]]) for i in range(len(mu))])
    eff_phase3_for_fit = np.copy(VE)

    eth_fit, _ = curve_fit(fit_eff, xdata=x_for_fit,
                        ydata=eff_phase3_for_fit, p0=e0)

    # compute fitted RMSE
    predict_eff = []
    for i in range(len(mu)):
        predict_eff.append(efficacy(mu[i], sigma[i], eth_fit))
    rmse = root_mean_square_error(VE, np.array(predict_eff))
    print('Fitted threshold=%.3f, with RMSE=%.2f%%' % (eth_fit, rmse*100))

    return eth_fit


# -------- Draw fitting results -------------------------------------------
def two_norm_sample(mu, sigma, size):
    cov = [[sigma[0]**2, 0], [0, sigma[1]**2]]
    x, y = np.random.multivariate_normal(mu, cov, size).T
    return (x, y)
