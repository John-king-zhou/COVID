import numpy as np
from scipy.integrate import odeint
from scipy.optimize import curve_fit
from scipy.signal import find_peaks

def Hill(x,k,n):
    return x**n/(k**n+x**n)

def func1(x0,time,Para):
    v=x0[0]
    If=x0[1]

    ek=Para[0]
    ec=1

    dvdt=1500*0.4*If-ec*v
    dIfdt=1.2*10**(-4)*50*v-ek*If

    return np.array([dvdt,dIfdt])

def func2(x0,time,Para):
    v=x0[0]
    If=x0[1]

    if len(Para)==1:
        ec=Para[0]
        ek=3
    else:
        ec=Para[0]
        ek=Para[1]

    dvdt=1500*0.4*If-ec*v
    dIfdt=1.2*10**(-4)*50*v-ek*If

    return np.array([dvdt,dIfdt])

def Imm1(t,lgv,ek):
    time=np.arange(0,100,0.1)
    results=odeint(func1,[10**lgv,6e-3/ek*10**lgv],time,args=((ek,),))
    v=results[:,0]
    v[v<1e2]=1e2
    indices=[int(t[i]*10) for i in range(len(t))]
    return np.log10(v)[indices]

def Imm2(t,lgv,ec):
    time=np.arange(0,100,0.1)
    results=odeint(func2,[10**lgv,6e-3/3*10**lgv],time,args=((ec,),))
    v=results[:,0]
    v[v<1e2]=1e2
    indices=[int(t[i]*10) for i in range(len(t))]
    return np.log10(v)[indices]


'''
FITTING function, input: time, log viral load, whether there is a rising phase in viral dynamics (default=1),
if not, simply fit the data to an exponential decay
'''
def FIT2(time,y,rising=1):
    e1, e2 = np.nan, np.nan
    max_id=np.argmax(y)
    time = time - time[0]
    bound=[0,np.inf]
    '''
    if a peak exist: try fit e1 and e2 separately
    '''
    if 0<max_id<len(y)-1:
        '''
        calculate e1 and e2 accordingly
        '''
        peak=time[max_id]
        def Imm3(t, lgv, ek, ec):
            t0=peak
            time1 = np.arange(0, t0, 0.1)
            results1 = odeint(func1, [10 ** lgv, 6e-3 / ek * 10 ** lgv], time1, args=((ek,),))
            time2 = np.arange(t0, 100, 0.1)
            results2 = odeint(func2, results1[-1, :], time2, args=((ec,),))
            results = np.vstack((results1, results2))
            time = np.hstack((time1, time2))
            v = results[:, 0]
            v[v < 1e2] = 1e2
            indices = [int(t[i] * 10) for i in range(len(t))]
            return np.log10(v)[indices]
        popt,pcov=curve_fit(Imm3, xdata=time, ydata=y, p0=(y[0],0.01,2),maxfev=2000, method='lm')
        if not ((popt>=bound[0]) & (popt<bound[1])).all():
            print('Parameters Out of Bound')
            popt, pcov = curve_fit(Imm3, xdata=time, ydata=y, p0=(y[0], 0.01, 1), maxfev=1e3, bounds=[0, np.inf], method='trf')
        e1=popt[1]*1
        e2=popt[2]*3
        if e1<e2 and e1<=3.6:
            return True, popt, pcov, Imm3
        if e1<3.6:#rising monotonically
            max_id=len(y)-1
        else:#decaying monotonically
            max_id=0
    '''
    if fitting two phase is not possible: fitting to one phase
    '''
    if max_id==len(y)-1:
        popt, pcov = curve_fit(Imm1, xdata=time, ydata=y, p0=(y[0], 0.1), maxfev=1000,method='lm')
        if not ((popt>=bound[0]) & (popt<bound[1])).all():
            print('Parameters Out of Bound')
            popt, pcov = curve_fit(Imm1, xdata=time, ydata=y, p0=(y[0], 0.1), maxfev=1000, bounds=[0, np.inf],method='trf')
        return False, 1, popt, pcov
    else:
        popt, pcov = curve_fit(Imm2, xdata=time, ydata=y, p0=(y[0], 1), maxfev=1000,method='lm')
        if not ((popt>=bound[0]) & (popt<bound[1])).all():
            print('Parameters Out of Bound')
            popt, pcov = curve_fit(Imm2, xdata=time, ydata=y, p0=(y[0], 1), maxfev=1000, bounds=[0, np.inf],method='trf')
        return False, 2, popt, pcov
