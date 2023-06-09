#Latin Hypercube Sampling method for generating random parameter samples
import numpy as np
def LHSample( D, bounds, N):#D=len(Para); bounds=upper and lower bounds of the parameters 2D-array (len(Para)*2);
    #N=number of parameter sets generated
    result = np.empty([N, D])
    temp = np.empty([N])
    d = 1.0 / N

    for i in range(D):

        for j in range(N):
            temp[j] = np.random.uniform(
                low=j * d, high=(j + 1) * d, size = 1)[0]

        np.random.shuffle(temp)

        for j in range(N):
            result[j, i] = temp[j]

    b = np.array(bounds)
    lower_bounds = b[:,0]
    upper_bounds = b[:,1]
    if np.any(lower_bounds > upper_bounds):
        print('range error')
        return None

    #   sample * (upper_bound - lower_bound) + lower_bound
    np.add(np.multiply(result,
                       (upper_bounds - lower_bounds),
                       out=result),
           lower_bounds,
           out=result)
    return result


#to be called as：
#bound = np.array([[-5, 3],[-5, 3],[-5, 3],[-5, 3],[-5, 3],[-5, 3],[2, 5],[2, 5],[2, 5],[2, 5],[2, 5],[2, 5], [1.0*10**3, 1.0*10**5], [1.0*10**3, 1.0*10**5]])
#kk = LHSample(14, bound, n)