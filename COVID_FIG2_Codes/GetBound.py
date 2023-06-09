#generating boundary of sampling for selected parameters with indices in sample_indices1, 2 and 3
#output txt file to be used in Heuristic_Sampling.py
import numpy as np

sample_indices1 = np.arange(0,25,1)
sample_indices1 = np.delete(sample_indices1,[8,9,11,12,13,14,15,17,18,19,21,22,23,10,20,24])
#print(sample_indices1)
sample_indices2 = np.hstack((np.array([8,9,11,12,13,14,15,17,18,19,21,22,23]),[30,34,35],np.arange(36,52,1)))#,np.arange(36,52,1)))
sample_indices3 = [150,159,160]+list(np.arange(128,132,1))+list(np.arange(27,30,1))+list(np.arange(31,34,1))

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
    else:
        LogParaBound.append([a, a])
LogParaBound = np.array(LogParaBound)
np.savetxt('logparabound.txt', LogParaBound)