from random import randrange
import numpy as np
#Variables:
delta_t = 0.2
t = 0
t_max = 5
k = 3.15
a = 9.2
w = 1.3
n = 3
r = np.zeros(n, len(t_set))
phi = np.zeros(n, len(t_set))
def d(i,p):
    return np.sqrt((r[p][1]-r[i][1])**2 + (r[p][1]-r[i][1])**2)
t_set = np.arange(0, t_max, delta_t)

def w(i,p):
    return (a/(e**(w*d(i,p)+a))

for p in range(n):
    for t in range(len(t_set)):
        w_i = 
        phi[p][t+1] = phi[p][t] - delta_t*k/n #sum