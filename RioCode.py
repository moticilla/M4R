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
t_set = np.arange(0, t_max, delta_t)
r = np.zeros(2, n, len(t_set))
r[:][:][0] = randrange(1,3)
phi = np.zeros(n, len(t_set))
def d(i,p,t):
    return np.sqrt((r[1][p][t]-r[1][i][t])**2 + (r[2][p][t]-r[2][i][t])**2)
def w(i,p,t):
    return ( a/(np.exp(w*d(i,p,t)+a)) )

for p in range(n):
    for t in range(len(t_set)):
        for i in range(n):
            w_i = w(i,p)
            phi[p][t+1] = phi[p][t] - delta_t* k/n * w_i * np.sin(phi[i][t] - phi[p][t])
        r[1][p][t+1] = r[1][p][t] 