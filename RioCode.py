from random import randrange
import numpy as np
#Variables:
delta_t = 0.2
t_max = 5
n = 3
d_max = 4
angle_max = 90
k = 3.15
a = 9.2
w_var = 1.3
c = 3.61
t_set = np.arange(0, t_max, delta_t)
r = np.zeros((len(t_set), n, 2))
r[:][:][0] = randrange(1,3)
phi = np.zeros((len(t_set),n))
def d(i,p,t):
    return np.sqrt((r[t][p][0]-r[t][i][0])**2 + (r[t][p][1]-r[t][i][1])**2)
def angle(i,p,t):
    return phi[t][p]-phi[t][i]
def w(i,p,t):
    dis = d(i,p,t)
    return ( a/(np.exp(w_var * dis + a)) )

myArray = np.array([[[1,2],[1,2],[1,2]],[[1,2],[1,2],[1,2]],[[1,2],[1,2],[1,2]],[[1,2],[1,2],[1,2]]])
#print(np.shape(myArray))

for t in range(len(t_set)-1):
    for p in range(n):
        for i in range(n):
            dis = d(i,p,t)
            ang = angle(i,p,t)
            if (dis< d_max) and (np.absolute(ang)<angle_max):
                w_i = w(i,p,t)
                phi[t+1][p] = phi[t][p] - delta_t* (k/n) * w_i * np.sin(phi[t][i] - phi[t][p])