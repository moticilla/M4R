from random import randrange
import numpy as np
import matplotlib.pyplot as plt
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
position = np.zeros((len(t_set), n, 2))
for i in range(n):
    position[0,i,0] = randrange(1,3)
    position[0,i,1] = randrange(1,3)
phi = np.zeros((len(t_set),n))
r = np.zeros((len(t_set),n))
def d(i,p,t):
    return np.sqrt((position[t,p,0]-position[t,i,0])**2 + (position[t,p,1]-position[t,i,1])**2)
def angle(i,p,t):
    return phi[t,p]-phi[t,i]
def w(i,p,t):
    dis = d(i,p,t)
    return ( a/(np.exp(w_var * dis + a)) )

#print(r)

for t in range(1, len(t_set)):
    for p in range(n):
        r[t,p]   = r[t-1,p]
        phi[t,p] = phi[t-1,p]
        for i in range(n):
            localN = 0
            phi_sum = 0
            r_sum = 0
            dis = d(i,p,t-1)
            ang = angle(i,p,t-1)
            if (dis< d_max) and (np.absolute(ang)<angle_max) and (i != p):#the angle condition is wrong
                localN += 1
                w_i = w(i,p,t-1)
                phi_sum += w_i * np.sin(phi[t-1,i] - phi[t-1,p])
                r_sum           += w_i *(r[t-1,i]-r[t-1,p])
        if localN >0:
            r[t,p]   = r[t-1,p] + delta_t * (c/localN) * r_sum
            phi[t,p] = phi[t-1,p] - delta_t * (k/localN) * phi_sum
        distance = r * delta_t
        position[t,p,0] = position[t-1,p,0] + distance[t-1,p]*np.cos(phi[t-1,p])
        position[t,p,1] = position[t-1,p,1] + distance[t-1,p]*np.sin(phi[t-1,p])

#plt.plot(position[:,1,0],position[:,1,1] )
#plt.show()
print(position[:,0,:])
print(position[:,1,:])
print(position[:,2,:])