from random import random, choice, normalvariate
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import special
import matplotlib.animation as animation
#Variables:
t_max = 15
n = 31
d_max = 5
pi = np.pi
angle_max = pi/2
k = 3.15
delta_t = 0.15/k
a = 9.2
w_var = 1.3
c = 3.61
t_set    = np.arange(0, t_max, delta_t)
position = np.zeros((len(t_set), n, 2))
phi      = np.zeros((len(t_set),n))
r        = np.zeros((len(t_set),n)) #this is r dot not r
def d(i,p,t):
    return np.sqrt((position[t,p,0]-position[t,i,0])**2 + (position[t,p,1]-position[t,i,1])**2)
def angle(i,p,t):
    theta = np.arctan((position[t,i,1]-position[t,p,1])/ (position[t,i,0] - position[t,p,0]))
    return phi[t,p]-theta
def w(i,p,t):
    dis = d(i,p,t)
    return ( a/(np.exp(w_var * dis)+a) )
def ogive(mu, sigma,t):
    return 0.5 *(1+special.erf((t-mu)/sigma*np.sqrt(2)))
def perturb(S, t, sign = 1, Heading = True, Speed = False):
    #S is the set of virtual neighbours to be perturbed either by speed or by heading at time t
    for i in S:
        if Heading == True:
            for all_t in range(int(t//delta_t),len(t_set)):
                phi[all_t, i] = sign*(10/360)*2*pi
            for t_step in range(int(0.5//delta_t)):
                phi[int(t//delta_t)+t_step,i] = sign*(10/360)*2*pi*ogive(t,0.083,t_step*delta_t-0.25)
        if Speed   == True:
            for all_t in range(int(t//delta_t),len(t_set)):
                r[all_t, i] = r[int(t//delta_t),i]   + sign*0.3
            for t_step in range(int(0.5//delta_t)):
                r[int(t//delta_t)+t_step,i]   = r[int(t//delta_t),i]   + sign*0.3*ogive(t,0.083,t_step*delta_t-0.25)

#set position of 30 neighbours for experiment 1 in two circles with real person in the middle
for i in range(14):
    position[0,i,0] = normalvariate(1.5,0.15)*np.sin(normalvariate(i*pi/7,(8/360)*2*pi))
    position[0,i,1] = normalvariate(1.5,0.15)*np.cos(normalvariate(i*pi/7,(8/360)*2*pi))
for i in range(14,30):
    position[0,i,0] = normalvariate(3.5,0.15)*np.sin(normalvariate((i-14)*pi/8,(8/360)*2*pi))
    position[0,i,1] = normalvariate(3.5,0.15)*np.cos(normalvariate((i-14)*pi/8,(8/360)*2*pi))
position[0,30,0] = 0
position[0,30,1] = 0
plt.scatter(position[0,:,0], position[0,:,1])
for i in range(31):
    plt.annotate(angle(i,30,0)*180/pi, (position[0,i,0], position[0,i,1]))
    #print(angle(i,30,0)*180/pi)
plt.show()

def move_person():
    total_distance = 0
    t = 1
    global phi
    global r
    global position
    while total_distance<12 and t<len(t_set):
        p = 30
        r[t,p]   = r[t-1,p]
        phi[t,p] = phi[t-1,p]
        localN   = 0
        phi_sum  = 0
        r_sum    = 0
        for i in range(30):
            dis = d(i,p,t-1)
            ang = angle(i,p,t-1)
            if (dis< d_max) and (np.absolute(ang)<angle_max) and (i != p):
                localN  += 1
                w_i      = w(i,p,t-1)
                phi_sum += w_i * np.sin(phi[t-1,i] - phi[t-1,p])
                r_sum   += w_i *(r[t-1,i]-r[t-1,p])
        if localN > 0:
            r[t,p]   = r[t-1,p]   + delta_t * (c/localN) * r_sum
            phi[t,p] = phi[t-1,p] + delta_t * (k/localN) * phi_sum #need to decide if + or -
        distance        = r * delta_t
        total_distance= np.sum(distance[:,30])
        position[t,p,0] = position[t-1,p,0] + distance[t-1,p]*np.cos(phi[t-1,p])
        position[t,p,1] = position[t-1,p,1] + distance[t-1,p]*np.sin(phi[t-1,p])
        t += 1

#make all virtual people speed up to 1.3m/s in 3 s and then keep going
for i in range(30):
    for t in range(int(3//delta_t),len(t_set)):
        r[t,i] = 1.3
    for t in range(int(3//delta_t)):
        r[t,i] = 1.3*ogive(0,0.5,t*delta_t - 1.5)

my_exp1_final_heading_data = [[1,2,3,0,0,0,0,0,0,0],[1,2,3,0,0,0,0,0,0,0]]
def all_together(i,S, h=True, speed =False):
    global position
    global phi
    global r
    #make the perturbed set change heading up after total of 5 seconds
    perturb(S, 5,Heading = h, Speed = speed)
    #get position of virtual neighbours from the perturbed speed and heading
    for t in range(1, len(t_set)):
        for p in range(n):
                distance        = r * delta_t
                position[t,p,0] = position[t-1,p,0] + distance[t-1,p]*np.cos(phi[t-1,p])
                position[t,p,1] = position[t-1,p,1] + distance[t-1,p]*np.sin(phi[t-1,p])
    #move the "real" person (for 12m)
    move_person()
    #find the t that corresponds with stopping after 12m
    final_t = 3
    while r[final_t,30] >0 and final_t<1000:
        final_t+=1
    #update dataframe
    if h == True:
        my_exp1_final_heading_data[0][i]=phi[final_t-1,30]*360/(2*pi)
    if speed ==True:
        my_exp1_final_heading_data[1][i]=r[final_t-1,30]
    #un perturb virtual people by speed
    for i in range(30):
        for t in range(int(3//delta_t),len(t_set)):
            r[t,i] = 1.3
        for t in range(int(3//delta_t)):
            r[t,i] = 1.3*ogive(0,0.5,t*delta_t - 1.5)
    for i in range(len(t_set)):
        r[i,30]=0
    #un perturb by heading
    for i in range(30):
        phi[:,i] = 0

def animate_everyone():
    #animate the people moving over time:
    fig, ax = plt.subplots()
    scat1 = ax.scatter(position[0,0,0],position[0,0,1], c="y", s=5, label='virtual person')
    scat = [scat1]
    for i in range(0,8):
        scat1 = ax.scatter(position[0,i,0],position[0,i,1], c="y", s=5)
        scat.append(scat1)
    for i in range(8,15):
        scat1 = ax.scatter(position[0,i,0],position[0,i,1], c="b", s=5)
        scat.append(scat1)
    for i in range(15,23):
        scat1 = ax.scatter(position[0,i,0],position[0,i,1], c="y", s=5)
        scat.append(scat1)
    for i in range(23,n-1):
        scat1 = ax.scatter(position[0,i,0],position[0,i,1], c="b", s=5)
        scat.append(scat1)
    scat30 = ax.scatter(position[0,30,0],position[0,30,1], c="r", s=5, label='"real" person')
    scat.append(scat30)
    ax.set(xlim=[-5, 26], ylim=[-4, 4], xlabel='x (m)', ylabel='y (m)')
    ax.legend()

    def update(frame):
        for i in range(n-1):
            # for each frame, update the data stored on each artist.
            x1 = position[:frame,i,0]
            y1 = position[:frame,i,1]
            # update the scatter plot:
            data = np.stack([x1, y1]).T
            scat[i].set_offsets(data)
        x1 = position[:frame,30,0]
        y1 = position[:frame,30,1]
        # update the scatter plot:
        data = np.stack([x1, y1]).T
        scat30.set_offsets(data)
        return (scat30, (scat[i] for i in range(n-1)))
    ani = animation.FuncAnimation(fig=fig, func=update, frames=300, interval=12)
    plt.show()

#define subset S for experiment 1 Near 0,Far 0, heading
S = np.array([])
all_together(0,S)
all_together(1,S)
#define subset S for experiment 1 Near 3, heading
S = np.array([0,1,2])
all_together(2,S)
#define subset S for experiment 1 Near 6, heading
S = np.array([0,1,2,3,4,5])
all_together(3,S)
#animate_everyone()
#define subset S for experiment 1 Near 9, heading
S = np.array([0,1,2,3,4,5,6,7,8])
all_together(4,S)
#define subset S for experiment 1 Near 12, heading
S = np.array([0,1,2,3,4,5,6,7,8,9,10,11])
all_together(5,S)
#define subset S for experiment 1 Far 3, heading
S = np.array([14,15,16])
all_together(6,S)
#define subset S for experiment 1 Far 6, heading
S = np.array([14,15,16,17,18,19])
all_together(7,S)

#define subset S for experiment 1 Far 9, heading
S = np.array([14,15,16,17,18,19,20,21,22])
all_together(8,S)
animate_everyone()
#define subset S for experiment 1 Far 12, heading
S = np.array([14,15,16,17,18,19,20,21,22,23,24,25,26])
all_together(9,S)
#define subset S for experiment 1 Near 0,Far 0, speed
S = np.array([])
all_together(0,S,h=False, speed =True)
all_together(1,S,h=False, speed =True)
#define subset S for experiment 1 Near 3, speed
S = np.array([0,1,2])
all_together(2,S,h=False, speed =True)
#define subset S for experiment 1 Near 6, speed
S = np.array([0,1,2,3,4,5])
all_together(3,S,h=False, speed =True)
#define subset S for experiment 1 Near 9, speed
S = np.array([0,1,2,3,4,5,6,7,8])
all_together(4,S,h=False, speed =True)
#define subset S for experiment 1 Near 12, speed
S = np.array([0,1,2,3,4,5,6,7,8,9,10,11])
all_together(5,S,h=False, speed =True)
#define subset S for experiment 1 Far 3, speed
S = np.array([14,15,16])
all_together(6,S,h=False, speed =True)
#define subset S for experiment 1 Far 6, speed
S = np.array([14,15,16,17,18,19])
all_together(7,S,h=False, speed =True)
#define subset S for experiment 1 Far 9, speed
S = np.array([14,15,16,17,18,19,20,21,22])
all_together(8,S,h=False, speed =True)
#define subset S for experiment 1 Far 12, speed
S = np.array([14,15,16,17,18,19,20,21,22,23,24,25,26])
all_together(9,S,h=False, speed =True)


#make dataFrame
my_exp1_final_heading = pd.DataFrame(my_exp1_final_heading_data, columns=['Near0', 'Far0','Near3','Near6','Near9',
                                                                          'Near12','Far3','Far6','Far9','Far12'])
#for i in range(n):
#    plt.scatter(position[:,i,0],position[:,i,1], s=5)
#plt.show()

#my_quiver = np.zeros((len(t_set), n, 2))
#for t in range(len(t_set)-1):
#    for p in range(n):
#        my_quiver[t,p,0] = -position[t,p,0] + position[t+1,p,0]
#        my_quiver[t,p,1] = -position[t,p,1] + position[t+1,p,1]
#my_quiver = my_quiver/200
#for i in range(n):
#    plt.quiver(position[::4,i,0], position[::4,i,1],my_quiver[::4,i,0], my_quiver[::4,i,1], units = 'width')

#plt.show()

#read in data from the paper
exp1_LateralDeviation_strip = pd.read_csv('/home/ellie/Documents/M4R/exp1/Exp1_LateralDeviation_strip.txt', 
                                          sep="\s+", skip_blank_lines=True, 
                                          header=None, skiprows = [0,1,12,15,16,17,18, 29,32,33,34,35,46])
exp1_LateralDeviation_strip.columns = ["Subject", "Near 0", "Far 0", "Near 3", "Near 6", "Near 9", "Near 12", 
                                       "Far 3", "Far 6", "Far 9", "Far 12"]
exp1_LateralDeviation_strip.drop([10,13,14,25,28,29,40], axis = 0)


exp1_final_heading = pd.read_csv('/home/ellie/Documents/M4R/exp1/Simulation_Exp1_FinalHeading.txt', 
                                          sep="\s+", skip_blank_lines=True, 
                                          header=None, skiprows = [0,1,12,15,17])
exp1_final_heading.columns = ["Subject", "Near0", "Far0", "Near3", "Near6", "Near9", "Near12", 
                                       "Far3", "Far6", "Far9", "Far12"]

exp1_final_speed = pd.read_csv('/home/ellie/Documents/M4R/exp1/Simulation_Exp1_FinalSpeed.txt', 
                                          sep="\s+", skip_blank_lines=True, 
                                          header=None, skiprows = [0,1,12,15,17])
exp1_final_speed.columns = ["Subject", "Near0", "Far0", "Near3", "Near6", "Near9", "Near12", 
                                       "Far3", "Far6", "Far9", "Far12"]

x = [0,3,6,9,12]
y_near_speed = [exp1_final_speed.Near0[11],exp1_final_speed.Near3[11],exp1_final_speed.Near6[11],exp1_final_speed.Near9[11],exp1_final_speed.Near12[11]]
y_far_speed = [exp1_final_speed.Far0[11],exp1_final_speed.Far3[11],exp1_final_speed.Far6[11],exp1_final_speed.Far9[11],exp1_final_speed.Far12[11]]
y_near_heading = [exp1_final_heading.Near0[11],exp1_final_heading.Near3[11],exp1_final_heading.Near6[11],exp1_final_heading.Near9[11],exp1_final_heading.Near12[11]]
y_far_heading = [exp1_final_heading.Far0[11],exp1_final_heading.Far3[11],exp1_final_heading.Far6[11],exp1_final_heading.Far9[11],exp1_final_heading.Far12[11]]

plt.subplot(2, 1, 1)
plt.plot(x,y_near_heading,color = 'b',label = 'Paper near')
plt.plot(x,y_far_heading,color = 'r', label = 'paper far')
plt.title(f'Final headings and speeds from paper, k={k}, a={a}')
plt.ylabel('Final lateral deviation')

plt.subplot(2, 1, 2)
plt.plot(x,y_near_speed,color = 'b', label = 'Paper near')
plt.plot(x,y_far_speed,color = 'r', label = 'Paper far')
plt.ylabel('Final change in speed')
plt.xlabel('Number of perturbed neighours')

#plt.show()

y_near_heading = [my_exp1_final_heading_data[0][i] for i in [0,2,3,4,5]]
y_far_heading = [my_exp1_final_heading_data[0][i] for i in [1,6,7,8,9]]
y_near_speed = [my_exp1_final_heading_data[1][i] for i in [0,2,3,4,5]]
y_far_speed = [my_exp1_final_heading_data[1][i] for i in [1,6,7,8,9]]

plt.subplot(2, 1, 1)
plt.plot(x,y_near_heading,color = 'g',label = 'My near')
plt.plot(x,y_far_heading,color = 'y',label = 'My far')
plt.title(f'My Final headings and speeds, k={k}, a={a}')
plt.ylabel('Final lateral deviation')

plt.subplot(2, 1, 2)
plt.plot(x,y_near_speed,color = 'g', label = 'My near')
plt.plot(x,y_far_speed,color = 'y', label = 'My far')
plt.ylabel('Final change in speed')
plt.xlabel('Number of perturbed neighours')

plt.show()

