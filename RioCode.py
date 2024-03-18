from random import random, choice, normalvariate
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import special
import matplotlib.animation as animation
import time
from itertools import combinations

start = time.time()
#Variables:
t_max = 15
n = 31
d_max = 5
pi = np.pi
angle_max = pi/2
k = 100
delta_t = 0.6/k
a = 9.2
w_var = 1.3
c = 250
t_set    = np.arange(0, t_max, delta_t)
position = np.zeros((len(t_set), n, 2))
phi      = np.zeros((len(t_set),n))
r        = np.zeros((len(t_set),n)) #this is r dot not r
copy_phi = phi
copy_r = r
copy_position = position
stage =0
def d(pos,position,p,t):
    return np.sqrt((position[t,p,0]-pos[t,:,0])**2 + (position[t,p,1]-pos[t,:,1])**2)
def angle(pos,position,p,t):
    theta = np.arctan(np.absolute(pos[t,:,1]-position[t,p,1])/ (pos[t,:,0] - position[t,p,0]))
    return (np.absolute(phi[t,p]-theta)-pi/2)*np.sign(theta) + pi/2
def w(pos,p,t):
    dis = d(pos,position,p,t)
    return ( a/(np.exp(w_var * dis)+a) )
def ogive(mu, sigma,t):
    return 0.5 *(1+special.erf((t-mu)/sigma*np.sqrt(2)))
def perturb(S, t, phi,r,sign = 1, Heading = True, Speed = False):
    #S is the set of virtual neighbours to be perturbed either by speed or by heading at time t
    for i in S:
        if Heading == True:
            #for all_t in range(int(t//delta_t),len(t_set)):
            phi[int(t//delta_t):, i] = sign*(10/360)*2*pi
            for t_step in range(int(0.5//delta_t)):
                phi[int(t//delta_t)+t_step,i] = sign*(10/360)*2*pi*ogive(t,0.083,t_step*delta_t-0.25)
        if Speed   == True:
            #for all_t in range(int(t//delta_t),len(t_set)):
            r[int(t//delta_t):, i] += sign*0.3
            for t_step in range(int(0.5//delta_t)):
                r[int(t//delta_t)+t_step,i]   = r[int(t//delta_t),i]   + sign*0.3*ogive(t,0.083,t_step*delta_t-0.25)
    return (phi,r)

#set position of 30 neighbours for experiment 1 in two circles with real person in the middle
for i in range(14):
    position[0,i,0] = normalvariate(1.5,0.15)*np.cos(normalvariate(i*pi/7,(8/360)*2*pi))
    position[0,i,1] = normalvariate(1.5,0.15)*np.sin(normalvariate(i*pi/7,(8/360)*2*pi))
for i in range(14,30):
    position[0,i,0] = normalvariate(3.5,0.15)*np.cos(normalvariate((i-14)*pi/8,(8/360)*2*pi))
    position[0,i,1] = normalvariate(3.5,0.15)*np.sin(normalvariate((i-14)*pi/8,(8/360)*2*pi))
position[0,30,0] = 0
position[0,30,1] = 0
'''plt.scatter(position[0,:,0], position[0,:,1])
for i in range(31):
    plt.annotate((i,angle(position,position,30,0)[i]*180/pi), (position[0,i,0], position[0,i,1]))
plt.show()'''

def move_person(phi,r,position):
    total_distance = 0
    t = 1
    p = 30
    while total_distance<12 and t<len(t_set):
        r[t,p]   = r[t-1,p]
        phi[t,p] = phi[t-1,p]
        localN   = 0
        phi_sum  = 0
        r_sum    = 0
        reduced_position = position[:,:30,:]
        dis = d(reduced_position,position,p,t-1)
        ang = angle(reduced_position,position,p,t-1)
        w_i = w(reduced_position,p,t-1)
        for i in range(30):
            if (dis[i]< d_max) and (np.absolute(ang[i])<angle_max) and (i != p):
                localN  += 1
                phi_sum += w_i[i] * np.sin(phi[t-1,i] - phi[t-1,p])
                r_sum   += w_i[i] *(r[t-1,i]-r[t-1,p])
        if localN > 0:
            r[t,p]   = r[t-1,p]   + delta_t * (c/localN) * r_sum
            phi[t,p] = phi[t-1,p] + delta_t * (k/localN) * phi_sum #need to decide if + or -
        distance        = r * delta_t
        total_distance= np.sum(distance[:,30])
        position[t,p,0] = position[t-1,p,0] + distance[t-1,p]*np.cos(phi[t-1,p])
        position[t,p,1] = position[t-1,p,1] + distance[t-1,p]*np.sin(phi[t-1,p])
        t += 1
    return phi,r,position

#make all virtual people speed up to 1.3m/s in 3 s and then keep going
r[int(3//delta_t):len(t_set),:30] = 1.3
for t in range(int(3//delta_t)):
    r[t,:30] = 1.3*ogive(0,0.5,t*delta_t - 1.5)

my_exp1_final_heading_data          = np.array([[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0]])
my_exp1_final_heading_data_opposite = np.array([[0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0]])

def all_together(i,S,position, phi, r, copy_position, copy_phi, copy_r, 
                 data = my_exp1_final_heading_data, h=True, s =False, sign2 = 1):
    global stage
    stage +=1
    print(stage)
    #make the perturbed set change heading up after total of 5 seconds
    phi,r = perturb(S, 5, phi, r, Heading = h, Speed = s, sign = sign2)
    #get position of virtual neighbours from the perturbed speed and heading
    for t in range(1, len(t_set)):
        #for p in range(n):
            distance        = r * delta_t
            position[t,:30,0] = position[t-1,:30,0] + distance[t-1,:30]*np.cos(phi[t-1,:30])
            position[t,:30,1] = position[t-1,:30,1] + distance[t-1,:30]*np.sin(phi[t-1,:30])
    #move the "real" person (for 12m)
    phi, r, position = move_person(phi,r,position)
    #find the t that corresponds with stopping after 12m
    final_t = 3
    while r[final_t,30] >0 and final_t<1000:
        final_t+=1
    #update dataframe
    if h == True:
        data[0][i]+=phi[final_t-1,30]*360/(2*pi)
    if s ==True:
        data[1][i]+=r[final_t-1,30]
    copy_phi = phi.copy()
    copy_r = r.copy()
    copy_position = position.copy()
    #un perturb virtual people by speed
    #for i in range(30):
    #for t in range(int(3//delta_t),len(t_set)):
    r[int(3//delta_t):,:] = 1.3
    for t in range(int(3//delta_t)):
        r[t,:] = 1.3*ogive(0,0.5,t*delta_t - 1.5)
    #for t in range(len(t_set)):
    r[:,30]=0
    #un perturb by heading
    #for i in range(30):
    phi[:,:30] = 0
    return(position, phi, r, copy_position, copy_phi, copy_r, data)

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

    plt.annotate(i, (position[0,i,0], position[0,i,1]))
    ani = animation.FuncAnimation(fig=fig, func=update, frames=300, interval=0.2/delta_t)
    plt.show()

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

y_near_speed_paper = [exp1_final_speed.Near0[11],exp1_final_speed.Near3[11],exp1_final_speed.Near6[11],exp1_final_speed.Near9[11],exp1_final_speed.Near12[11]]
y_far_speed_paper = [exp1_final_speed.Far0[11],exp1_final_speed.Far3[11],exp1_final_speed.Far6[11],exp1_final_speed.Far9[11],exp1_final_speed.Far12[11]]
y_near_heading_paper = [exp1_final_heading.Near0[11],exp1_final_heading.Near3[11],exp1_final_heading.Near6[11],exp1_final_heading.Near9[11],exp1_final_heading.Near12[11]]
y_far_heading_paper = [exp1_final_heading.Far0[11],exp1_final_heading.Far3[11],exp1_final_heading.Far6[11],exp1_final_heading.Far9[11],exp1_final_heading.Far12[11]]

def create_final_speed_heading(position, phi, r, copy_position, copy_phi, copy_r, 
                 data3=my_exp1_final_heading_data, sign3=1):
    #define subset S for experiment 1 Near 0,Far 0, heading
    data3[:][:] = 0
    S = np.array([])
    position, phi, r, copy_position, copy_phi, copy_r, data3 = all_together(0,S,position, phi, r, copy_position, copy_phi, copy_r, data = data3, sign2 =sign3)
    position, phi, r, copy_position, copy_phi, copy_r, data3 = all_together(1,S,position, phi, r, copy_position, copy_phi, copy_r, data = data3, sign2 =sign3)
    #define subset S for experiment 1 Near 3, heading
    #S = np.random.choice([0,1,2,12,13],3, replace=False)
    for S in combinations([0,1,2,12,13], 3):
        position, phi, r, copy_position, copy_phi, copy_r, data3 = all_together(2,S,position, phi, r, copy_position, copy_phi, copy_r, data = data3, sign2 =sign3)
    data3[0][2]=data3[0][2]/10
    #define subset S for experiment 1 Near 6, heading
    #S = np.concatenate((np.array([0,1,2,12,13]),np.random.choice([14,15,16,17,27,28,29],2, replace=False)))
    for i in combinations([14,15,16,17,27,28,29],1):
        S = np.concatenate((np.array([0,1,2,12,13]), np.array(i)))
        position, phi, r, copy_position, copy_phi, copy_r, data3 = all_together(3,S,position, phi, r, copy_position, copy_phi, copy_r, data = data3, sign2 =sign3)
    data3[0][3]=data3[0][3]/7
    #define subset S for experiment 1 Near 9, heading
    #S = np.concatenate((np.array([0,1,2,12,13]),np.random.choice([14,15,16,17,27,28,29],4, replace=False)))
    for i in combinations([14,15,16,17,27,28,29],4):
        S = np.concatenate((np.array([0,1,2,12,13]), np.array(i)))
        position, phi, r, copy_position, copy_phi, copy_r, data3 = all_together(4,S,position, phi, r, copy_position, copy_phi, copy_r, data = data3, sign2 =sign3)
    data3[0][4]=data3[0][4]/35
    #define subset S for experiment 1 Near 12, heading
    S = np.concatenate((np.array([0,1,2,12,13]),np.array([14,15,16,17,27,28,29])))
    position, phi, r, copy_position, copy_phi, copy_r, data3 = all_together(5,S,position, phi, r, copy_position, copy_phi, copy_r, data = data3, sign2 =sign3)
    #define subset S for experiment 1 Far 3, heading
    #S = np.random.choice([14,15,16,17,27,28,29],3, replace=False)
    for S in combinations([14,15,16,17,27,28,29],3):
        position, phi, r, copy_position, copy_phi, copy_r, data3 = all_together(6,S,position, phi, r, copy_position, copy_phi, copy_r, data = data3, sign2 =sign3)
    data3[0][6]=data3[0][6]/35
    #define subset S for experiment 1 Far 6, heading
    #S = np.random.choice([14,15,16,17,27,28,29],6, replace=False)
    for S in combinations([14,15,16,17,27,28,29],6):
        position, phi, r, copy_position, copy_phi, copy_r, data3 = all_together(7,S,position, phi, r, copy_position, copy_phi, copy_r, data = data3, sign2 =sign3)
    data3[0][7]=data3[0][7]/7
    #define subset S for experiment 1 Far 9, heading
    #S = np.concatenate((np.array([14,15,16,17,27,28,29]),np.random.choice([0,1,2,12,13],2, replace=False)))
    for i in combinations([0,1,2,12,13],2):
        S = np.concatenate((np.array([14,15,16,17,27,28,29]), np.array(i)))
        position, phi, r, copy_position, copy_phi, copy_r, data3 = all_together(8,S,position, phi, r, copy_position, copy_phi, copy_r, data = data3, sign2 =sign3)
    data3[0][8]=data3[0][8]/10
    #define subset S for experiment 1 Far 12, heading
    S = np.concatenate((np.array([14,15,16,17,27,28,29]),np.array([0,1,2,12,13])))
    position, phi, r, copy_position, copy_phi, copy_r, data3 = all_together(9,S,position, phi, r, copy_position, copy_phi, copy_r, data = data3, sign2 =sign3)
    '''print(np.max(np.abs(r[:,30])))
    plt.plot(t_set,(180/pi)*copy_phi[:,30], color = 'b')
    for i in range(30):
        plt.plot(t_set,(180/pi)*copy_phi[:,i], color = 'r')
    plt.title("Heading for near 12 heading perturbed")
    plt.show()
    plt.plot(t_set,copy_r[:,30], color = 'b')
    for i in range(30):
        plt.plot(t_set,copy_r[:,i], color = 'r')
    plt.title("speed for near 12 heading perturbed")
    plt.show()'''

    #define subset S for experiment 1 Near 0,Far 0, speed
    S = np.array([])
    position, phi, r, copy_position, copy_phi, copy_r, data3 = all_together(0,S,position, phi, r, copy_position, copy_phi, copy_r, h=False, s =True, data = data3, sign2 =sign3)
    position, phi, r, copy_position, copy_phi, copy_r, data3 = all_together(1,S,position, phi, r, copy_position, copy_phi, copy_r, h=False, s =True, data = data3, sign2 =sign3)
    #define subset S for experiment 1 Near 3, speed
    #S = np.random.choice([0,1,2,12,13],3, replace=False)
    for S in combinations([0,1,2,12,13],3):
        position, phi, r, copy_position, copy_phi, copy_r, data3 = all_together(2,S,position, phi, r, copy_position, copy_phi, copy_r, h=False, s =True, data = data3, sign2 =sign3)
    data3[1][2]=data3[1][2]/10
    #define subset S for experiment 1 Near 6, speed
    #S = np.concatenate((np.array([0,1,2,12,13]),np.random.choice([14,15,16,17,27,28,29],1, replace=False)))
    for i in combinations([14,15,16,17,27,28,29],1):
        S = np.concatenate((np.array([0,1,2,12,13]), np.array(i)))
        position, phi, r, copy_position, copy_phi, copy_r, data3 = all_together(3,S,position, phi, r, copy_position, copy_phi, copy_r, h=False, s =True, data = data3, sign2 =sign3)
    data3[1][3]=data3[1][3]/7
    #define subset S for experiment 1 Near 9, speed
    #S = np.concatenate((np.array([0,1,2,12,13]),np.random.choice([14,15,16,17,27,28,29],4, replace=False)))
    for i in combinations([14,15,16,17,27,28,29],4):
        S = np.concatenate((np.array([0,1,2,12,13]), np.array(i)))
        position, phi, r, copy_position, copy_phi, copy_r, data3 = all_together(4,S,position, phi, r, copy_position, copy_phi, copy_r, h=False, s =True, data = data3, sign2 =sign3)
    data3[1][4]=data3[1][4]/35
    #define subset S for experiment 1 Near 12, speed
    S = np.concatenate((np.array([0,1,2,12,13]),np.array([14,15,16,17,27,28,29])))
    position, phi, r, copy_position, copy_phi, copy_r, data3 = all_together(5,S,position, phi, r, copy_position, copy_phi, copy_r, h=False, s =True, data = data3, sign2 =sign3)
    '''plt.plot(t_set,(180/pi)*copy_phi[:,30], color = 'b')
    for i in range(30):
        plt.plot(t_set,(180/pi)*copy_phi[:,i], color = 'r')
    plt.title("Heading for near 12 speed perturbed")
    plt.show()
    plt.plot(t_set,copy_r[:,30], color = 'b')
    for i in range(30):
        plt.plot(t_set,copy_r[:,i], color = 'r')
    plt.title("speed for near 12 speed perturbed")
    plt.show()
    '''
    #define subset S for experiment 1 Far 3, speed
    #S = np.random.choice([14,15,16,17,27,28,29],3, replace=False)
    for S in combinations([14,15,16,17,27,28,29],3):
        position, phi, r, copy_position, copy_phi, copy_r, data3 = all_together(6,S,position, phi, r, copy_position, copy_phi, copy_r, h=False, s =True, data = data3, sign2 =sign3)
    data3[1][6]=data3[1][6]/35
    #define subset S for experiment 1 Far 6, speed
    #S = np.random.choice([14,15,16,17,27,28,29],6, replace=False)
    for S in combinations([14,15,16,17,27,28,29],6):
        position, phi, r, copy_position, copy_phi, copy_r, data3 = all_together(7,S,position, phi, r, copy_position, copy_phi, copy_r, h=False, s =True, data = data3, sign2 =sign3)
    data3[1][7]=data3[1][7]/7
    #define subset S for experiment 1 Far 9, speed
    #S = np.concatenate((np.array([14,15,16,17,27,28,29]),np.random.choice([0,1,2,12,13],2, replace=False)))
    for i in combinations([0,1,2,12,13],2):
        S = np.concatenate((np.array([14,15,16,17,27,28,29]), np.array(i)))
        position, phi, r, copy_position, copy_phi, copy_r, data3 = all_together(8,S,position, phi, r, copy_position, copy_phi, copy_r, h=False, s =True, data = data3, sign2 =sign3)
    data3[1][8]=data3[1][8]/10
    #define subset S for experiment 1 Far 12, speed
    S = np.concatenate((np.array([14,15,16,17,27,28,29]),np.array([0,1,2,12,13])))
    position, phi, r, copy_position, copy_phi, copy_r, data3 = all_together(9,S,position, phi, r, copy_position, copy_phi, copy_r, h=False, s =True, data = data3, sign2 =sign3)
    
    global y_near_heading_my
    global y_near_speed_my
    global y_far_heading_my
    global y_far_speed_my
    y_near_heading_my = [data3[0][i] for i in [0,2,3,4,5]]
    y_far_heading_my = [data3[0][i] for i in [1,6,7,8,9]]
    y_near_speed_my = [data3[1][i]-0.3 for i in [0,2,3,4,5]]
    y_far_speed_my = [data3[1][i]-0.3 for i in [1,6,7,8,9]]
    return(y_near_heading_my,y_near_speed_my,y_far_heading_my,y_far_speed_my,data3)

def create_comparison_graph(y_near_heading_my, y_near_speed_my, y_far_heading_my, y_far_speed_my):
    plt.subplot(2, 1, 1)
    plt.plot(x,y_near_heading_paper,color = 'b',label = 'Paper near')
    plt.plot(x,y_far_heading_paper,color = 'r', label = 'paper far')
    plt.title(f'Final headings and speeds from paper, k={k}, a={a}')
    plt.ylabel('Final lateral deviation')

    plt.subplot(2, 1, 2)
    plt.plot(x,y_near_speed_paper,color = 'b', label = 'Paper near')
    plt.plot(x,y_far_speed_paper,color = 'r', label = 'Paper far')
    plt.ylabel('Final change in speed')
    plt.xlabel('Number of perturbed neighours')

    plt.subplot(2, 1, 1)
    plt.plot(x,y_near_heading_my,color = 'g',label = 'My near')
    plt.plot(x,y_far_heading_my,color = 'y',label = 'My far')
    plt.title(f'My Final headings and speeds, k={k}, a={a}, c={c}')
    plt.ylabel('Final lateral deviation')

    plt.subplot(2, 1, 2)
    plt.plot(x,y_near_speed_my,color = 'g', label = 'My near')
    plt.plot(x,y_far_speed_my,color = 'y', label = 'My far')
    plt.ylabel('Final change in speed')
    plt.xlabel('Number of perturbed neighours')

def similarity_func(y_near_speed_paper, y_near_speed_my, y_far_speed_paper, y_far_speed_my, y_near_heading_paper, 
                    y_near_heading_my, y_far_heading_paper,  y_far_heading_my, speed = True, heading = False):
    #This measures the similarity of the graphs and is big if they are very different
    sum = 0
    if speed == True:
        for i in range(5):
            sum += np.absolute(y_near_speed_paper[i] - y_near_speed_my[i])
            sum += np.absolute(y_far_speed_paper[i]  - y_far_speed_my[i])
    if heading == True:
        for i in range(5):
            sum += np.absolute(y_near_heading_paper[i] - y_near_heading_my[i])
            sum += np.absolute(y_far_heading_paper[i]  - y_far_heading_my[i])
    return sum

'''k,c = 150,200
y_near_heading_my,y_near_speed_my,y_far_heading_my,y_far_speed_my,my_exp1_final_heading_data = create_final_speed_heading(position, phi, r, copy_position, copy_phi, copy_r, data3=my_exp1_final_heading_data, sign3=1)
#y_near_heading_my,y_near_speed_my,y_far_heading_my,y_far_speed_my,my_exp1_final_heading_data_opposite = create_final_speed_heading(position, phi, r, copy_position, copy_phi, copy_r, data3=my_exp1_final_heading_data_opposite, sign3=-1)
create_comparison_graph(y_near_heading_my, y_near_speed_my, y_far_heading_my, y_far_speed_my)
plt.show()
#print(similarity_func(y_near_speed_paper, y_near_speed_my,y_far_speed_paper, y_far_speed_my, 
 #                                                   y_near_heading_paper, y_near_heading_my,y_far_heading_paper,  y_far_heading_my,
  #                                                  speed = False, heading = True))
k_range = np.array([100,200])
kc_results_speed = np.zeros(len(k_range))
kc_results_direction = np.zeros(len(k_range))

for i in range(len(k_range)):
    k,c = k_range[i], k_range[i]
    
    y_near_heading_my,y_near_speed_my,y_far_heading_my,y_far_speed_my,my_exp1_final_heading_data = create_final_speed_heading(position, phi, r, copy_position, copy_phi, copy_r, data3=my_exp1_final_heading_data, sign3=1)
    kc_results_speed[i] = similarity_func(y_near_speed_paper, y_near_speed_my, y_far_speed_paper, y_far_speed_my, 
                                            y_near_heading_paper, y_near_heading_my,y_far_heading_paper,  y_far_heading_my,
                                            speed = True, heading = False)
    kc_results_direction[i] = similarity_func(y_near_speed_paper, y_near_speed_my,y_far_speed_paper, y_far_speed_my, 
                                                y_near_heading_paper, y_near_heading_my,y_far_heading_paper,  y_far_heading_my,
                                                speed = False, heading = True)
    create_comparison_graph(y_near_heading_my, y_near_speed_my, y_far_heading_my, y_far_speed_my)
plt.show()
'''

number = 4
k_range = np.linspace(10,50,num=number)
c_range = np.linspace(10,50,num=number)
k_matrix = np.vstack((k_range,k_range,k_range,k_range,k_range,k_range))
c_matrix = np.transpose(k_matrix)
kc_results_speed = np.zeros((len(k_range),len(c_range)))
kc_results_direction = np.zeros((len(k_range),len(c_range)))
print(k_matrix)
print(c_matrix)
for i in range(number):
    for j in range(number):
        k = k_matrix[i,j]
        c = c_matrix[i,j]
        y_near_heading_my,y_near_speed_my,y_far_heading_my,y_far_speed_my,my_exp1_final_heading_data = create_final_speed_heading(position, phi, r, copy_position, copy_phi, copy_r, data3=my_exp1_final_heading_data, sign3=1)
        kc_results_speed[i,j] = similarity_func(y_near_speed_paper, y_near_speed_my, y_far_speed_paper, y_far_speed_my, 
                                                y_near_heading_paper, y_near_heading_my,y_far_heading_paper,  y_far_heading_my,
                                                speed = True, heading = False)
        kc_results_direction[i,j] = similarity_func(y_near_speed_paper, y_near_speed_my,y_far_speed_paper, y_far_speed_my, 
                                                    y_near_heading_paper, y_near_heading_my,y_far_heading_paper,  y_far_heading_my,
                                                    speed = False, heading = True)
        create_comparison_graph(y_near_heading_my, y_near_speed_my, y_far_heading_my, y_far_speed_my)
plt.show()

print("K, C matrix for speed")
print(kc_results_speed)
print("K, C Matrix for direction")
print(kc_results_direction)
#create_final called 36 times and within that all_toegther is called 216 times so count gets to 7776

end = time.time()
print("The time of execution of above program is :",
      (end-start), "s",(end-start)/60, "minutes" )