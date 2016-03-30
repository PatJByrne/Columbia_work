import numpy as np
import matplotlib.pyplot as plt

# set global settings
def init_plotting():
    plt.rcParams['font.size'] = 10
    #plt.rcParams['font.family'] = 'Times New Roman'
    plt.rcParams['axes.labelsize'] = 1.4*plt.rcParams['font.size']
    plt.rcParams['axes.titlesize'] = 2*plt.rcParams['font.size']
    plt.rcParams['legend.fontsize'] = 1.2*plt.rcParams['font.size']
    plt.rcParams['xtick.labelsize'] = 1.4*plt.rcParams['font.size']
    plt.rcParams['ytick.labelsize'] = 1.4*plt.rcParams['font.size']
    plt.rcParams['xtick.major.size'] = 3
    plt.rcParams['xtick.minor.size'] = 3
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['xtick.minor.width'] = 1
    plt.rcParams['ytick.major.size'] = 3
    plt.rcParams['ytick.minor.size'] = 3
    plt.rcParams['ytick.major.width'] = 1
    plt.rcParams['ytick.minor.width'] = 1
    plt.rcParams['legend.frameon'] = False
    plt.rcParams['legend.loc'] = 'upper right'
    plt.rcParams['axes.linewidth'] = 1


init_plotting()

MR = .92
mr = .15
sensor_rad = .155
u_0 = 4*np.pi*10**-7
R_shift = .018
dtheta = 11.25
mu0 = 4*np.pi*10**-7
t = np.linspace(0,0.02,15000)


def sensor_locations():
    sensor_nums = np.append(np.arange(7,7+18),
                            np.append(np.arange(7+18,32),np.arange(0,7)))
    segments = [(np.arange(-9,9)+.5)*dtheta/180.*np.pi,
                (np.arange(9,23)+.5)*dtheta/180.*np.pi]
    sensor_R = np.append(MR+sensor_rad*np.cos(segments[0]), 
                         (MR-R_shift)+sensor_rad*np.cos(segments[1]))
    sensor_Z = np.append(sensor_rad*np.sin(segments[0]), 
                         sensor_rad*np.sin(segments[1]))
    s_R = [None]*(len(sensor_R)+2)
    s_Z = [None]*(len(sensor_Z)+2)

    s_R[0] = sensor_R[sensor_nums.argsort()[-1]]
    s_Z[0] = sensor_Z[sensor_nums.argsort()[-1]]
    for num, i in enumerate(sensor_nums.argsort()):
        s_R[num+1] = sensor_R[i]
        s_Z[num+1] = sensor_Z[i]
    s_R[-1] = sensor_R[sensor_nums.argsort()[0]]
    s_Z[-1] = sensor_Z[sensor_nums.argsort()[0]]

    sensor_angles = np.linspace(-180,180-360/32.,32)+11.25*.5
    sensor_angles = np.append(sensor_angles[0]-11.25,sensor_angles)
    sensor_angles = np.append(sensor_angles, sensor_angles[-1]+11.25)
    return(sensor_nums,sensor_angles,np.asarray(s_R),np.asarray(s_Z))

def field_strength(plasma_R,I):
    r = np.sqrt((sensor_R - plasma_R)**2 + (sensor_Z)**2)
    B_p = u_0/2/np.pi*np.outer(1/r,I)
    return(B_p)

def eddy(sig_lev):
    eddy_strength = np.zeros(np.shape(sensor_R))

    eddy_strength[8] = .1*sig_lev
    eddy_strength[16] = .5*sig_lev
    eddy_strength[17] = .35*sig_lev
    eddy_strength[25] = .75*sig_lev

    eddy_decay = np.zeros(np.shape(sensor_R))

    eddy_decay[8] = 1/.002
    eddy_decay[16] = 1/.0005
    eddy_decay[17] = 1/.006
    eddy_decay[25] = 1/.009

    eddies = np.zeros((len(sensor_R),len(t)))

    pre_time = (len(t)//15)
    eddy_signal = np.exp(-np.outer(eddy_decay,(t[pre_time:]-t[pre_time])))
    for i in range(len(sensor_R)):
        eddy_signal[i] *= eddy_strength[i]
    eddies[:,pre_time:] += eddy_signal
    return(eddies)

def equilibrium(I,Maj_R = MR):
    pre_time = len(t)//15
    post_time = -len(t)//25
    signal = np.zeros((len(sensor_R),len(t)))

    signal[:,pre_time:post_time] = field_strength(Maj_R,I+t[pre_time:post_time])
    return(signal)

def plasma_edge(Maj_R = MR,min_R = mr,m=0,n=0,mag = 0):
    if Maj_R > .92:
        min_R = .15+(.92-Maj_R)
    elif Maj_R < .903:
        min_R = .15 + (Maj_R-.903)
    
    points = np.linspace(0,1,100)
    R = Maj_R + min_R*(np.cos(points*np.pi*2)) + mag*np.cos(m*points*np.pi*2)
    Z = min_R*(np.sin(points*np.pi*2)) + mag*np.sin(m*points*np.pi*2)
    return(Maj_R,min_R,R,Z)
 
def mode(m,n,Maj_R = MR,min_R = mr):
    phi = 0
    freq = 2*np.pi*1000
    theta = np.zeros((len(sensor_R),len(t)))
    rad = np.zeros((len(sensor_R),len(t)))
    signal = np.zeros(np.shape(theta))
    
    pre_time = len(t)//15
    post_time = -len(t)//25

    #ampl = .00025*(np.exp((t[pre_time:post_time]-t[pre_time])/.01)-1)
    ampl = .005*np.sqrt(t[pre_time:post_time]-t[pre_time])
    

    for i in range(len(sensor_R)):
        theta[i,:] += np.arctan2(sensor_Z[i], sensor_R[i]-Maj_R)
        rad[i,:] += np.sqrt((sensor_R[i]-Maj_R)**2 + sensor_Z[i]**2)
    
        signal[i,pre_time:post_time] = ampl*(min_R/rad[i,pre_time:post_time])**(2*m)*np.sin(n*phi+m*theta[i,pre_time:post_time] - freq*t[pre_time:post_time])
    return(theta,signal)


def signals(Maj_R = MR,shift = 0,pre_time=len(t)/15,post_time = -len(t)//25):
    I = 1500

    (sensor_num,sensor_angles,sensor_R,sensor_Z) = sensor_locations()
    signal = equilibrium(I,Maj_R)
    eddies = eddy(max(signal.reshape(-1)))
    (Maj_R,min_R,R,Z) = plasma_edge(Maj_R = Maj_R)
 
    (theta,fluct) = mode(3,1,Maj_R)
    tot = signal+eddies + fluct

    (topo,mag,chrono) = np.linalg.svd(tot[:,(pre_time+shift):post_time],full_matrices = 0)

    return(tot,topo,mag,chrono,signal,eddies,fluct,t[pre_time+shift:post_time])

(sensor_nums,sensor_angles,sensor_R,sensor_Z) = sensor_locations()
(tot,topo,mag,chrono,signal,eddies,fluct,bdt) = signals(Maj_R = .95, shift = 10000,pre_time = 0,post_time = -1000)
fig = plt.figure()
plt.plot(mag,'o')
fig = plt.figure()
plt.plot(topo[:,0],'o')
plt.plot(topo[:,1],'o')
plt.plot(topo[:,2],'o')
plt.plot(topo[:,3],'o')
fig = plt.figure()
plt.plot(sensor_R,sensor_Z,'o')
plt.plot(sensor_R[1],sensor_Z[1],'o')
plt.plot(.95,0,'o')

fig = plt.figure(figsize = (30,15))
ax = fig.add_subplot(4,2,1)
plt.plot(t,signal[17])
plt.yticks([])
plt.ylim(-.0005,.0045)
plt.ylabel('arb. units')
plt.xticks([])
plt.xlim(0,.02)
ax = fig.add_subplot(4,2,2)
plt.plot(t,eddies[25]+.0025,'b')
plt.plot(t,eddies[17]+.00175,'g--')
plt.plot(t,eddies[16]+.00075,'r-.')
plt.plot(t,eddies[8] +.0005,'k:')
plt.xlim(0,.02)
plt.ylim(0,.005)
plt.yticks([])
plt.xticks([])
ax = fig.add_subplot(4,2,3)
plt.contourf(t,sensor_angles,signal,50,cmap = plt.cm.gray)
#plt.yticks([])
plt.ylabel('Poloidal Angle')
plt.xticks([])
plt.ylim(-180,180)
#plt.colorbar()
ax = fig.add_subplot(4,2,4)
plt.contourf(t,sensor_angles,eddies,50,cmap = plt.cm.gray)
plt.yticks([])
plt.ylim(-180,180)
plt.xticks([])
#plt.colorbar()
ax = fig.add_subplot(4,2,5)
plt.plot(t,fluct[17])
plt.yticks([])
plt.ylim(-.0025,.0025)
plt.ylabel('arb. units')
plt.xticks([])
plt.xlim(0,.02)
ax = fig.add_subplot(4,2,6)
plt.plot(t,tot[17])
plt.yticks([])
plt.xticks([])
plt.ylim(-.0005,.0045)
plt.xlim(0,.02)
ax = fig.add_subplot(4,2,7)
plt.contourf(t*1000,sensor_angles,fluct,50,cmap = plt.cm.gray)
#plt.yticks([])
plt.ylim(-180,180)
plt.ylabel('Poloidal Angle')
plt.xlabel('time (ms)')
#plt.colorbar()
ax = fig.add_subplot(4,2,8)
plt.contourf(t*1000,sensor_angles,tot,50,cmap = plt.cm.gray)
plt.yticks([])
plt.ylim(-180,180)
plt.xlabel('time (ms)')
#plt.colorbar()
fig.subplots_adjust(left = .04,right = .99,top = .99, bottom = .05,wspace = .01,hspace = .02)

fig = plt.figure(figsize = (30,15))
ax = fig.add_subplot(4,2,1)
plt.plot(t,signal[17])
plt.plot(bdt[::5],(topo[17,0]*mag[0]*chrono[0])[::5],linewidth = 3)
plt.yticks([])
plt.ylim(-.0005,.004)
plt.ylabel('arb. units')
plt.xticks([])
plt.xlim(0,.02)
ax = fig.add_subplot(4,2,2)
plt.plot(t,eddies[25]+.0025,'b')
plt.plot(t,eddies[17]+.00175,'g--')
plt.plot(t,eddies[16]+.00075,'r-.')
plt.plot(t,eddies[8] +.0005,'k:')
plt.plot(bdt[::5],(topo[25,3]*mag[3]*chrono[3]+.0025)[::5],'b',linewidth = 3)
plt.plot(bdt[::5],(topo[17,3]*mag[3]*chrono[3]+.00175)[::5],'g--',linewidth = 3)
plt.plot(bdt[::5],(topo[16,3]*mag[3]*chrono[3]+.00075)[::5],'r-.',linewidth = 3)
plt.plot(bdt[::5],(topo[8,3]*mag[3]*chrono[3]+.0005)[::5],'c:',linewidth = 3)

plt.xlim(0,.02)
plt.ylim(0,.005)
plt.yticks([])
plt.xticks([])
ax = fig.add_subplot(4,2,3)
plt.contourf(t,sensor_angles,signal,50,cmap = plt.cm.gray)
#plt.yticks([])
plt.ylabel('Poloidal Angle')
plt.xticks([])
plt.ylim(-180,180)
#plt.colorbar()
ax = fig.add_subplot(4,2,4)
plt.contourf(t,sensor_angles,eddies,50,cmap = plt.cm.gray)
plt.yticks([])
plt.ylim(-180,180)
plt.xticks([])
#plt.colorbar()
ax = fig.add_subplot(4,2,5)
#plt.plot(bdt,mag[1]*chrono[1])
#plt.plot(bdt,mag[2]*chrono[2])
plt.plot(t,fluct[17])
plt.plot(bdt[::5],(topo[17,1]*mag[1]*chrono[1]+topo[17,2]*mag[2]*chrono[2])[::5],linewidth = 3)
plt.yticks([])
plt.ylim(-.0025,.0025)
plt.ylabel('arb. units')
plt.xticks([])
plt.xlim(0,.02)
ax = fig.add_subplot(4,2,6)
plt.plot(t,tot[17])
plt.plot(bdt[::5],(topo[17,0]*mag[0]*chrono[0]+topo[17,3]*mag[3]*chrono[3]+topo[17,1]*mag[1]*chrono[1]+topo[17,2]*mag[2]*chrono[2])[::5],linewidth = 3)
plt.yticks([])
plt.xticks([])
plt.ylim(-.0005,.0045)
plt.xlim(0,.02)
ax = fig.add_subplot(4,2,7)
plt.contourf(t*1000,sensor_angles,fluct,50,cmap = plt.cm.gray)
#plt.yticks([])
plt.ylim(-180,180)
plt.ylabel('Poloidal Angle')
plt.xlabel('time (ms)')
#plt.colorbar()
ax = fig.add_subplot(4,2,8)
plt.contourf(t*1000,sensor_angles,tot,50,cmap = plt.cm.gray)
plt.yticks([])
plt.ylim(-180,180)
plt.xlabel('time (ms)')
#plt.colorbar()
fig.subplots_adjust(left = .04,right = .99,top = .99, bottom = .05,wspace = .02,hspace = .02)

fig = plt.figure()

plt.show()

