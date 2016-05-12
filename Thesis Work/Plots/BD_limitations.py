import numpy as np

MR = .92
mr = .15
u_0 = 4*np.pi*10**-7
R_shift = .018
dtheta = 11.25
mu0 = 4*np.pi*10**-7
points = np.append(np.arange(-9,9)+.5,np.arange(9,23)+.5)*11.25/180.*np.pi

sensor_R = np.append(MR+mr*np.cos((np.arange(-9,9)+0.5)*dtheta*np.pi/180.),
                     (MR-R_shift)+mr*np.cos((np.arange(9,23)+0.5)*dtheta*np.pi/180.))
sensor_Z = np.append(mr*np.sin((np.arange(-9,9)+0.5)*dtheta*np.pi/180.),
                     mr*np.sin((np.arange(9,23)+0.5)*dtheta*np.pi/180.))

def field_strength(plasma_R,I):
    r = np.sqrt((sensor_R - plasma_R)**2 + (sensor_Z)**2)
    B_p = u_0/2/np.pi*np.outer(1/r,I)
    return(B_p)

def eddy(t,sig_lev):
    eddy_strength = np.zeros(np.shape(sensor_R))

    eddy_strength[0] = .1*sig_lev
    eddy_strength[8] = .5*sig_lev
    eddy_strength[9] = .35*sig_lev
    eddy_strength[17] = .75*sig_lev

    eddy_decay = np.zeros(np.shape(sensor_R))

    eddy_decay[0] = 1/.002
    eddy_decay[8] = 1/.0005
    eddy_decay[9] = 1/.006
    eddy_decay[17] = 1/.015

    eddies = np.zeros((len(sensor_R),len(t)))

    pre_time = (len(t)//15)
    eddy_signal = np.exp(-np.outer(eddy_decay,(t[pre_time:]-t[pre_time])))
    for i in range(len(sensor_R)):
        eddy_signal[i] *= eddy_strength[i]
    eddies[:,pre_time:] += eddy_signal
    return(eddies)

def equilibrium(I,t):
    pre_time = len(t)//15
    post_time = -len(t)//25
    signal = np.zeros((len(sensor_R),len(t)))

    signal[:,pre_time:post_time] = field_strength(MR,I+t[pre_time:post_time])
    return(signal)

 
def mode(m,n):
    theta = np.arctan(sensor_Z, sensor_R)
    rad = np.sqrt(sensor_R**2 + sensor_Z**2)

def signal_eddies():
    I = 1500
    t = np.linspace(0,0.2,15000)
    signal = equilibrium(I,t)
    eddies = eddy(t,max(signal.reshape(-1)))

    tot = signal + eddies

    (topo,mag,chrono) = np.linalg.svd(tot,full_matrices = 0)

    return(tot,topo,mag,chrono)
