#!/usr/bin/python

<<<<<<< HEAD
#from Chamber_Shell_flux_and_PA_sig import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d
import tokamac

def pull_TokaMac():
    (X,Z,Psi) = tokamac.get_flux('./')
    [x,z,p] = [np.asarray(i) for i in [X,Z,Psi]]
    interp_func = interp2d(x,z,p)
    return(x,z,p,interp_func)

def axes():
    results = tokamac.get_run_results('./')
    mag_ax = np.array([results.mag_axis_r,results.mag_axis_z])
    rad = np.arange(-.17,.17,.00005)
    theta = np.arange(0,90)*np.pi/180.
    #coord1 is R for axis 1, Z for axis 2.
    #coord2 is Z for axis 1, -R for axis 2.
    coord1 = np.outer(np.cos(theta),rad)
    coord2 = np.outer(np.sin(theta),rad)

    return(theta,rad,coord1,coord2,mag_ax)

def consecutive(array):
    consec = [i == n for i,n in enumerate(array,array[0])]
    if all(consec):
        return(array[0],array[-1])
    if sum(consec) >= 0.5*len(array):
        return (array[0], array[consec.index(False)-1])
    if sum(consec) < 0.5*len(array):
        consec = [array[-1]-i == n for i,n in enumerate(array[::-1])]
        return(array[-consec.index(False)],array[-1])

def flux_boundary(f):
    (theta,rad,c1,c2,[mag_ax_r,mag_ax_z]) = axes()
    #(x,z,p,f) = pull_TokaMac()
    flux1 = np.zeros(np.shape(c1))
    flux2 = np.zeros(np.shape(c1))
    elong = np.zeros((90))

    edge_X_1 = []
    edge_X_2 = []
    edge_X_3 = []
    edge_X_4 = []
    edge_Z_1 = []
    edge_Z_2 = []
    edge_Z_3 = []
    edge_Z_4 = []
    (sep,lim)  = tokamac.get_separatrix('./')
    sep_flux = sep[0][2]
    for i in range(len(theta)):
        for j in range(len(c1[0])):
            flux1[i,j] = f(c1[i,j]+mag_ax_r,c2[i,j]+mag_ax_z)
            flux2[i,j] = f(-c2[i,j]+mag_ax_r,c1[i,j]+mag_ax_z)
    
    for i in range(len(theta)):
        good1 = np.where(flux1[i] < sep_flux)[0]
        good2 = np.where(flux2[i] < sep_flux)[0]
        
        (st1,end1) = consecutive(good1)
        (st2,end2) = consecutive(good2)

        dx1 = c1[i,st1] - c1[i,end1]
        dz1 = c2[i,st1] - c2[i,end1]
        edge_X_1.append(c1[i,end1]+mag_ax_r)
        edge_Z_1.append(c2[i,end1]+mag_ax_z)

        edge_X_2.append(-c2[i,end2]+mag_ax_r)
        edge_Z_2.append(c1[i,end2]+mag_ax_z)

        edge_X_3.append(c1[i,st1]+mag_ax_r)
        edge_Z_3.append(c2[i,st1]+mag_ax_z)

        edge_X_4.append(-c2[i,st2]+mag_ax_r)
        edge_Z_4.append(c1[i,st2]+mag_ax_z)

        dx2 = c2[i,st2] - c2[i,end2]
        dz2 = c1[i,st2] - c1[i,end2]

        d1 = np.sqrt(dx1**2 + dz1**2)
        d2 = np.sqrt(dx2**2 + dz2**2)
        elong[i] = d1/d2

    edge_X = edge_X_1+edge_X_2+edge_X_3+edge_X_4
    edge_Z = edge_Z_1+edge_Z_2+edge_Z_3+edge_Z_4

    return(edge_X,edge_Z,sep_flux)

def points(edge_X,edge_Z):

    MR = 0.5*(max(edge_X)+min(edge_X))
    a = 0.5*(max(edge_X)-min(edge_X))

    zoff = edge_Z[np.argmax(edge_X)]
    extrema = [(min(edge_X),edge_Z[np.argmin(edge_X)]),
               (max(edge_X),edge_Z[np.argmax(edge_X)]),
               (edge_X[np.argmin(edge_Z)],min(edge_Z)),
               (edge_X[np.argmax(edge_Z)],max(edge_Z))]

    return(MR,a,zoff,extrema)

def elongation(edge_Z,a,zoff):#mag_ax = None, edge_X=None,edge_Z=None):
    #if edge_X == None and edge_Z == None:
    #    (edge_X,edge_Z) = flux_boundary()

    #MR = 0.5*(max(edge_X)+min(edge_X))
    #a = 0.5*(max(edge_X)-min(edge_X))

    #zoff = edge_Z[np.argmax(edge_X)]

    elong_up = abs((max(edge_Z) - zoff)/a)
    elong_down = abs((min(edge_Z) - zoff)/a)

    return(elong_up,elong_down)

def triangularity(edge_X,edge_Z,MR,a):#edge_X = None, edge_Z = None):
    #if edge_X == None and edge_Z == None:
    #    (edge_X,edge_Z) = flux_boundary()
        
    #MR = 0.5*(max(edge_X)+min(edge_X))
    #a = 0.5*(max(edge_X)-min(edge_X))

    idx_zmax = np.argmax(edge_Z)
    idx_zmin = np.argmin(edge_Z)

    triang_up = (MR - edge_X[idx_zmax])/a
    triang_down = (MR - edge_X[idx_zmin])/a

    return(triang_up,triang_down)


def main():
    (x,z,p,f) = pull_TokaMac()
    (edge_X,edge_Z,sep_flux) = flux_boundary(f)
    (MR,a,zoff,extrema) = points(edge_X,edge_Z)
    (eup,edown) = elongation(edge_Z,a,zoff)
    (tup,tdown) = triangularity(edge_X,edge_Z,MR,a)
    return(eup,edown,tup,tdown,MR,a,zoff,edge_X,edge_Z,sep_flux,x,z,p)


def plot_flux_points():
    (eup,edown,tup,tdown,MR,a,zoff,edge_X,edge_Z,sep_flux,x,z,p) = main()
    fig = plt.figure()
    ax = fig.add_subplot(111,aspect = 'equal')
    flux_plt = plt.contour(x,z,p,200)
    sep_plt = plt.contour(x,z,p,[sep_flux])
    sep_plt_check = plt.plot(edge_X,edge_Z)
    plt.plot(max(edge_X),zoff,'o')
    plt.plot(edge_X[np.argmax(edge_Z)],max(edge_Z),'o')
    plt.plot(edge_X[np.argmin(edge_Z)],min(edge_Z),'o')
    plt.plot(min(edge_X),zoff,'o')

=======
from Chamber_Shell_flux_and_PA_sig import *
from scipy.interpolate import interp2d

def pull_TokaMac():
    (X,Z,Psi) = tokamac.get_flux('/home/byrne/Thesis_Work/mode_analysis/TokaMac/TokaMac_tutorial')
    (sep,lim)  = tokamac.get_separatrix('/home/byrne/Thesis_Work/mode_analysis/TokaMac/TokaMac_tutorial')
    results = tokamac.get_run_results('/home/byrne/Thesis_Work/mode_analysis/TokaMac/TokaMac_tutorial')
    [x,z,p] = [np.asarray(i) for i in [X,Z,Psi]]
    xx,zz = np.meshgrid(x,z)
    interp_func = interp2d(xx,zz,p)
    mag_ax = np.array([results.mag_axis_r,results_mag_axis_z])

    return(interp_func)

def 
>>>>>>> 04d5324d7e4cda4028c6deec3ddfd6dc74a4a004
