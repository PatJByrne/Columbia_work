#!/usr/bin/python

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
