#!/usr/bin/python

import numpy as np
from MDSplus import *
import cPickle as pickle
import sys
import pylab
import os
import matplotlib
from matplotlib import pyplot as plt
sys.path.append('/home/byrne/Thesis_Work/mode_analysis')
from subtraction_mode_analysis import crash_finder
sys.path.append('/home/byrne/Thesis_Work/Simulations/Bfields')
import Synthetic_sensors
sys.path.append('/home/byrne/APS_Work/2011_Shaping_Coil/Circuit_simulations/Open_Positions')
sys.path.append('/home/byrne/APS_Work/2011_Shaping_Coil/flux_surface_plots')
sys.path.append('/home/byrne/APS_Work/2011_Shaping_Coil/Circuit_simulations')
sys.path.append('/home/byrne/repos/hbt/python/hbtep')
import major_radius
import Ip
import Vf
import Oh
import SH
import saddle_finder
import miscmod
import tokamac
#import Decay_index
#import magnetic_surfaces

numpoints = 10000
Radius = np.linspace(60,160,numpoints)
Zed    = np.linspace(-50,50,numpoints)

def plot2(IS,IV,IO,IP,shotnum,time):
    fig = pylab.figure()
    ax = fig.add_subplot(111,aspect = 'equal')

    ((R,Z),((chamberr,chamberz),(shellr,shellz),(flanger,flangez),(limiterr,limiterz),(PA_in_R,PA_in_Z),(SHR,SHZ)),(s,p,o,v),(sc,pc,oc,vc)) = coil_config(shotnum,time)

    pylab.plot(SHR,SHZ,'og')
    pylab.plot(flanger,flangez,'k')
    pylab.plot(shellr,shellz,'k')
    pylab.plot(PA_in_R,PA_in_Z,'k')
    pylab.plot(limiterr,limiterz,'k')
    flux = (v*IV+s*IS+o*IO+IP*p)
    (r_point,z_point,saddle_flux) = saddle_finder.saddle_finder(flux[500:600,0:160],(R[0:160],Z[500:600]))
    pylab.contour(R,Z,flux,np.append((vc*IV+sc*IS+IO*oc+IP*pc)[3:10],saddle_flux))

def Currents_fluxes(shotnum,time,real_sig = True):
    tree = Tree('hbtep2',shotnum)
    t1 = tree.getNode('timing.banks.OH_Bias_st').data()*1e-6
    #t2 = 6

    #(t_MR,mr,MR,q) = major_radius.getMajorRadius_Q(shotnum)
    #MR_smooth = major_radius.boxcar(MR,5)
    #mr = miscmod.calc_r_minor(MR_smooth*.01)*100
    #q = miscmod.calc_q(tree,t_MR,r_major = MR_smooth)

    (t_o,io) = Oh.Oh(shotnum)
    (t_v,iv) = Vf.Vf(shotnum)
    (t_i,ip) = Ip.Ip(shotnum)
    (t_s,ish) = SH.Sh(shotnum)
 

    tquench = crash_finder(shotnum,t_i)[0]
    t2 = tquench*1000+.2

    MR_t = t_i[np.where(t_i<tquench+.001)]
    #print np.shape(MR_t)
    (q,MR) = miscmod.calc_q(tree, times = t_i[np.where(t_i<tquench+.001)], 
                            byrne = True)
    #print(len(MR))
    MR_ma = np.ma.masked_array(MR,MR_t<.0009)
    IO = io[np.argmin(abs(t_o-time))]
    IV = iv[np.argmin(abs(t_v-time))]
    IS = .01*ish[np.argmin(abs(t_s-time))]
    IP =  ip[np.argmin(abs(t_i-time))]
    MR_point = MR[np.argmin(abs(t_i-time))]# MR_smooth[np.argmin(abs(t_MR-time))]
    currents  = (IV,IS,IO,IP)

    (X,Z,Psi) = tokamac.get_flux('/home/byrne/Thesis_Work/mode_analysis/TokaMac/TokaMac_tutorial')
    (sep,lim)  = tokamac.get_separatrix('/home/byrne/Thesis_Work/mode_analysis/TokaMac/TokaMac_tutorial')
    Psi = np.asarray(Psi)
    center = np.round(X[np.argmin(Psi[128].reshape(-1))],2)
    contour_locations = np.arange(center,1.07,.01)
    contour_index = []
    for i in contour_locations:
        contour_index.append(np.argmin(np.abs(X-i)))
    edge = lim
    for i in sep:
        if i[2] < edge[0]:
            edge = [i[2]]

    (flux_nom,saddle_flux_nom,flux_nomc,((R,z),((chamberr,chamberz),(shellr,shellz),(flanger,flangez),((limiterr,limiterz),(top_lim_r,top_lim_z),(bot_lim_r,bot_lim_z),(out_lim_r,out_lim_z)),(PA_in_R,PA_in_Z),(SHR,SHZ)))) = APS_plot2(shotnum,time,currents,X,Z,Psi)

    for i in range(len(flanger)):
       d_x = X-flanger[i]
       if len(d_x[d_x>0]) == 0:
           continue
       x_i = list(d_x).index(min(d_x[d_x>0]))
       d_z = (Z-flangez[i])*np.sign(flangez[i])
       if len(d_z[d_z>0]) == 0:
           continue
       z_i = list(d_z).index(min(d_z[d_z>0]))
       
       if Psi[z_i,x_i]<edge[0]:
           edge = [Psi[z_i,x_i]]
       
       
    #print contour_index
    #print('currents '), currents

    #((R,Z),((chamberr,chamberz),(shellr,shellz),(flanger,flangez),((limiterr,limiterz),(top_lim_r,top_lim_z),(bot_lim_r,bot_lim_z),(out_lim_r,out_lim_z)),(PA_in_R,PA_in_Z),(SHR,SHZ)),(s,p,o,v),(sc,pc,oc,vc))
    #(flux_sxr,saddle_flux_sxr,flux_sxrc,((R,Z),((chamberr,chamberz),(shellr,shellz),(flanger,flangez),(PA_in_R,PA_in_Z),(SHR_sxr,SHZ_sxr)))) = APS_plot2_sxr(shotnum,time,currents)

    fig7 = pylab.figure(7,figsize = (22,18))
    pylab.figtext(.035,.925,'a)',fontsize = 'large')
    #ax1 = fig7.add_subplot(331)
    pylab.figtext(.03,.975,'Shot number {}, time {}ms  PA1  Bfield (G)  PA2  Bfield (G)  '.format(shotnum,time*1e3))
    ax1 = fig7.add_subplot(421)
    pylab.xlim([t1,t2])
    pylab.ylim([0,15])
    pylab.plot(t_i*1000,ip*.001,label = 'Ip')
    pylab.plot(time*1000,IP*.001,'oy')
    pylab.plot(t_s*1000,ish*.001,'--',label = 'Ish')
    pylab.plot(time*1000,IS*.1,'oy')
    #pylab.title('Plasma Current')    
    pylab.ylabel('Current (kA)')
    pylab.legend(loc =3)
    ax1.set_xticklabels([])

    #ax3 = fig7.add_subplot(334)
    ax3 = fig7.add_subplot(423)
    pylab.xlim([t1,t2])
    pylab.ylim([88,98])
    #t_MR = t_i[np.where(t_i>.001)]
    #pylab.plot(t_MR*1000,MR[np.where(t_i>.001)]*100)
    pylab.plot(t_i[np.where(t_i<tquench+.001)]*1000,MR_ma*100)#_smooth)
    #pylab.plot(t_MR*1000,MR_smooth*.01)   
    #pylab.plot(t_MR*1000,MR_smooth*100)
    pylab.plot(time*1000,MR_point*100,'oy')
    pylab.plot([t1,t2],[92,92],'k:')
    pylab.ylabel('MR (cm)')
    pylab.xlabel('Time (ms)')
    
    fig7.subplots_adjust(left = .032, right = 1, top = .95, bottom = .075,
                         hspace = .1,wspace = .025)
    
    #ax2 = fig7.add_subplot(337)
    #ax2 = fig7.add_subplot(122)
    #pylab.xlim([t1,t2])
    #pylab.ylim([2,4])
    #pylab.plot(t_i[np.where(t_i<tquench+.001)]*1000,q)
    #pylab.plot(time*1000,q[np.argmin(abs(t_i-time))],'oy')
    #pylab.plot([t1,t2],[3,3],'k:')
    #pylab.ylabel('Safety Factor (q*)')
    #pylab.xlabel('Time (ms)')
    
    #fig6 = pylab.figure(6)
    pylab.figtext(.535,.925,'b)',fontsize = 'large')
    #ax4 = fig7.add_subplot(132,aspect = 'equal')
    ax2 = fig7.add_subplot(122,aspect = 'equal')
    
    #pylab.title('Shot number {}, time {}ms'.format(shotnum,time*1e3))
    ((PA_r,PA_z,PA_phi),(TA_r,TA_z,TA_phi)) = Synthetic_sensors.positions()
    for i in range(6):
        pylab.plot(PA_r[-1-i],PA_z[-1-i],'s',markersize = 7)
    pylab.plot(SHR,SHZ,'og')
    pylab.plot(flanger,flangez,'k')
    pylab.plot(limiterr,limiterz,'k')
    pylab.plot(top_lim_r,top_lim_z,'k')
    pylab.plot(bot_lim_r,bot_lim_z,'k')
    pylab.plot(out_lim_r,out_lim_z,'k')
    pylab.plot(shellr,shellz,'k')
    pylab.plot(PA_in_R,PA_in_Z,'k')
    #pylab.contour(R,Z,flux_nom, np.append(flux_nomc[3:10],saddle_flux_nom) )
    #pylab.contour(R,Z,flux_nom, [saddle_flux_nom],colors = 'r',linestyles = 'dashed',linethickness = 5)
    #print Psi[128,contour_index]
    #print len(contour_index),len(Psi[128,contour_index])

    pylab.contour(X,Z,Psi,Psi[128,contour_index])# np.append(Psi[128,contour_index],sep[0][2]) )
    for i in sep:
        if edge[0] >= i[2]:
            pylab.contour(X,Z,Psi,[i[2]],colors = 'r',linestyles = 'dashed',linethickness = 5)
    #pylab.contour(X,Z,Psi,[saddle_flux_nom],colors = 'r',linestyles = 'dashed',linethickness = 5)
    pylab.contour(X,Z,Psi,edge,colors = 'b',linestyles = 'dashed',linethickness = 5)
    ax2.set_xticklabels([])
    ax2.set_yticklabels([])
    #pylab.contour(R,Z,flux_nom, [saddle_flux_nom,saddle_flux_sxr],colors = 'g',linestyles = 'dashed')
    #pylab.xlabel('R (m)')
    pylab.xlim([.675,1.08])
    #pylab.ylabel('Z (m)')
    pylab.ylim([-.22,.22])
    fig7.subplots_adjust(left = .032, right = 1, top = .95, bottom = .075,hspace = .1,wspace = .025)

    ylimit = [-250,400]
    ((t_o,Br_PA,Bth_PA,Br_TA,Bth_TA),
     ((PA1_sensor_p, PA1_sensor_r, PA2_sensor_p, PA2_sensor_r),
      (TA_sensor_p, TA_sensor_r))) = Synthetic_sensors.main(shotnum)

    #fig5 = pylab.figure(5)
    pylab.figtext(.035,.475,'c)',fontsize = 'large')
    ax1 = fig7.add_subplot(425)
    #ax1 = fig5.add_subplot(211)
    #fig5.subplots_adjust(left = .032, right = 1, top = .95, bottom = .075,hspace = .1,wspace = .025)
    colors = ['b','g','r','c','m','y']
    for i in range(6):
        pylab.plot(t_o*1e3,Bth_PA[31-i],colors[i]+':')
        #pylab.plot(t_o,Br_PA[-2-i*2],colors[i])
        if real_sig:
            pylab.plot(t_o*1e3,PA1_sensor_p[-1-i]*1e4,colors[i])
    pylab.plot([t1,t2],[0,0],'k:')
    pylab.fill_between([time*1e3-.05,time*1e3+.05],-200,400,facecolor = 'k',alpha=.2)

    pylab.xlim([t1,t2])
    pylab.xticks([])
    #pylab.title('PA1')
    pylab.ylabel('PA1  Bfield (G)')
    pylab.ylim([-200,400])
    ax2 = fig7.add_subplot(427)
    #ax2 = fig5.add_subplot(212)
    pylab.plot([time,time],[-2,4],'k:')
    for i in range(5):
        pylab.plot(t_o*1e3,Bth_PA[31-i],colors[i]+':')
        #pylab.plot(t_o,Br_PA[-2-i*2],colors[i])
        if real_sig:
            pylab.plot(t_o*1e3,PA2_sensor_p[-1-i]*1e4,colors[i])
    pylab.plot([t1,t2],[0,0],'k:')
    pylab.fill_between([time*1e3-.05,time*1e3+.05],-200,400,facecolor = 'k',alpha=.2)
    pylab.xlim([t1,t2])
    #pylab.title('PA2')
    pylab.ylabel('PA2   Bfield (G)')
    pylab.ylim([-200,400])
    pylab.xlabel('Time (ms)')
    #for i in range(4):
    #    pylab.plot(t_o,PA2_sensor_p[31-i],':')
      #pylab.plot(t_o,PA2_sensor_r[15-i]*1e4,':')
        
    #pylab.ylim(ylimit)
    #pylab.xlim(-.001,.025)
    #pylab.figtext(.25,.922,'solid simulation, dots measured')

   #fig2 = pylab.figure(figsize = (11,8))
   #for i in range(6):
   #   pylab.plot(t_o,PA1_sensor_p[31-i]*1e4)
   #pylab.ylim(ylimit)
   #pylab.twinx()
   #for i in range(6):
   #   pylab.plot(t_o,PA2_sensor_p[31-i]*1e4,':')
   #pylab.ylim(ylimit)
   #pylab.xlim(-.001,.025)   
   #pylab.figtext(.25,.922,'solid PA1, dots PA2')


    pylab.savefig('flux_surfaces_during_shaping'+str(shotnum))

def APS_plot2(shotnum,time,(IV,IS,IO,IP),X,Z,Flux):
    

    ((R,z),((chamberr,chamberz),(shellr,shellz),(flanger,flangez),((limiterr,limiterz),(top_lim_r,top_lim_z),(bot_lim_r,bot_lim_z),(out_lim_r,out_lim_z)),(PA_in_R,PA_in_Z),(SHR,SHZ)),(s,p,o,v),(sc,pc,oc,vc)) = coil_config(shotnum,time)
    #flux = (v*IV+s*IS+o*IO+IP*p)
    contours = (vc*IV+sc*IS+oc*IO+IP*pc)
    center =  int(np.round(len(X)*.5))
    print center
    #(r_point,z_point,saddle_flux) = saddle_finder.saddle_finder(flux[500:600,0:160],(R[0:160],Z[500:600]))
    (r_point,z_point,saddle_flux) = saddle_finder.saddle_finder(Flux[center:,0:center],(X[0:center],Z[center:]))
    print 'Saddle point @: ',r_point,z_point
    coords = ((R,Z),((chamberr,chamberz),(shellr,shellz),(flanger,flangez),((limiterr,limiterz),(top_lim_r,top_lim_z),(bot_lim_r,bot_lim_z),(out_lim_r,out_lim_z)),(PA_in_R,PA_in_Z),(SHR,SHZ)))
    #pylab.contour(R,Z,flux,saddle_flux,'g:')
    return (Flux,saddle_flux,contours,coords)
       #if not i%4:
       #    pylab.yticks(fontsize = 6)
       #else:
       #    ax.set_yticklabels([])
       #if i/4 == 4:
       #    pylab.xticks(fontsize = 6)

def APS_plot2_sxr(shotnum,time,(IV,IS,IO,IP)):

    ((R,Z),((chamberr,chamberz),(shellr,shellz),(flanger,flangez),(limiterr,limiterz),(PA_in_R,PA_in_Z),(SHR,SHZ)),(s,p,o,v),(sc,pc,oc,vc)) = coil_config_sxr(shotnum,time)
    flux = (v*IV+s*IS+o*IO+IP*p)
    contours = (vc*IV+sc*IS+oc*IO+IP*pc)
    (r_point,z_point,saddle_flux) = saddle_finder.saddle_finder(flux[500:600,0:160],(R[0:160],Z[500:600]))
    coords = ((R,Z),((chamberr,chamberz),(shellr,shellz),(flanger,flangez),(PA_in_R,PA_in_Z),(SHR,SHZ)))
    return (flux,saddle_flux,contours,coords)

def MR(shotnum,t):
    tree = Tree('hbtep2',shotnum)
    (time,r_major) = miscmod.calc_r_major(tree,t_start = 0,t_end = .01)
    index = np.argmin(abs(time-t))
    return r_major[index]*1000

def coil_config(shotnum,time):
    r_major = MR(shotnum,time)
    mr = '{0:03d}'.format(int(r_major))

    f = open('/home/byrne/APS_Work/2011_Shaping_Coil/flux_surface_plots/plots/Final_Design/shaping_0100A_ellip_flux','r')

    (sa,s,R,Z,sc) = pickle.load(f)

    f.close()

    f = open('/home/byrne/APS_Work/2011_Shaping_Coil/flux_surface_plots/plots/plasma_ellip_flux_1A_'+mr+'mm','r')

    (pa,p,R,Z,pc) = pickle.load(f)

    f.close()

    f = open('/home/byrne/APS_Work/2011_Shaping_Coil/flux_surface_plots/plots/oh_ellip_flux_1A','r')   

    (oa,o,R,Z,oc) = pickle.load(f)

    f.close()

    f = open('/home/byrne/APS_Work/2011_Shaping_Coil/flux_surface_plots/plots/vf_ellip_flux_1A','r')   
    
    (va,v,R,Z,vc) = pickle.load(f) 

    f.close()

    f = open('/home/byrne/APS_Work/2011_Shaping_Coil/flux_surface_plots/plots/Final_Design/limiter_R_Z','r')
    
    ((limiterr,limiterz),(top_lim_r,top_lim_z),(bot_lim_r,bot_lim_z),(out_lim_r,out_lim_z)) = pickle.load(f)

    f.close()

    f = open('/home/byrne/APS_Work/2011_Shaping_Coil/flux_surface_plots/plots/Final_Design/chamber_R_Z','r')
    
    (chamberr,chamberz) = pickle.load(f)

    f.close()

    f = open('/home/byrne/APS_Work/2011_Shaping_Coil/flux_surface_plots/plots/Final_Design/shell_R_Z','r')
    
    (shellr,shellz) = pickle.load(f)

    f.close()

    f = open('/home/byrne/APS_Work/2011_Shaping_Coil/flux_surface_plots/plots/Final_Design/flange_R_Z','r')

    (flanger,flangez) = pickle.load(f)

    f.close()

    f = open('/home/byrne/APS_Work/2011_Shaping_Coil/flux_surface_plots/plots/Final_Design/PA_in_R_Z','r')
    
    (PA_in_R,PA_in_Z) = pickle.load(f)

    f.close()

    f = open('/home/byrne/APS_Work/2011_Shaping_Coil/flux_surface_plots/plots/Final_Design/coil_R_Z','r') 
    
    ((OHR,OHZ),(VFR,VFZ),(SHR,SHZ)) = pickle.load(f)

    f.close()

    return((R,Z),((chamberr,chamberz),(shellr,shellz),(flanger,flangez),((limiterr,limiterz),(top_lim_r,top_lim_z),(bot_lim_r,bot_lim_z),(out_lim_r,out_lim_z)),(PA_in_R,PA_in_Z),(SHR,SHZ)),(s,p,o,v),(sc,pc,oc,vc))

def coil_config_sxr(shotnum,time):
    r_major = MR(shotnum,time)
    mr = '{0:03d}'.format(int(r_major))

    f = open('/home/byrne/APS_Work/2011_Shaping_Coil/flux_surface_plots/plots/Final_Design/shaping_SXR_displaced_0100A_ellip_flux','r')

    (sa,s,R,Z,sc) = pickle.load(f)

    f.close()

    f = open('/home/byrne/APS_Work/2011_Shaping_Coil/flux_surface_plots/plots/plasma_ellip_flux_1A_'+mr+'mm','r')
    
    (pa,p,R,Z,pc) = pickle.load(f)

    f.close()

    f = open('/home/byrne/APS_Work/2011_Shaping_Coil/flux_surface_plots/plots/oh_ellip_flux_1A','r')   

    (oa,o,R,Z,oc) = pickle.load(f)

    f.close()

    f = open('/home/byrne/APS_Work/2011_Shaping_Coil/flux_surface_plots/plots/vf_ellip_flux_1A','r')   
    
    (va,v,R,Z,vc) = pickle.load(f) 

    f.close()

    f = open('/home/byrne/APS_Work/2011_Shaping_Coil/flux_surface_plots/plots/Final_Design/chamber_R_Z','r')
    
    (chamberr,chamberz) = pickle.load(f)

    f.close()

    f = open('/home/byrne/APS_Work/2011_Shaping_Coil/flux_surface_plots/plots/Final_Design/limiter_R_Z','r')
    
    (limiterr,limiterz) = pickle.load(f)

    f.close()
    f = open('/home/byrne/APS_Work/2011_Shaping_Coil/flux_surface_plots/plots/Final_Design/shell_R_Z','r')
    
    (shellr,shellz) = pickle.load(f)

    f.close()

    f = open('/home/byrne/APS_Work/2011_Shaping_Coil/flux_surface_plots/plots/Final_Design/flange_R_Z','r')

    (flanger,flangez) = pickle.load(f)

    f.close()

    f = open('/home/byrne/APS_Work/2011_Shaping_Coil/flux_surface_plots/plots/Final_Design/PA_in_R_Z','r')
    
    (PA_in_R,PA_in_Z) = pickle.load(f)

    f.close()

    f = open('/home/byrne/APS_Work/2011_Shaping_Coil/flux_surface_plots/plots/Final_Design/coil_R_Z_sxrmod','r') 
    
    ((OHR,OHZ),(VFR,VFZ),(SHR,SHZ)) = pickle.load(f)

    f.close()

    return((R,Z),((chamberr,chamberz),(shellr,shellz),(flanger,flangez),(limiterr,limiterz),(PA_in_R,PA_in_Z),(SHR,SHZ)),(s,p,o,v),(sc,pc,oc,vc))
