#!/usr/bin/python

#Stepping through the field
############################################################################
import numpy as np
import pylab
import cPickle as pickle
from MDSplus import *
import os
import sys
from scipy.interpolate import interp2d
sys.path.append('/home/byrne/APS_Work/2011_Shaping_Coil/flux_surface_plots/plots/Final_Design')
import plotter

sys.path.append('/home/byrne/APS_Work/2011_Shaping_Coil/flux_surface_plots')
import magnetic_surfaces

sys.path.append('/home/byrne/repos/hbt/python/hbtep')
import major_radius
sys.path.append('/home/byrne/APS_Work/2011_Shaping_Coil/Circuit_simulations')
import Ip
import Oh
import Vf
import SH
import Sensors
sys.path.append('/home/byrne/APS_Work/2011_Shaping_Coil/decay_index')
import Bfield
sys.path.append('/home/byrne/repos/hbt/python/hbtep')
import major_radius
import miscmod

#Eqilibium Fields
#####################################################################

numpoints = 2000
#Radius = np.linspace(.60,1.60,numpoints)
#Zed    = np.linspace(-.50,.50,numpoints)
Radius = np.linspace(.72,1.12,numpoints)
#Radius = np.linspace(-0.000001,1.12,numpoints)
Zed    = np.linspace(-.20,.20,numpoints)
Gauss  = 10000
degrees = np.pi/180.
inches = .0254
(R,Z) = np.meshgrid(Radius,Zed)

def main(shotno,time):

   ((PA_R,PA_Z,PA_phi),(TA_R,TA_Z,TA_phi)) = positions()

   ((brpa,bzpa),(brta,bzta)) = fields(shotno,(PA_R,PA_Z),(TA_R,TA_Z),time)

   (timepa, PA1_sensor_p, PA1_sensor_r,
    PA2_sensor_p, PA2_sensor_r) = Sensors.PA(shotno)

   (timeta, TA_sensor_p, TA_sensor_r) = Sensors.TA(shotno)

   #The convention for HBT-EP sensors is r-hat is into the plasma, and 
   #theta-hat is clockwise (-r and -theta wrt plasma)
   #TA - positive r is in the R direction in RZphi space
   #TA - positive theta is in the Z direction in RZphi space

   fig1 = pylab.figure(1)
   pylab.title('Ideal Vacuum Fields (R) vs. TA measurements @ {0:2}ms, shot# {1:5}'.format(time*1000,shotno))
   pylab.xlabel('Sensor position - 0: TA1_S1, 1: TA1_s2, 3:TA2_s1')
   pylab.ylabel('Bfield (G)')
   pylab.plot([0,30],[brta,brta])

   i = 1
   for value in TA_sensor_r.values():
      pylab.figure(1)
      pylab.plot(i,value[np.argmin(abs(timeta-time))]*Gauss,'o')
      i += 3
   #need to figure out why I need the negative sign!!!!
   fig3 = pylab.figure(2)   
   pylab.title('Ideal Vacuum Fields (Z) vs. TA measurements @ {0:2}ms, shot# {1:5}'.format(time*1000,shotno))
   pylab.xlabel('Sensor position - 0: TA1_S1, 1: TA1_s2, 3:TA2_s1')
   pylab.ylabel('Bfield (G)')
   pylab.plot([0,30],[bzta,bzta])
   pylab.ylim([min(0,2*bzta),max([0,2*bzta])])
   i = 0
   for value in TA_sensor_p.values():
      pylab.figure(2)
      pylab.plot(i,value[np.argmin(abs(timepa-time))]*Gauss,'o')
      i += 1
 
   return ( (brpa,bzpa), (brta,bzta) )

def fields(shotno,(PA_R,PA_Z),(TA_R,TA_Z),time):

   if shotno == 81796:

      (t,I) = SH.Sh(shotno)
      (r,z,j) = magnetic_surfaces.shaping_coil_2_4_2('1_0',1)

   if shotno == 80825 or shotno == 80241 or shotno == 82065:

      (t,I) = Oh.Oh(shotno)
      (r,z,j) = magnetic_surfaces.OH()

   if shotno == 81077 or shotno == 82066:

      (t,I) = Vf.Vf(shotno)
      (r,z,j) = magnetic_surfaces.VF()

   I_time = I[np.argmin(abs(t-time))]

   j_time = I_time*j
   print j_time

   fig3 = pylab.figure(3)
   pylab.plot(r,z,'og')

   BR_PA = np.zeros((len(PA_R)))
   BZ_PA = np.zeros((len(PA_R)))

   BR_TA = 0
   BZ_TA = 0
      
   for i in range(len(r)):
      (BR,BZ) = Bfield.Dipole((PA_R,PA_Z),r[i],z[i],j_time[i])
      BR_PA += BR
      BZ_PA += BZ
      
      (BR,BZ)= Bfield.Dipole((TA_R[0],TA_Z[0]),r[i],z[i],j_time[i])
      print BZ/Gauss
      BR_TA += BR
      BZ_TA += BZ
 
   return((BR_PA,BZ_PA),(BR_TA,BZ_TA))


def positions():
   '''Returns the locations of all HD sensors in R,Z,Phi space. 
See P. 61 of book 6 for explanation of numbers used
'''

   if os.path.isfile('/home/byrne/APS_Work/2011_Shaping_Coil/Circuit_simulations/Eddy_currents/sensor_locations'):
      f = open('/home/byrne/APS_Work/2011_Shaping_Coil/Circuit_simulations/Eddy_currents/sensor_locations','r')
      ((PA_R,PA_Z,PA_phi),(TA_R,TA_Z,TA_phi)) = pickle.load(f)
      f.close()

   else:
      PA_R1 = 35.549*inches
      PA_R2 = 36.220*inches
      PA_r = .1575

      PA_theta1 = (np.arange(-7,7)+.5)*11.63*degrees + 180*degrees
      PA_theta2 = (np.arange(-8,8)+.5)*11.63*degrees

      PA_R_noncirc = PA_R2-.975*inches
      PA_Z_noncirc = PA_r
      
      PA_R = np.append(PA_R2 + PA_r * np.cos(PA_theta2) , PA_R_noncirc)
      PA_R = np.append(PA_R, PA_R1+PA_r*np.cos(PA_theta1))
      PA_R = np.append(PA_R,PA_R_noncirc)

      PA_Z = np.append(PA_r*np.sin(PA_theta2),PA_Z_noncirc)
      PA_Z = np.append(PA_Z,PA_r*np.sin(PA_theta1))
      PA_Z = np.append(PA_Z,-PA_Z_noncirc)

      PA_phi = np.ones(32)*(-41*degrees)
      PA_phi = np.append(PA_phi,PA_phi+180*degrees)

      TA_R = PA_R1*np.ones(30)-PA_r
      TA_Z = -1.5*inches*np.ones(30)
      TA_phi = np.zeros((30))

      for i in range(30):
         TA_sub_array = i//3
         TA_phi[i] = 9*(i%3-1)*degrees+36*TA_sub_array*degrees

      f = open('/home/byrne/APS_Work/2011_Shaping_Coil/Circuit_simulations/Eddy_currents/sensor_locations','w')
      pickle.dump( ((PA_R,PA_Z,PA_phi),(TA_R,TA_Z,TA_phi)), f)
      f.close()

   fig3 = pylab.figure(3)
   ax = fig3.add_subplot(111,aspect = 'equal')
   pylab.plot(TA_R[0],TA_Z[0],'ob')
   return ((PA_R,PA_Z,PA_phi),(TA_R,TA_Z,TA_phi))
