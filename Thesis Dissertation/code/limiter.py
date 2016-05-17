#!/usr/bin/python
from __future__ import division
import numpy as np
import pylab
import sys
sys.path.append('/home/byrne/byrne_edited_code')
import scanner

degrees = np.pi/180.
r_flange = 8.75*2.54
R_chamber = 76.16
flange_angle = 18*degrees
    
offset = 22.26
PA_tor_angle = 5*degrees
PA_radius = 15
#PA hemispheres are offset from each other.  Inboard centroid @ 90.3
PA_center = 92#90.294 

def flange_to_RZ(R= None,Z=None):
    '''Converts the flange of a chamber to its projection at the 
    HD Poloidal Array'''
    
    PA_pol_angle = np.linspace(0.5*np.pi, 1.5*np.pi,200)

    psi = np.linspace(0,2*np.pi,200)
    
    
    Z = r_flange*np.sin(psi)
    
    R = np.sqrt( 
        ( R_chamber*np.sin(flange_angle) + 
          r_flange*np.sin(flange_angle)*np.cos(psi)
          )**2 + 
        ( offset + R_chamber*np.cos(flange_angle) + 
          r_flange*np.cos(flange_angle)*np.cos(psi)
          )**2
        )

    x = R*np.cos(flange_angle)+offset
    y = R*np.sin(flange_angle)

    #r = np.sqrt(x**2+y**2)
    

    z_pa = 15*np.sin(PA_pol_angle)
    
    r_pa = PA_center+PA_radius*np.cos(PA_pol_angle) 

    #((Rprime,Zprime),(Rprime_flange,Zprime_flange)) = RZ_to_flange()

    #print Rprime,Zprime
    #xprime_PA = Rprime*np.cos(18*degrees)+offset
    #yprime_PA = Rprime*np.sin(18*degrees)

    #dubl_translate_PA_R = np.sqrt(xprime_PA**2+yprime_PA**2)
    
    fig = pylab.figure()
    ax = fig.add_subplot(111,aspect = 'equal')
    pylab.plot(R,Z)
    pylab.plot(r_pa,z_pa)
    #pylab.plot(dubl_translate_PA_R,Zprime)
    pylab.show()
    return((r_pa,z_pa),(r,Z))

def RZ_to_flange():

    PA_pol_angle = np.linspace(0.5*np.pi,np.pi,100)

    R = PA_radius*np.cos(PA_pol_angle) + PA_center
    Z = PA_radius*np.sin(PA_pol_angle)


    A = offset*np.tan(flange_angle)**2
    B = PA_center+PA_radius*np.cos(PA_pol_angle)
    C = 1 + np.tan(flange_angle)**2
    
    term1 = A/(B*C)

    D = offset*np.tan(flange_angle)

    term2 = ( 1 - (D/B)**2) /C
    
    print term1[-1]
    print term2[-1]

    print (term1 + np.sqrt(term1**2 + term2))[-1]

    phi = np.arccos( term1 + np.sqrt(term1**2 + term2) )
    print phi/degrees

    x = R*np.cos(phi)
    y = R*np.sin(phi)

    x_prime = x - offset
    
    
    R_prime = np.sqrt( x_prime**2 + y**2)
    
    psi = np.arctan2( Z, R_prime )
    return ( (R_prime, Z),
             (R_chamber + r_flange*np.cos(np.linspace(0,2*np.pi)), 
              r_flange*np.sin(np.linspace(0,2*np.pi)) )
             )
    #fig = pylab.figure()
    #ax = fig.add_subplot(111,aspect = 'equal')
    #pylab.plot(R_chamber + r_flange*np.cos(np.linspace(0,2*np.pi)), 
    #           r_flange*np.sin(np.linspace(0,2*np.pi)))

    #pylab.plot(R_prime,Z)
    #pylab.show()

def RZ_to_flange2():

    PA_pol_angle = np.linspace(-0.5*np.pi,1.5*np.pi,1000)

    R = PA_radius*np.cos(PA_pol_angle) + PA_center
    Z = PA_radius*np.sin(PA_pol_angle)


    term1 = offset*np.cos(flange_angle)
    term2 = R**2-offset**2
    
    r = -term1+ np.sqrt(term1**2 + term2)
    fig = pylab.figure()
    ax = fig.add_subplot(111,aspect = 'equal')
    pylab.title('PA (red) projected onto flange (blue). Circular approximation in green')
    pylab.xlabel('R (cm)')
    pylab.ylabel('Z (cm)')

    pylab.plot(r_flange*np.cos(PA_pol_angle)+R_chamber,r_flange*np.sin(PA_pol_angle))
    pylab.plot(r,Z)
    
    r_mod = Z[0]*np.cos(PA_pol_angle)
    pylab.plot(r_mod+r[0],Z)

    (PA_array,shadow) = flange_to_RZ(r,Z)
    (PA_array,shadow_mod) = flange_to_RZ(r_mod+r[0],Z)

    fig2 = pylab.figure()
    ax = fig.add_subplot(111,aspect = 'equal')
    pylab.title('limiter (blue) projected onto PA (red). Circular approximation in green')
    pylab.xlabel('R (cm)')
    pylab.ylabel('Z (cm)')
    pylab.plot(shadow[0],shadow[1])
    pylab.plot(shadow_mod[0],shadow_mod[1])
    pylab.plot(PA_array[0],PA_array[1])

    pylab.figure()
    ax = fig.add_subplot(111,aspect = 'equal')
    pylab.title('Difference in ')
    pylab.xlabel('angle (deg)')
    pylab.ylabel('difference (cm)')
    pylab.plot( PA_pol_angle/degrees,np.sqrt(( (shadow_mod[0]-shadow_mod[0][0])**2 + shadow_mod[1]**2))-
                np.sqrt(( (shadow[0]-shadow[0][0])**2 + shadow[1]**2)))

def Chamber_to_TF():
    '''Draws the cross-section of the vacuum chamber as it appears in the TF bore.  Also calculates the distance between the chamber wall and the TF inner diameter.'''

    psi = np.linspace(0,2*np.pi,600)
    chamber_angle = np.zeros((600))
    r_chamber = 21*2.54/2.
    r_TF_bore = 24.63*2.54/2.
    R_TF_bore = 38.251*2.54

    Z = r_chamber*np.sin(psi)

    a = ( 1 + np.tan(9*degrees)**2 )

    b = ( (2*offset*np.tan(9*degrees)**2) / 
          ( R_chamber + r_chamber*np.cos(psi) )
          )

    c = ( ( offset*np.tan(9*degrees) ) / 
          ( R_chamber + r_chamber*np.cos(psi) )
          )**2 -1
    #quadratic formula gives two roots.  If answers don't make sense,
    #check both roots to be sure you're using the correct one!
    cos_theta = quad_form(a,b,c)
    
    theta = np.arccos(cos_theta)
    
    x = (R_chamber + r_chamber*np.cos(psi))*np.cos(theta) + offset
    y = (R_chamber + r_chamber*np.cos(psi))*np.sin(theta)

    R = np.sqrt(x**2+y**2)
    dist  = np.zeros((150))

    #I only care about the arc from 90-180 degrees:
    for i in range(150):
        #(index,value) = scanner.scan1d(Z[:450]-R[1,:450]*np.tan(psi[150+i]),
        #                               -97.33*np.tan(psi[150+i]))
        #dist[i] = np.sqrt( (r_TF_bore*np.cos(psi[index])+R_TF_bore - 
        #                    R[1,index])**2 +
        #                   (r_TF_bore*np.sin(psi[index]) - 
        #                    Z[index])**2 
        #                   )
        #print (r_TF_bore*np.cos(psi[index])+R_TF_bore - R[1,index])
        #print (r_TF_bore*np.sin(psi[index]) - Z[index])
        angle = np.arctan2(Z[299-i]/(R[299-i]-R_TF_bore))
        dist[i] = np.sqrt( ()**2 + ()**2)

    return (R,Z,dist)

def quad_form(a,b,c):
    
   sign = np.array([-1,1])
   output = np.zeros((2,len(b)))
   output[0,:] = (-b + sign[0]*np.sqrt(b**2 - 4*a*c)) / (2.*a)
   output[1,:] = (-b + sign[1]*np.sqrt(b**2 - 4*a*c)) / (2.*a)
   
   return output


def general_projection((x,y,z),(X,Y,Z)):
    ''' Plugging in a point cloud (x,y,z) of the object in question,
    then the Offset (X,Y,Z), this function will give the cylindrical 
    coordinates in the coordinate system centered at (-X,-Y,-Z) in 
    the old system.'''
    
    R = np.sqrt(x**2+y**2)
    phi = np.arctan2(y,x)

    x -= X
    y -= Y
    z -= Z

    Rproj = np.sqrt(x**2+y**2)
    phiproj = np.arctan2(y,x)

    return ((R,phi,z),(Rproj,phiproj,z))

def point_cloud(flange = None, PA = None,Rc = None,phic = None,zc = None):
    angle = np.linspace(-0.5*np.pi,1.5*np.pi,400)
    if flange:
        x = ( R_chamber*np.sin(flange_angle) +
              r_flange*np.sin(flange_angle)*np.cos(angle) )

        y = ( R_chamber*np.cos(flange_angle) +
              r_flange*np.cos(flange_angle)*np.cos(angle) )

        z = r_flange*np.sin(angle)
        
        offset = -22.26

    if PA:
        x = ( PA_center*np.sin(PA_tor_angle) +
              PA_radius*np.sin(PA_tor_angle)*np.cos(angle[200:]) )

        y = ( PA_center*np.cos(PA_tor_angle) +
              PA_radius*np.cos(PA_tor_angle)*np.cos(angle[200:]) )

        z = PA_radius*np.sin(angle[200:])

        offset = 22.26
        
    if  phic :
        x = Rc*np.sin(flange_angle)
        y = Rc*np.cos(flange_angle)
        z = zc

        offset = -22.26
    return((x,y,z),(0,offset,0))

def RZ_to_flange3():

    #PA_pol_angle = np.linspace(0.5*np.pi,1.5*np.pi,10000)
    PA_pol_angle = np.linspace(-0.5*np.pi,0.5*np.pi,10000)
    R = PA_radius*np.cos(PA_pol_angle) + PA_center
    Z = PA_radius*np.sin(PA_pol_angle)

    ShellflangeR = np.linspace(-2.5*2.54,0,10)
    ShellflangeZ = np.linspace(15,15,10)

    R = np.append(R,ShellflangeR)
    Z = np.append(Z,ShellflangeZ)
    #R = PA_radius*np.cos(180*degrees) + PA_center
    #Z = PA_radius*np.sin(180*degrees)
    a = 1 + np.tan(18*degrees)**2
    b = 2*offset
    c = offset**2-R**2

    y_flange = (-b + np.sqrt( b**2 - 4*c*a))/(2*a)
    x_flange = y_flange*np.tan(18*degrees)
    r = np.sqrt(x_flange**2+y_flange**2)

    #R_shadow = np.sqrt(x_flange**2+(y_flange+offset)**2)
    
    #fig = pylab.figure()
    #ax = fig.add_subplot(111,aspect = 'equal')
    #pylab.title(
    #    'PA (red) projected onto flange (blue). Circular approximation in green'
    #    )

    #pylab.xlabel('R (cm)')
    #pylab.ylabel('Z (cm)')
    #pylab.plot(r_flange*np.cos(PA_pol_angle)+R_chamber,
    #           r_flange*np.sin(PA_pol_angle))
    #pylab.plot(r,Z)
    
    #a  = (max(Z) - min(Z))/2.
    #b  = (r[0] - r[4999])
    #dz = (max(Z) + min(Z))/2.
    #dr = (r[0] + r[9999])/2.

    #print r[0], r[9999]
    #print Z[0], Z[9999]
    #print r[0], r[4999]
    #print Z[0], Z[4999]
    #print a, dz, b, dr

    #flange_center = R[0]
    #limiter_r = Z[0]*np.cos(PA_pol_angle)+flange_center
    #limiter_z = Z[0]*np.sin(PA_pol_angle)
    
    #shadow_z = limiter_z
    #shadow_x = limiter_r*np.sin(18*degrees)
    #shadow_y = limiter_r*np.cos(18*degrees)

    #shadow_R = np.sqrt(shadow_x**2+(shadow_y)**2)
    #shadow_Z = shadow_z

    #pylab.plot((b)*np.cos(PA_pol_angle)+dr,(a)*np.sin(PA_pol_angle))

    #fig = pylab.figure()
    #ax = fig.add_subplot(111,aspect = 'equal')
    #pylab.plot(shadow_R,shadow_Z)
    #pylab.plot(PA_radius*np.cos(PA_pol_angle)+PA_center,
    #           PA_radius*np.sin(PA_pol_angle))
    #pylab.plot(np.sqrt(x_flange**2+(y_flange+offset)**2),Z)
    #fig = pylab.figure()
    #pylab.plot(np.sqrt(x_flange**2+(y_flange+offset)**2)-PA_radius*np.cos(PA_pol_angle)+PA_center)
    
    return r,Z,R
