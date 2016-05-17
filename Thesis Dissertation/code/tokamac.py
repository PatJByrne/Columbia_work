#!/usr/bin/env python
# -*- coding: UTF-8 -*-
'''
Copyright (C) Nikolaus Rath <Nikolaus@rath.org>

This program can be distributed under the terms of the GNU LGPL.
'''

from __future__ import division, print_function, absolute_import

import logging
import subprocess
import os 
import h5py
import shutil
import re
import numpy as np
from collections import namedtuple
import numbers

log = logging.getLogger()

# Define measurement classes
def _make_measurement(name, mesType, fields):
    mes = namedtuple(name, fields + [ 'value', 'stddev'] )
    mes.mType = mesType
    return mes
PlasmaCurrent = _make_measurement('PlasmaCurrent', 9, [])
CoilCurrent = _make_measurement('CoilCurrent', 8, ['coil'])
PoloidalField = _make_measurement('PoloidalField', 1, [ 'angle', 'r', 'z' ])
Pressure = _make_measurement('Pressure', 2, [ 'r', 'z' ])
CosRogowski = _make_measurement('CosRogowski', 7, [ 'r_major', 'r_minor' ])
PoloidalFlux = _make_measurement('PoloidalFlux', 5, [ 'r', 'z' ])
      
class Config(object):
    '''Holds TokaMac configuration data'''

    def __init__(self):
        self.PsiGrid = None
        self.Plasma = None
        self.limiters = None
        self.separatrix_areas = None
        self.coils = None
        self.measurements = None
        self.MaxIterations = 30


def get_pressure(workdir):
    f = h5py.File(os.path.join(workdir, 'TokOut_Plasma.h5'), 'r')
    
    pressure = f['Data-Set-6']
    Y = f[pressure.attrs['DIMENSION_NAMELIST'][0]]
    X = f[pressure.attrs['DIMENSION_NAMELIST'][1]]  
    
    return (X, Y, pressure)
    
def clean(workdir):
    '''Remove all output files created by TokaMac'''
    
    for name in ['TokOut_Plasma.h5', 'TokOut_BndMomts.out',
                 'TokLog.out', 'TokOut_Conductors.out',
                 'TokOut_DCON.dat', 'TokOut_EQGRUM.stb', 'TokOut_FluxPr.out',
                 'TokOut_InValues.out', 'TokOut_Meas.out',
                 'TokOut_Plasma.out', 'TokOut_PsiGrid.h5', 'TokOut_PsiGrid.out',
                 'TokOut_SVDFit.out', 'TokOut_PsiGrid.HDF', 'TokOut_Plasma.HDF',
                 'HGreens.bin', 'MGreen.bin' ]:
#TokIn.dat *was* included in name, but was removed, because it is obviously
#not output, and it is unclear to the modifier (P.Byrne) how it is regenerated
             
        path = os.path.join(workdir, name)
        if os.path.exists(path):
            os.unlink(name)
    
def get_current(workdir):
    f = h5py.File(os.path.join(workdir, 'TokOut_PsiGrid.h5'), 'r')
    
    Current = f['Data-Set-2']
    Y = f[Current.attrs['DIMENSION_NAMELIST'][0]]
    X = f[Current.attrs['DIMENSION_NAMELIST'][1]]  
    
    return (X, Y, Current)
        
def get_flux(workdir):
    f = h5py.File(os.path.join(workdir, 'TokOut_PsiGrid.h5'), 'r')
    
    Psi = f['Data-Set-3']
    Z = f[Psi.attrs['DIMENSION_NAMELIST'][0]]
    X = f[Psi.attrs['DIMENSION_NAMELIST'][1]]  
    
    return (X, Z, Psi)

              
def rcopy(src, dest):
    '''Copy all files in `src` to `dest`'''
    
    for name in os.listdir(src):
        shutil.copy2(os.path.join(src, name), dest)
                  
def run_tokamac(workdir):
    '''Run TokaMac in `workdir`'''
    
    log.debug('Running tokamac')
#changed 1/13/15, P.Byrne.  location of TokaMac executable seems to have moved
#    subprocess.check_call(['/opt/hbt/TokaMac/bin/TokaMac'], stderr=subprocess.STDOUT, cwd=workdir,
#                          stdout=open(os.path.join(workdir, 'TokaMac.log'), 'w')) 
                          
    subprocess.check_call(['/opt/hbt/TokaMac/TokaMac'], 
                          stderr=subprocess.STDOUT, cwd=workdir, 
                          stdout=open(os.path.join(workdir, 'TokaMac.log'), 
                                      'w'))       
    log.debug('Running h4toh5')
    for i in ('TokOut_Plasma.h5', 'TokOut_PsiGrid.h5'):
        if os.path.exists(os.path.join(workdir, i)):
            os.unlink(os.path.join(workdir, i))
    subprocess.check_call(['h4toh5', 'TokOut_Plasma.HDF'],
                          cwd=workdir)
    os.unlink(os.path.join(workdir, 'TokOut_Plasma.HDF'))
    subprocess.check_call(['h4toh5', 'TokOut_PsiGrid.HDF'],
                          cwd=workdir)
    os.unlink(os.path.join(workdir, 'TokOut_PsiGrid.HDF'))


# Attributes defined outside __init__
#pylint: disable=W0201
    
class RunResults(object):
    '''Stores the results of a TokaMac run'''
    
    

    __slots__ = [ 'beta', 'beta_p', 'beta_normal', 'separatrix', 'li',
                 'measurements', 'chi1', 'chi2', 'q_0', 'F_profile',
                 'q_edge', 'q_profile', 'psi_profile', 'P_profile',
                 'r_major', 'mag_axis_r', 'mag_axis_z', 'Q_prob' ]
    
    def __init__(self):
        for name in self.__slots__:
            setattr(self, name, None)

        
    
NUMBER = r'-?[0-9]+(?:.[0-9]+)?(?:[eE]-?[0-9]+)?'    
def get_run_results(workdir):
    
    res = RunResults()
    
    with open(os.path.join(workdir, 'TokOut_Plasma.out'), 'r') as fd:
        for line in fd:
            hit = re.search(r'\sMagnetic axis \(x, z\) = \((%s), (%s)\)'
                              % (NUMBER, NUMBER), line)
            if hit:
                res.mag_axis_r = float(hit.group(1))
                res.mag_axis_z = float(hit.group(2))
                      
            hit = re.search(r'\sCurrent centroid R = (%s) ' % NUMBER, line)
            if hit:
                res.r_major = float(hit.group(1))
                
            hit = re.search(r'\sbeta poloidal = (%s), '
                             'li = (%s)' % (NUMBER, NUMBER), line)
            if hit:
                res.beta_p = float(hit.group(1))
                res.li = float(hit.group(2))
                
            hit = re.search(r'\sVolume-averaged beta = (%s). '
                             'BetaN = (%s)' % (NUMBER, NUMBER), line)
            if hit:
                res.beta = float(hit.group(1))
                res.beta_normal = float(hit.group(2))
    
            hit = re.search(r'\sSafety factor: q0 is (%s); '
                            'q at the %s flux surface is (%s)'
                            % (NUMBER, NUMBER, NUMBER), line)
            if hit:
                res.q_0 = float(hit.group(1))
                res.q_edge = float(hit.group(2))
                
    assert res.li is not None
    assert res.beta_normal is not None
    assert res.beta is not None
    assert res.beta_p is not None
    assert res.q_0 is not None
    assert res.q_edge is not None
    assert res.r_major is not None
    
    with open(os.path.join(workdir, 'TokOut_FluxPr.out'), 'r') as fd:
        q_idx = 3
        psi_idx = 2
        P_idx = 4
        F_idx = 6
        for _ in range(6):
            fd.next()
            
        fields = fd.next().split()
        assert fields[q_idx] == 'q'
        assert fields[psi_idx] == 'Psi'
        assert fields[P_idx] == 'P'
        assert fields[F_idx] == 'F'
        
        qs = []
        psis = []
        Ps = []
        Fs = []
        for line in fd:
            fields = line.split()
            if not fields:
                break
            qs.append(float(fields[q_idx]))
            psis.append(float(fields[psi_idx]))
            Ps.append(float(fields[P_idx]))
            Fs.append(float(fields[F_idx]))
            
        res.q_profile = np.asarray(qs)
        res.psi_profile = np.asarray(psis)
        res.P_profile = np.asarray(Ps)
        res.F_profile = np.asarray(Fs)
            
    res.measurements = dict()
    with open(os.path.join(workdir, 'TokOut_Meas.out'), 'r') as fh:
        for line in fh:
            hit = re.search(r'\sThe chi-square probability, Q = (%s)\.'
                              % NUMBER, line)
            if hit:
                res.Q_prob = float(hit.group(1))
                      
            fields = line.split()
            if fields and fields[0] == 'Name':
                break
        else:
            assert False     
            
        assert res.Q_prob is not None
        for line in fh:
            fields = line.split()
            if fields[0] == '---------':
                break
            res.measurements[fields[0]] = float(fields[4])
  
  
    with open(os.path.join(workdir, 'TokOut_SVDFit.out'), 'r') as fh:
        for _ in range(3):
            fh.next()
            
        line = fh.next().strip()
        match = re.match(r'^Chisq1 = (%s), Chisq2 = (%s)' % (NUMBER, NUMBER), line)
        
    if not match:
        raise RuntimeError('Cannot parse SVDFit.out')
        
    res.chi1 = float(match.group(1))
    res.chi2 = float(match.group(2))
                 
    res.separatrix = get_separatrix(workdir)
    
    return res
    

def get_separatrix(workdir):
    '''Get separatrix coordinates from TokaMac run in `workdir`
    
    Returns list if (X, Z, Psi) coordinates. If there is no separatrix,
    or if the plasma is limited, return the last closed flux surface.
    '''
    
    fh = open(os.path.join(workdir, 'TokLog.out'))
    
    while True:
        line = fh.readline()
        if not line:
            break
        if line.strip().startswith('[Free iteration'):
            pos = fh.tell()
    fh.seek(pos)
    
    sep = list()
    lim = list()
    while True:
        line = fh.readline().strip()
        if not line:
            break
        if line.startswith('[Limited plasma,'):
            return []

        match = re.match(r'\[Sep  found at \(X = (-?[0-9.]+), Z = (-?[0-9.]+)\) Psi = (-?[0-9.]+)e?(-?[0-9.]+)?',line)
        if match:
            x = float(match.group(1))
            z = float(match.group(2))
            psi = float(match.group(3))
            if match.group(4) != None:
                psi *= 10**(float(match.group(4)))
            sep.append((x, z, psi))
        if line.strip()[:9] =='[PsiLim =':
            lim.append(float(line.strip().split()[2][:-1]))
    return (sep,lim)

        
def write_tokamac_config(config, workdir, restart=False):
    log.debug('Writing TokaMac configuration')
    
    PsiGrid = config.PsiGrid
    limiters = config.limiters
    measurements = config.measurements
    coils = config.coils
    Plasma = config.Plasma
    separatrix_areas = config.separatrix_areas
    
    CodeControl = {
        'MaxIterFixed': 0,
        'MaxIterFree': config.MaxIterations,
        'NumShells': 0,
        'Oname': 'TokOut',
        'RSname': 'Restart.bin',
        'MGname': 'MGreen.bin',
        'LHname': 'HGreens.bin',
        'SMname': 'SInduct.bin',
        'SGname': 'SGreen.bin',
        }
    
    if restart:
        CodeControl['RestartStatus'] = 1
        CodeControl['SInductStatus'] = 1
        CodeControl['LHGreenStatus'] = 1
        CodeControl['MGreenStatus'] = 1
        CodeControl['SGreenStatus'] = 1
    else:
        CodeControl['RestartStatus'] = 0
        CodeControl['SInductStatus'] = 0
        CodeControl['LHGreenStatus'] = 0    
        CodeControl['MGreenStatus'] = 0
        CodeControl['SGreenStatus'] = 0
                
    PsiGrid['Symmetric'] = 0
    CodeControl['NumLimiters'] = len(limiters)
    CodeControl['NumMeasures'] = len(measurements)
    CodeControl['NumCoils'] = len(coils)
    CodeControl['NumSeps'] = len(separatrix_areas)  
    
    # Determine Plasma Current
    for mes in measurements.itervalues():
        if isinstance(mes, PlasmaCurrent):
            Plasma['Ip0'] = mes.value
            break
    else:
        raise RuntimeError('No plasma current defined') 
     
    fh = open(os.path.join(workdir, 'TokIn.dat'), 'w')
    
    fh.write('K_PsiGrid\n')
    for (key, val) in PsiGrid.iteritems():
        fh.write('    %s = %s\n' % (key, val))
        
    fh.write('K_CodeControl\n')
    for (key, val) in CodeControl.iteritems():
        fh.write('    %s = %s\n' % (key, val))
        
    fh.write('K_Plasma\n')
    for (key, val) in Plasma.iteritems():
        fh.write('    %s = %s\n' % (key, val))
        
    for tuple_ in limiters:
        fh.write('K_Limiter\n')
        fh.write('    Name = %s\n'
                 '    X1 = %.5g\n'
                 '    Z1 = %.5g\n'
                 '    X2 = %.5g\n'
                 '    Z2 = %.5g\n'
                 '    Enabled = 1\n' % tuple_)
        
    for (i, tuple_) in enumerate(separatrix_areas):
        fh.write('K_Separatrix\n')
        fh.write('    Name = Separatrix_%d\n'
                 '    X1 = %.5g\n'
                 '    Z1 = %.5g\n'
                 '    X2 = %.5g\n'
                 '    Z2 = %5.g\n'
                 '    Enabled = 1\n' % ((i,) + tuple_))
                
    coilCurrents = dict()
    for mes in measurements.itervalues():
        if isinstance(mes, CoilCurrent):
            coilCurrents[mes.coil] = mes.value
    
    coilNo = dict()
    for (no, (name, subcoils)) in enumerate(coils.iteritems()):
        if name not in coilCurrents:
            log.info('No current measurement defined for coil %s' % name)
            initial = 0
        else:
            initial = coilCurrents[name]
            
        coilNo[name] = no+1
        
        fh.write('K_Coil\n')
        fh.write('    Name = %s\n' % name)
        fh.write('    Enabled = 1\n')
        fh.write('    InitialCurrent = %s\n' % initial)
        
        for (no, tuple_) in enumerate(subcoils):
            fh.write('K_SubCoil\n')
            fh.write('    Name = %s_%s\n' % (name, no)) 
            fh.write('    X = %.5g\n'
                     '    Z = %.5g\n'
                     '    CurrentFraction = %.5g\n' % tuple(tuple_))
                    
    for (name, mes) in sorted(measurements.iteritems(), key=lambda x: x[0]):
        fh.write('K_Measure\n')
        fh.write('    mType = %d\n' % mes.mType)
        fh.write('    Name = %s\n' % name)
        
        # Access to protected member
        #pylint: disable=W0212
        fields = mes._asdict()
        if isinstance(mes, CosRogowski):
            fields['CircleType'] = 1
            fields['Z'] = 0.0
            fields['X'] = mes.r_major
            fields['Radius'] = mes.r_minor
            fields['Number'] = 100
            del fields['r_major']
            del fields['r_minor']
            
        elif isinstance(mes, CoilCurrent):
            fields['CoilNum'] = coilNo[fields['coil']]
            del fields['coil']
            
        elif isinstance(mes, (PoloidalField, Pressure)):
            fields['x'] = fields['r']
            del fields['r']
            
        for (key, val) in sorted(fields.iteritems(), key=lambda x: x[0]):
            if isinstance(val, numbers.Real):
                val = '%.5g' % val
            fh.write('    %s = %s\n' % (key, val))
