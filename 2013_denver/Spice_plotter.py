#!/usr/bin/python
import pylab
import sys
sys.path.append('/home/byrne/APS_Work/2011_Shaping_Coil/Circuit_simulations')
import SH
import numpy as np

def main():
    currentstr = '3.5kA'

    spicefile = open("/home/byrne/APS_Work/LateX/2013_denver/new_plots/Final_rev_banks_APS_2013.txt")
    data  = spicefile.readlines()
    spicefile.close()

    
    N = len(data)
    print N
    print len(range(1,N))
    time = np.zeros((N-1), dtype=np.float)
    I_coil = np.zeros((N-1), dtype=np.float)
    print len(time)
    for i in range(1,N-1):
        entries = data[i].split()
        time[i] = float(entries[0])
        I_coil[i] = float(entries[1])
    
    (t_real,I_coil_real) = SH.Sh(81767)
    pylab.plot(time*1000,I_coil*.001, label = 'SPICE Simulation')
    pylab.plot(t_real*1000-2,I_coil_real*.001,label = 'Measured Current')
    pylab.title('Simulation Vs Actual Discharge - Shot 81767')
    pylab.xlabel('Time (ms)')
    pylab.ylabel('Coil Current (kA)')
    pylab.xlim([0,10])
    pylab.ylim([0,.001*1.1*np.max([np.max(I_coil),np.max(I_coil_real[np.where(t_real>.002)])])])
    pylab.legend()
    #pylab.show()
    pylab.savefig('shaping_current_sim_vs_real_81797')
