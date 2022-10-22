#%%
# This file integrates a single Resonant object for 10 Myr and declares it as
# TRUE or FALSE regarding whether it is in the 3:2 MMR. It also saves a PDF plot
# of the 10 Myr history of the 3:2 angle.
from matplotlib import pyplot as plt
import aa_utilities as ut
import numpy as np
import pandas as pd
import rebound
import time
t0 = time.time()
THIS_INSTANCE = 1
# Parameters for MMR angle to plot, ie the 3:2 angle for the Plutinos.
big = 3
small = 2
# Settings for integration interval
tmax_yrs = 1e7
tstep_yrs = 100
# Settings for resonance check
settings =  {}
settings['bigmax'] = [3]
settings['smallmax'] = [2]
settings['ordermax'] = [1]
settings['min_stick_length'] = [1000] # how many points in a window
settings['smatol'] = [0.05] # not used
settings['shift_amount'] = [100] # how many time steps to shift overlapping resonance search windows
settings['mmr_tolerance'] = [0.1] # how far away from the nominal period ratio to look for a resonance
settings['libration_threshold'] = [175] # maximum accepted resonance amplitude
settings['points_to_average'] = [10] # used in checking libration within a single window
settings['averaging_proportion'] = [2e-2] # used in recomputing centers and amplitudes of continuous sticks
settings['long_or_small_threshold_pt_count'] = [2000] # how many points makes a stick long or short when recomputing amplitude
settings['small_division_count'] = [6] # how many divisions to split a short stick into
settings['long_division_count'] = [20] # how many divisions to split a long stick into
settings['length_threshold'] = [0.5] # what total fraction of the time an object must librate to count as Resonant
settings['scattering_tolerance_delta_P'] = [0.2] # if period ratio tno/neptune changes by this much, it's Scattering
settings['scattering_tolerance_delta_a'] = [1.5] # if tno semimajor axis changes by this many au, it's Scattering
# Load heliocentric orbital elements of Resonant object.
resonant_file = 'resonant_objects.csv'
df = pd.read_csv(resonant_file)
des = df['resonant_des'][THIS_INSTANCE-1]
ePh = df['ePh'][THIS_INSTANCE-1] # eccentricity, Plutino, heliocentric
qPh_au = df['qPh_au'][THIS_INSTANCE-1] # perihelion distance, Plutino, heliocentric, au
tpPh_jd = df['tpPh_jd'][THIS_INSTANCE-1] # time of perihelion passage, Plutino, heliocentric, Julian date TDB
WPh_deg = df['WPh_deg'][THIS_INSTANCE-1] # longitude of ascending node, Plutino, heliocentric, degrees
wPh_deg = df['wPh_deg'][THIS_INSTANCE-1] # argument of perihelion, Plutino, heliocentric, degrees
iPh_deg = df['iPh_deg'][THIS_INSTANCE-1] # inclination, Plutino, heliocentric, degrees
# Load masses of outer planets.
GMdict = ut.get_GMdict()
# Load heliocentric orbital elements of outer planets.
planet_file = 'planets_for_' + des + '.csv'
df = pd.read_csv(planet_file)
epochP_jd = df['epochP_jd'][0]
eJh = df['eJh'][0] # Jupiter
qJh_au = df['qJh_au'][0]
tpJh_jd = df['tpJh_jd'][0]
WJh_deg = df['WJh_deg'][0]
wJh_deg = df['wJh_deg'][0]
iJh_deg = df['iJh_deg'][0]
eSh = df['eSh'][0] # Saturn
qSh_au = df['qSh_au'][0]
tpSh_jd = df['tpSh_jd'][0]
WSh_deg = df['WSh_deg'][0]
wSh_deg = df['wSh_deg'][0]
iSh_deg = df['iSh_deg'][0]
eUh = df['eUh'][0] # Uranus
qUh_au = df['qUh_au'][0]
tpUh_jd = df['tpUh_jd'][0]
WUh_deg = df['WUh_deg'][0]
wUh_deg = df['wUh_deg'][0]
iUh_deg = df['iUh_deg'][0]
eNh = df['eNh'][0] # Neptune
qNh_au = df['qNh_au'][0]
tpNh_jd = df['tpNh_jd'][0]
WNh_deg = df['WNh_deg'][0]
wNh_deg = df['wNh_deg'][0]
iNh_deg = df['iNh_deg'][0]
# Build simulation
sim = rebound.Simulation()
sim.add(m = 1,hash = '0')
sim.integrator = 'ias15'
epochobj = epochP_jd
# Build simulation, outer planets first.
eobj = eJh # Jupiter
qobj = qJh_au
tpobj = tpJh_jd
Wobj = np.radians(np.mod(WJh_deg,360))
wobj = np.radians(np.mod(wJh_deg,360))
iobj = np.radians(iJh_deg)
aobj = qobj/(1-eobj)
dt = epochobj - tpobj # time since pericenter passage in days
dt = dt * (2*np.pi/365.25) # convert days to yr/(2pi)
n = np.sqrt(1/aobj**3) # radians / (yr/2pi)
Mobj = np.mod(n*dt,2*np.pi) # radians
sim.add(primary=sim.particles[0],m=GMdict['Jupiter'],hash='Jupiter',\
        a=aobj,e=eobj,inc=iobj,omega=wobj,Omega=Wobj,M=Mobj)
eobj = eSh # Saturn
qobj = qSh_au
tpobj = tpSh_jd
Wobj = np.radians(np.mod(WSh_deg,360))
wobj = np.radians(np.mod(wSh_deg,360))
iobj = np.radians(iSh_deg)
aobj = qobj/(1-eobj)
dt = epochobj - tpobj # time since pericenter passage in days
dt = dt * (2*np.pi/365.25) # convert days to yr/(2pi)
n = np.sqrt(1/aobj**3) # radians / (yr/2pi)
Mobj = np.mod(n*dt,2*np.pi) # radians
sim.add(primary=sim.particles[0],m=GMdict['Saturn'],hash='Saturn',\
        a=aobj,e=eobj,inc=iobj,omega=wobj,Omega=Wobj,M=Mobj)
eobj = eUh # Uranus
qobj = qUh_au
tpobj = tpUh_jd
Wobj = np.radians(np.mod(WUh_deg,360))
wobj = np.radians(np.mod(wUh_deg,360))
iobj = np.radians(iUh_deg)
aobj = qobj/(1-eobj)
dt = epochobj - tpobj # time since pericenter passage in days
dt = dt * (2*np.pi/365.25) # convert days to yr/(2pi)
n = np.sqrt(1/aobj**3) # radians / (yr/2pi)
Mobj = np.mod(n*dt,2*np.pi) # radians
sim.add(primary=sim.particles[0],m=GMdict['Uranus'],hash='Uranus',\
        a=aobj,e=eobj,inc=iobj,omega=wobj,Omega=Wobj,M=Mobj)
eobj = eNh # Neptune
qobj = qNh_au
tpobj = tpNh_jd
Wobj = np.radians(np.mod(WNh_deg,360))
wobj = np.radians(np.mod(wNh_deg,360))
iobj = np.radians(iNh_deg)
aobj = qobj/(1-eobj)
dt = epochobj - tpobj # time since pericenter passage in days
dt = dt * (2*np.pi/365.25) # convert days to yr/(2pi)
n = np.sqrt(1/aobj**3) # radians / (yr/2pi)
Mobj = np.mod(n*dt,2*np.pi) # radians
sim.add(primary=sim.particles[0],m=GMdict['Neptune'],hash='Neptune',\
        a=aobj,e=eobj,inc=iobj,omega=wobj,Omega=Wobj,M=Mobj)
# Add the Plutino
eobj = ePh
qobj = qPh_au
tpobj = tpPh_jd
Wobj = np.radians(np.mod(WPh_deg,360))
wobj = np.radians(np.mod(wPh_deg,360))
iobj = np.radians(iPh_deg)
aobj = qobj/(1-eobj)
dt = epochobj - tpobj # time since pericenter passage in days
dt = dt * (2*np.pi/365.25) # convert days to yr/(2pi)
n = np.sqrt(1/aobj**3) # radians / (yr/2pi)
Mobj = np.mod(n*dt,2*np.pi) # radians
sim.add(primary=sim.particles[0],m=0,hash='Plutino',\
        a=aobj,e=eobj,inc=iobj,omega=wobj,Omega=Wobj,M=Mobj)
# Prepare to integrate simulation.
sim.N_active = 5
sim.move_to_com()
# Integrate simulation and return orbels_dict for use in Yu2018 resonance check.
orbels_dict = ut.integrate_single_tno(sim,tmax_yrs,tstep_yrs)
# Check whether object is Resonant (TRUE or FALSE) and which resonance it
# occupies (hopefully big=3, small=2).
resonant,big,small = ut.check_resonance(orbels_dict,settings)
# Plot the identified resonant angle (hopefully 3:2) or, if no resonant angle
# is identified, plot the 3:2 angle anyway.
try:
    bigint = int(float(big))
    smallint = int(float(small))
    l_tno = orbels_dict['l']
    l_Neptune = orbels_dict['l_Neptune']
    pomega_tno = orbels_dict['pomega']
    l_tno = np.array(l_tno)
    l_Neptune = np.array(l_Neptune)
    pomega_tno = np.array(pomega_tno)
    yvar = bigint*l_tno - smallint*l_Neptune - (bigint-smallint)*pomega_tno
    yvar = np.mod(yvar,360)
    ylabel = 'theta_' + big + '_' + small
except:
    big = 3
    small = 2
    l_tno = orbels_dict['l']
    l_Neptune = orbels_dict['l_Neptune']
    pomega_tno = orbels_dict['pomega']
    l_tno = np.array(l_tno)
    l_Neptune = np.array(l_Neptune)
    pomega_tno = np.array(pomega_tno)
    yvar = big*l_tno - small*l_Neptune - (big-small)*pomega_tno
    yvar = np.mod(yvar,360)
    ylabel = 'theta_' + str(big) + '_' + str(small)
xvar = orbels_dict['t']
xvar = np.array(xvar)
xvar = xvar/1e6 # Myr
xvar = xvar[0:-1:10]
yvar = yvar[0:-1:10]
xlabel = 'Myr'
title = 'plot_' + str(resonant) + '_' + str(big) + '_' + str(small) + '_' + des
plt.plot(xvar,yvar,'b.',markersize=0.5)
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.title(title)
plt.savefig(title+'.pdf',format='pdf',transparent=True)
plt.show()
t1 = time.time()
print('elapsed hours',(t1-t0)/3600)
