#%%
# This file classifies all the clones of a single Kuiper belt object.
THIS_INSTANCE = 1
import aa_utilities as ut
import time
import numpy as np
import pandas as pd
import rebound
t0 = time.time()
# Prime the machine learning classifier.
classifier, int_dict, types_dict = ut.get_sv20classifier()
# Get list of objects to classify.
sbdb_file = 'sbdb_reduce_output.csv'
df = pd.read_csv(sbdb_file)
des_list = df['packed_designation'].tolist()
des = des_list[THIS_INSTANCE-1]
# Read heliocentric orbital elements of clones.
clones_input_file = 'clones_' + des + '.csv'
df = pd.read_csv(clones_input_file)
Nclones = df.shape[0]
ePh     = df['ePh'].tolist() # eccentricity, Plutino, heliocentric
qPh_au  = df['qPh_au'].tolist() # perihelion distance, Plutino, heliocentric, au
tpPh_jd = df['tpPh_jd'].tolist() # time of perihelion passage, Plutino, heliocentric, Julian date TDB
WPh_deg = df['WPh_deg'].tolist() # longitude of ascending node, Plutino, heliocentric, degrees
wPh_deg = df['wPh_deg'].tolist() # argument of perihelion, Plutino, heliocentric, degrees
iPh_deg = df['iPh_deg'].tolist() # inclination, Plutino, heliocentric, degrees
# Read heliocentric orbital elements of planets.
planets_file = 'planets_for_' + des + '.csv'
df = pd.read_csv(planets_file)
epochP  = df['epochP_jd'][0] # epoch of orbital elements, Julian date
eJh     = df['eJh'][0] # Jupiter
qJh_au  = df['qJh_au'][0]
tpJh_jd = df['tpJh_jd'][0]
WJh_deg = df['WJh_deg'][0]
wJh_deg = df['wJh_deg'][0]
iJh_deg = df['iJh_deg'][0]
eSh     = df['eSh'][0] # Saturn
qSh_au  = df['qSh_au'][0]
tpSh_jd = df['tpSh_jd'][0]
WSh_deg = df['WSh_deg'][0]
wSh_deg = df['wSh_deg'][0]
iSh_deg = df['iSh_deg'][0]
eUh     = df['eUh'][0] # Uranus
qUh_au  = df['qUh_au'][0]
tpUh_jd = df['tpUh_jd'][0]
WUh_deg = df['WUh_deg'][0]
wUh_deg = df['wUh_deg'][0]
iUh_deg = df['iUh_deg'][0]
eNh     = df['eNh'][0] # Neptune
qNh_au  = df['qNh_au'][0]
tpNh_jd = df['tpNh_jd'][0]
WNh_deg = df['WNh_deg'][0]
wNh_deg = df['wNh_deg'][0]
iNh_deg = df['iNh_deg'][0]
# Get masses of outer planets.
GMdict = ut.get_GMdict()
# Create lists to store probabilities of class membership for each clone.
clone_class_list = []
classical_list = []
resonant_list = []
scattering_list = []
detached_list = []
# Run a Rebound simulation for each clone and classify it.
for iclone in range(Nclones):
    sim = rebound.Simulation()
    sim.add(m = 1,hash = '0')
    sim.integrator = 'ias15'
    epochobj = epochP
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
    # Add the Plutino clone.
    eobj = ePh[iclone]
    qobj = qPh_au[iclone]
    tpobj = tpPh_jd[iclone]
    Wobj = np.radians(np.mod(WPh_deg[iclone],360))
    wobj = np.radians(np.mod(wPh_deg[iclone],360))
    iobj = np.radians(iPh_deg[iclone])
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
    time_outs = np.linspace(0,100E3,101)*2*np.pi # 100 kyr
    data = []
    for i,t in enumerate(time_outs):
        if t>0:
            sim.move_to_com()
            # Integrate to next output.
            sim.integrate(t, exact_finish_time=True)
        orbits = sim.calculate_orbits(primary=sim.calculate_com())
        o = orbits[-1] # take KBO
        # Save t, a, e, i, Omega, omega. Time in data needs to be in years, so divide by 2pi.
        step = np.array([t/2/np.pi, o.a, o.e, np.degrees(o.inc), np.degrees(o.Omega)%360, np.degrees(o.omega)%360])
        # Add step to data.
        if len(data)==0: data = step
        else: data = np.vstack((data,step))
    # Release memory so we don't accidentally keep integrating the same sim with a different clone.
    sim = None
    category,prediction = ut.parsedata_sv20classifier(data,classifier,int_dict)
    clone_class_list.append(category)
    detached_list,resonant_list,scattering_list,classical_list = \
        ut.add_category_lists_sv20classifer(int_dict,prediction,detached_list,\
                        resonant_list,scattering_list,classical_list)
# Make new, separate file of clone class lists.
df = pd.DataFrame()
df['category'] = clone_class_list
df['classical_probability'] = classical_list
df['resonant_probability'] = resonant_list
df['scattering_probability'] = scattering_list
df['detached_probability'] = detached_list
clones_output_file = 'class_lists_' + des + '.csv'
df.to_csv(clones_output_file,index=False)
t1 = time.time()
elapsed_time_hours = (t1-t0)/3600
print('took',elapsed_time_hours,'hrs',THIS_INSTANCE,des)
