import pandas as pd
import numpy as np
import os
import aa_utilities as ut
from astroquery.jplhorizons import Horizons
import matplotlib.pyplot as plt
import random
from scipy import special
from scipy import stats
import scipy.optimize as optimization
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
# First, delete any .out files (High Performance Computing output logs).
thisdir = os.getcwd()
files_in_dir = os.listdir(thisdir)
for this_file in files_in_dir:
    if this_file.endswith('.out'):
        os.remove(os.path.join(thisdir,this_file))
# Get epoch at which to retrieve barycentric orbital elements.
epoch_file = 'epoch_settings.csv'
df = pd.read_csv(epoch_file)
year = df['year'][0]
month = df['month'][0]
day = df['day'][0]
hour = df['hour'][0]
minute = df['minute'][0]
second = df['second'][0]
center = '500@0' # solar system barycenter
path = os.getcwd()
JD = ut.JDfun(year,month,day,hour,minute,second)
# Get list of Resonant objects to process.
res_files = 'resonant_objects.csv'
df = pd.read_csv(res_files)
Nres = df.shape[0]
print('Number of Resonant objects = ',Nres)
# Initialize list of Plutinos and their barycentric orbital elements.
plut_des = []
aPb_au = [] # semimajor axis, Plutino, barycentric, au
ePb = [] # eccentricity, Plutino, barycentric
iPb_deg = [] # inclination, Plutino, barycentric, degrees
WPb_deg = [] # longitude of ascending node, Plutino, barycentric, degrees
wPb_deg = [] # argument of perihelion, Plutino, barycentric, degrees
MPb_deg = [] # mean anomaly, Plutino, barycentric, degrees
for ires in range(Nres):
    des = df['resonant_des'][ires]
    plot_file = 'plot_True_3_2_' + des + '.pdf'
    if os.path.isfile(os.path.join(path,plot_file)):
        plut_des.append(des)
        unpacked = ut.unpack(des)
        if des == 'D4340':
            unpacked = '9' # We want the Pluto-Charon barycenter, not the Pluto body center.
        obj = Horizons(id=unpacked,location=center,epochs=str(JD))
        el = obj.elements()
        aPb_au.append(float(el['a'])) # au
        ePb.append(float(el['e']))
        iPb_deg.append(float(el['incl'])) # degrees
        WPb_deg.append(np.mod(float(el['Omega']),360))
        wPb_deg.append(np.mod(float(el['w']),360))
        MPb_deg.append(np.mod(float(el['M']),360))
# # make new dataframe with reordered columns
df2 = pd.DataFrame()
df2['Packed MPC designation'] = plut_des
df2['Semimajor axis au barycentric'] = aPb_au
df2['Eccentricity barycentric'] = ePb
df2['Inclination ecliptic J2000 barycentric degrees'] = iPb_deg
df2['Longitude of ascending node ecliptic J2000 barycentric degrees'] = WPb_deg
df2['Argument of perihelion ecliptic J2000 barycentric degrees'] = wPb_deg
df2['Mean anomaly ecliptic J2000 barycentric degrees'] = MPb_deg
df2['Epoch JD'] = [JD for i in range(len(aPb_au))]
df2['Barycentric element source'] = ['JPL Horizons via Astroquery' for i in range(len(aPb_au))]
df2.to_csv('plutinos_for_mnras.csv',index=False)
Nplut = df2.shape[0]
print('Number of Plutinos = ',Nplut)
# Prepare plots and figures.
name = '8' # neptune barycenter
center = '500@0' # solar system barycenter
obj = Horizons(id=name,location=center,epochs=JD)
el = obj.elements()
i_neptune = np.radians(float(el['incl']))
W_neptune = np.radians(float(el['Omega']))
# equation 1 in the paper
hx_neptune = np.sin(i_neptune)*np.sin(W_neptune)
hy_neptune = -np.sin(i_neptune)*np.cos(W_neptune)
hz_neptune = np.cos(i_neptune)
a_array = np.array(aPb_au)
e_array = np.array(ePb)
i_array = np.array(iPb_deg)
w_array = np.array(wPb_deg)
W_array = np.array(WPb_deg)
M_array = np.array(MPb_deg)
w_array = np.mod(w_array,360)
W_array = np.mod(W_array,360)
M_array = np.mod(M_array,360)
i_array = np.radians(i_array)
w_array = np.radians(w_array)
W_array = np.radians(W_array)
M_array = np.radians(M_array)
# equation 1 in paper
hx_array = np.sin(i_array)*np.sin(W_array)
hy_array = -np.sin(i_array)*np.cos(W_array)
hz_array = np.cos(i_array)
Nobj = len(hz_array)
# want to tell apart plutinos on the north hemisphere and plutinos on the south
# hemisphere in case it makes fig 1 and fig 2 in the paper a little easier to read
hx_array_north = []
hx_array_south = []
hy_array_north = []
hy_array_south = []
hz_array_north = []
hz_array_south = []
for iobj in range(Nobj):
    if hz_array[iobj] >= 0:
        hx_array_north.append(hx_array[iobj])
        hy_array_north.append(hy_array[iobj])
        hz_array_north.append(hz_array[iobj])
    else:
        hx_array_south.append(hx_array[iobj])
        hy_array_south.append(hy_array[iobj])
        hz_array_south.append(hz_array[iobj])
hx_array_north = np.array(hx_array_north)
hy_array_north = np.array(hy_array_north)
hz_array_north = np.array(hz_array_north)
hx_array_south = np.array(hx_array_south)
hy_array_south = np.array(hy_array_south)
hz_array_south = np.array(hz_array_south)
# equation 9 in paper
hxbar = np.mean(hx_array)
hybar = np.mean(hy_array)
hzbar = np.mean(hz_array)
# equation 10 in paper
R = np.sqrt(hxbar**2+hybar**2+hzbar**2)
# equation 11 in paper
muhat_hx = hxbar/R
muhat_hy = hybar/R
muhat_hz = hzbar/R
print('[muhat_hx,muhat_hy,muhat_hz] =',[muhat_hx,muhat_hy,muhat_hz])
# equation A1 in paper
i_mu = np.sqrt(muhat_hx**2+muhat_hy**2)
W_mu = np.mod(np.arctan2(muhat_hx/np.sin(i_mu),-muhat_hy/np.sin(i_mu)),2*np.pi)
print('i_mu =',np.degrees(i_mu),'deg')
print('W_mu =',np.degrees(W_mu),'deg')
# equation 12 in paper
kappa = R*(3-R**2)/(1-R**2)
diff = 1
tol = 1e-8
# equation 13 in paper
while abs(diff) > tol:
    top = special.iv(3/2, kappa)
    bottom = special.iv(1/2, kappa)
    # equation 14 in paper
    A3 = top/bottom
    diff =(A3-R)/(1-A3**2-2/kappa*A3)
    kappa = kappa - diff
print('kappa =',kappa)
print('sigma =',kappa**-0.5,'=',np.degrees(kappa**-0.5),'deg')
sigma = 1/np.sqrt(kappa)
del diff,tol,top,bottom,A3
i_invariable = np.radians(1.578694)
W_invariable = np.radians(107.582222)
# equation 1 in paper
hx_invariable = np.sin(i_invariable)*np.sin(W_invariable)
hy_invariable = -np.sin(i_invariable)*np.cos(W_invariable)
hz_invariable = np.cos(i_invariable)
hx_ecliptic = 0
hy_ecliptic = 0
hz_ecliptic = 1
polevec_invariable = np.array([hx_invariable,hy_invariable,hz_invariable])
polevec_ecliptic = np.array([hx_ecliptic,hy_ecliptic,hz_ecliptic])
polevec_vmf = np.array([muhat_hx,muhat_hy,muhat_hz])
# equation A3 in paper
a = np.cross(polevec_vmf,polevec_ecliptic)
a = a / np.linalg.norm(a)
# equation 18 in paper
theta_rot = np.arccos(np.dot(polevec_ecliptic,polevec_vmf))
# check that the vmf pole is rotated onto the ecliptic pole
rodrigues_check = ut.rodrigues_rotation(polevec_vmf, a, theta_rot)
# equation A5 in the paper
khat = ut.rodrigues_rotation(np.array([0,0,1]),-a,theta_rot)
ihat = ut.rodrigues_rotation(np.array([1,0,0]),-a,theta_rot)
jhat = np.cross(khat,ihat)
print('khat__rel =',khat)
print('ihat__rel =',ihat)
# rotate the plutino poles around so their meanvec is the ecliptic, then find relative Omega
hx_array_post = []
hy_array_post = []
hz_array_post = []
i_array_post = []
W_array_post = []
for iobj in range(Nobj):
    vold = np.array([hx_array[iobj],hy_array[iobj],hz_array[iobj]])
    vnew = ut.rodrigues_rotation(vold,a,theta_rot)
    hx_array_post.append(vnew[0])
    hy_array_post.append(vnew[1])
    hz_array_post.append(vnew[2])
    # equation A4 in the paper
    i_post = np.arcsin(np.sqrt(vnew[0]**2+vnew[1]**2))
    W_post = np.mod(np.arctan2(vnew[0]/np.sin(i_post),-vnew[1]/np.sin(i_post)),2*np.pi)
    i_array_post.append(i_post)
    W_array_post.append(W_post)
hx_array_post = np.array(hx_array_post)
hy_array_post = np.array(hy_array_post)
hz_array_post = np.array(hz_array_post)
i_array_post = np.array(i_array_post)
W_array_post = np.array(W_array_post)
# equation 16 in the paper
# method 1 (ie total,term,d,delta,theta) should be the same as method 2 (total2,term2,d2,delta2,theta2)
total  = 0
total2 = 0
for iobj in range(Nobj):
    term  = hz_array_post[iobj]
    term2 = muhat_hx*hx_array[iobj] + muhat_hy*hy_array[iobj] + muhat_hz*hz_array[iobj]
    total  = total  + term
    total2 = total2 + term2
d  = 1 - total /Nobj
d2 = 1 - total2/Nobj
# equation 15 in the paper
delta  = np.sqrt(d /(Nobj*R**2))
delta2 = np.sqrt(d2/(Nobj*R**2))
# equation 17 in the paper
alpha = 0.05 # 95% confidence region
theta  = np.arcsin(delta *np.sqrt(-np.log(alpha)))
print('theta =',np.degrees(theta),'deg')
theta2 = np.arcsin(delta2*np.sqrt(-np.log(alpha)))
del alpha
# define points on the 95% confidence circle centered on the origin
Npts_circle = 100
circle_clock_angles = np.linspace(start=0,stop=2*np.pi,num=Npts_circle,endpoint=True)
radius_circle = np.sin(theta)
hx_circle_pre = np.cos(circle_clock_angles)*radius_circle
hy_circle_pre = np.sin(circle_clock_angles)*radius_circle
hz_circle_pre = np.cos(theta)*np.ones(Npts_circle)
# rotate the confidence circle around to center on the VMF mean pole
hx_circle_post = []
hy_circle_post = []
hz_circle_post = []
for ipt in range(Npts_circle):
    vold = np.array([hx_circle_pre[ipt],hy_circle_pre[ipt],hz_circle_pre[ipt]])
    vnew = ut.rodrigues_rotation(vold, -a, theta_rot)
    hx_circle_post.append(vnew[0])
    hy_circle_post.append(vnew[1])
    hz_circle_post.append(vnew[2])
hx_circle_post = np.array(hx_circle_post)
hy_circle_post = np.array(hy_circle_post)
hz_circle_post = np.array(hz_circle_post)
#%% 3d rendering of orbit poles on the unit sphere
u, v = np.mgrid[0:np.radians(360):40j, 0:np.radians(90):40j]
x = np.cos(u) * np.sin(v)
y = np.sin(u) * np.sin(v)
z = np.cos(v)
scale = 1
x = x * scale
y = y * scale
z = z * scale
fig = plt.figure(figsize=(2.5,2.5))
plt.rcParams['font.size'] = 8
ax = fig.add_subplot(111,projection='3d')
ax.plot_surface(x,y,z,rstride=1,cstride=1,color='whitesmoke',alpha=0.1,edgecolor='gray',linewidth=0.25) # unit sphere
ax.scatter(hx_array_north,hy_array_north,hz_array_north,color='tomato',s=2) # plutino poles in northern hemisphere
ax.scatter(hx_array_south,hy_array_south,hz_array_south,color='cadetblue',s=2) # plutino poles in southern hemisphere
ax.set_box_aspect((1,1,0.5))
ax.set_xticks(ticks=[-0.5,0.5], minor=False)
ax.set_yticks(ticks=[-0.5,0,0.5], minor=False)
ax.set_zticks(ticks=[0,0.5,1], minor=False)
ax.view_init(60, 45)
plt.tight_layout()
ax.set_xlabel('$h_x$ = sin(i)sin(Ω)')
ax.set_ylabel('$h_y$ = -sin(i)cos(Ω)')
ax.set_zlabel('$h_z$ = cos(i)')
titlestr = 'fig1_' + str(Nplut)
plt.tight_layout()
plt.savefig(titlestr + '.pdf',dpi=300,bbox_inches='tight',pad_inches=0)
plt.show()
#%% plot orbit poles in (hx,hy) plane
fig = plt.figure()
ax = fig.add_subplot(121)
th = np.linspace(start=0,stop=2*np.pi,num=100,endpoint=True)
costh = np.cos(th)
sinth=  np.sin(th)
ax.plot(np.cos(np.radians(90))*costh,np.cos(np.radians(90))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(80))*costh,np.cos(np.radians(80))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(70))*costh,np.cos(np.radians(70))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(60))*costh,np.cos(np.radians(60))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(50))*costh,np.cos(np.radians(50))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(40))*costh,np.cos(np.radians(40))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(30))*costh,np.cos(np.radians(30))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(20))*costh,np.cos(np.radians(20))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(10))*costh,np.cos(np.radians(10))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(0))*costh,np.cos(np.radians(0))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(0*30))],[0,np.sin(np.radians(0*30))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(1*30))],[0,np.sin(np.radians(1*30))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(2*30))],[0,np.sin(np.radians(2*30))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(3*30))],[0,np.sin(np.radians(3*30))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(4*30))],[0,np.sin(np.radians(4*30))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(5*30))],[0,np.sin(np.radians(5*30))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(6*30))],[0,np.sin(np.radians(6*30))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(7*30))],[0,np.sin(np.radians(7*30))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(8*30))],[0,np.sin(np.radians(8*30))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(9*30))],[0,np.sin(np.radians(9*30))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(10*30))],[0,np.sin(np.radians(10*30))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(11*30))],[0,np.sin(np.radians(11*30))],color='gray',linestyle='-',linewidth=0.25)
ax.axhline(color='black',linestyle='-',linewidth=1)
ax.axvline(color='black',linestyle='-',linewidth=1)
ax.scatter(hx_array_north,hy_array_north,color='tomato',s=3) # plutino poles in the northern hemisphere
ax.scatter(hx_array_south,hy_array_south,color='cadetblue',s=3) # plutino poles in the southern hemisphere
plt.tight_layout()
plt.axis('equal')
ax.set_xlabel('$h_x$')
ax.set_ylabel('$h_y$')
ax.set_xlim([-1,1])
ax.set_ylim([-1,1])
ax.set_box_aspect(1)
# detail near the origin
ax = fig.add_subplot(122)
th = np.linspace(start=0,stop=2*np.pi,num=100,endpoint=True)
costh = np.cos(th)
sinth=  np.sin(th)
ax.plot(np.cos(np.radians(80))*costh,np.cos(np.radians(80))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(70))*costh,np.cos(np.radians(70))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(60))*costh,np.cos(np.radians(60))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(50))*costh,np.cos(np.radians(50))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(40))*costh,np.cos(np.radians(40))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(30))*costh,np.cos(np.radians(30))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(20))*costh,np.cos(np.radians(20))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(10))*costh,np.cos(np.radians(10))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(90-9))*costh,np.cos(np.radians(90-9))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(90-8))*costh,np.cos(np.radians(90-8))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(90-7))*costh,np.cos(np.radians(90-7))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(90-6))*costh,np.cos(np.radians(90-6))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(90-5))*costh,np.cos(np.radians(90-5))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(90-4))*costh,np.cos(np.radians(90-4))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(90-3))*costh,np.cos(np.radians(90-3))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(90-2))*costh,np.cos(np.radians(90-2))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot(np.cos(np.radians(90-1))*costh,np.cos(np.radians(90-1))*sinth,color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(10))],[0,np.sin(np.radians(10))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(20))],[0,np.sin(np.radians(20))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(30))],[0,np.sin(np.radians(30))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(40))],[0,np.sin(np.radians(40))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(50))],[0,np.sin(np.radians(50))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(60))],[0,np.sin(np.radians(60))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(70))],[0,np.sin(np.radians(70))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(80))],[0,np.sin(np.radians(80))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(90))],[0,np.sin(np.radians(90))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(100))],[0,np.sin(np.radians(100))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(110))],[0,np.sin(np.radians(110))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(120))],[0,np.sin(np.radians(120))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(130))],[0,np.sin(np.radians(130))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(140))],[0,np.sin(np.radians(140))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(150))],[0,np.sin(np.radians(150))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(160))],[0,np.sin(np.radians(160))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(170))],[0,np.sin(np.radians(170))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(180))],[0,np.sin(np.radians(180))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(190))],[0,np.sin(np.radians(190))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(200))],[0,np.sin(np.radians(200))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(210))],[0,np.sin(np.radians(210))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(220))],[0,np.sin(np.radians(220))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(230))],[0,np.sin(np.radians(230))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(240))],[0,np.sin(np.radians(240))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(250))],[0,np.sin(np.radians(250))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(260))],[0,np.sin(np.radians(260))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(270))],[0,np.sin(np.radians(270))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(280))],[0,np.sin(np.radians(280))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(290))],[0,np.sin(np.radians(290))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(300))],[0,np.sin(np.radians(300))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(310))],[0,np.sin(np.radians(310))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(320))],[0,np.sin(np.radians(320))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(330))],[0,np.sin(np.radians(330))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(340))],[0,np.sin(np.radians(340))],color='gray',linestyle='-',linewidth=0.25)
ax.plot([0,np.cos(np.radians(350))],[0,np.sin(np.radians(350))],color='gray',linestyle='-',linewidth=0.25)
ax.axhline(color='black',linestyle='-',linewidth=1)
ax.axvline(color='black',linestyle='-',linewidth=1)
ax.scatter(hx_array_north,hy_array_north,color='tomato',s=5) # plutino poles on the northern hemisphere
ax.scatter(hx_array_south,hy_array_south,color='cadetblue',s=5) # plutino poles on the southern hemisphere
ax.plot(hx_circle_post,hy_circle_post,color='blue',linestyle='-',linewidth=1) # confidence circle
ax.scatter(muhat_hx,muhat_hy,color='blue',s=50,marker='o') # vmf midplane of plutinos
ax.scatter(hx_invariable,hy_invariable,color='red',s=50,marker='v') # invariable pole of the solar system
ax.scatter(hx_neptune,hy_neptune,color='blue',s=50,marker='>') # orbit pole of neptune
plt.tight_layout()
plt.axis('equal')
ax.set_xlabel('$h_x$')
ax.set_ylabel('$h_y$')
ax.set_xlim([-0.01,0.09])
ax.set_ylim([-0.01,0.09])
ax.set_box_aspect(1)
titlestr = 'fig2_' + str(Nplut)
plt.tight_layout()
plt.savefig(titlestr + '.eps',dpi=300,bbox_inches='tight',pad_inches=0)
plt.show()
#%%
'''
make separate plot of plutino inclinations relative to vMF midplane
histogram of inclinations, vmf inclination distribution, on same plot
'''
vmf_relative_inclinations_degrees = []
for iobj in range(Nobj):
    polevec_obj = np.array([hx_array[iobj],hy_array[iobj],hz_array[iobj]])
    dotted = np.dot(polevec_vmf,polevec_obj)
    normalized = dotted/np.linalg.norm(polevec_vmf)/np.linalg.norm(polevec_obj)
    inclination_radians = np.arccos(normalized)
    inclination_degrees = np.degrees(inclination_radians)
    vmf_relative_inclinations_degrees.append(inclination_degrees)
vmf_relative_inclinations_degrees = np.array(vmf_relative_inclinations_degrees)
vmf_relative_inclinations_radians = np.radians(vmf_relative_inclinations_degrees)
Nhere = 1000000
rand_vmf = ut.rand_von_mises_fisher(mu=polevec_vmf,kappa=kappa,N=Nhere)
rand_vmf_relative_inclinations_degrees = []
rand_vmf_relative_inclinations_radians = []
for iobj in range(Nhere):
    polevec_obj = rand_vmf[iobj]
    dotted = np.dot(polevec_vmf,polevec_obj)
    normalized = dotted/np.linalg.norm(polevec_vmf)/np.linalg.norm(polevec_obj)
    inclination_radians = np.arccos(normalized)
    inclination_degrees = np.degrees(inclination_radians)
    rand_vmf_relative_inclinations_degrees.append(inclination_degrees)
    rand_vmf_relative_inclinations_radians.append(inclination_radians)
rand_vmf_relative_inclinations_degrees = np.array(rand_vmf_relative_inclinations_degrees)
rand_vmf_relative_inclinations_radians = np.array(rand_vmf_relative_inclinations_radians)
# Anderson-Darling test mentioned in first paragraph of section 3.4
AD_result = stats.anderson_ksamp([vmf_relative_inclinations_radians,\
                                  rand_vmf_relative_inclinations_radians])
degrees_vec = np.linspace(0,np.max(rand_vmf_relative_inclinations_degrees),num=1000000,endpoint=True)
radians_vec = np.linspace(0,np.max(rand_vmf_relative_inclinations_radians),num=1000000,endpoint=True)
# equation 22 in paper
curve_vec = kappa/(np.exp(kappa)-np.exp(-kappa))*\
    np.exp(kappa*np.cos(radians_vec))*np.sin(radians_vec)
curve_vec = np.radians(curve_vec) # x-axis is stretched by radians->degrees for plotting,
# so y axis needs to be squeezed by that amount to have same area under plot
fig = plt.figure(figsize=(2,2))
ax1  = fig.add_subplot(111)
ax1.set_xlabel('Relative inclination (degrees)')
# ax1.set_xticks(np.radians([0,10,20,30,40]))
# ax1.set_xticklabels(['0','10','20','30','40'])
plt.rcParams['font.size'] = 6
data = vmf_relative_inclinations_degrees
myHist = ax1.hist(data, bins=20,density=True,histtype='bar',ec='black',alpha=0.25,zorder=3)
# plot exact relative inclination pdf of vmf distribution using kappa from sample
g = ax1.plot(degrees_vec[0:-1:100],curve_vec[0:-1:100],color='gold',zorder=3,lw=1)
# plot truncated rayleigh pdf distribution using least-squares kappa
initial_guess = kappa
xdata = myHist[1]
xdata = xdata[0:len(xdata)-1]
ydata = myHist[0]
dx = xdata[1]-xdata[0]
print('bin_width = ',dx,'deg')
xdata = np.radians(xdata)
ydata = np.degrees(ydata)
popt,pcov = optimization.curve_fit(f=ut.truncated_rayleigh_pdf, xdata=xdata, ydata=ydata, \
       p0=initial_guess)
best_fit_kappa = popt[0]
print('best_fit_kappa =',best_fit_kappa)
print('best_fit_sigma =',np.degrees(best_fit_kappa**-0.5),'deg')
x = np.linspace(0,np.max(vmf_relative_inclinations_radians),num=1000)
curve_vec_2 = ut.truncated_rayleigh_pdf(x,best_fit_kappa)
curve_vec_2 = np.radians(curve_vec_2) # x-axis is stretched by radians->degrees for plotting,
# so y axis needs to be squeezed by that amount to have same area under plot
h = ax1.plot(np.degrees(x[0:-1:10]),curve_vec_2[0:-1:10],color='yellowgreen',zorder=15,lw=1)
maxx = np.max([np.max(vmf_relative_inclinations_degrees),np.max(np.degrees(i_array_post))])
maxy = np.max([np.max(curve_vec),np.max(curve_vec_2),np.max(myHist[0])])
maxx = 1.1*maxx
maxy = 1.1*maxy
ax1.set_xlim([0,maxx])
ax1.set_ylim([0,maxy])
titlestr = 'fig3_' + str(Nres)
plt.savefig(titlestr + '.eps',dpi=300,bbox_inches='tight',pad_inches=0)
plt.show()
#%%
'''
check for uniformity of W using Rayleigh z test
http://webspace.ship.edu/pgmarr/Geo441/Lectures/Lec%2016%20-%20Directional%20Statistics.pdf
https://www.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics
https://arxiv.org/pdf/1310.5457.pdf
'''
sinWbar = (np.sum(np.sin(W_array)))/Nobj
cosWbar = (np.sum(np.cos(W_array)))/Nobj
r_Rayleighztest = np.sqrt(sinWbar**2 + cosWbar**2)
R_Rayleighztest = Nobj * r_Rayleighztest
pval = np.exp(np.sqrt(1+4*Nobj+4*(Nobj**2-R_Rayleighztest**2))-(1+2*Nobj))
print('Rayleigh z-test pval =',pval)
# Construct figure and axis to plot on
fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
# Visualise by area of bins
n,bins,patches=ut.circular_hist(ax,W_array,bins=10)
titlestr = 'fig4_' + str(Nplut)
plt.savefig(titlestr + '.eps',dpi=300,bbox_inches='tight',pad_inches=0)
#%%
'''
get a simple bootstrap confidence interval for kappa
'''
Nboot = 1000
bootstrapped_kappas = []
for iboot in range(Nboot):
    hx_array_boot = []
    hy_array_boot = []
    hz_array_boot = []
    for iobj in range(Nobj):
        rand = random.randint(0,Nobj-1)
        hx_array_boot.append(hx_array[rand])
        hy_array_boot.append(hy_array[rand])
        hz_array_boot.append(hz_array[rand])
    hx_array_boot = np.array(hx_array_boot)
    hy_array_boot = np.array(hy_array_boot)
    hz_array_boot = np.array(hz_array_boot)
    # equation 9 in thee paper
    hxbar_boot = np.mean(hx_array_boot)
    hybar_boot = np.mean(hy_array_boot)
    hzbar_boot = np.mean(hz_array_boot)
    # equation 10 in the paper
    R_boot = np.sqrt(hxbar_boot**2+hybar_boot**2+hzbar_boot**2)
    # equation 11 in the paper
    muhat_hx_boot = hxbar_boot/R_boot
    muhat_hy_boot = hybar_boot/R_boot
    muhat_hz_boot = hzbar_boot/R_boot
    # equation 12 in the paper
    kappa_boot = R_boot*(3-R_boot**2)/(1-R_boot**2)
    diff_boot = 1
    tol_boot = 1e-8
    while abs(diff_boot) > tol_boot:
        top_boot = special.iv(3/2, kappa_boot)
        bottom_boot = special.iv(1/2, kappa_boot)
        # equation 14 in the papeer
        A3_boot = top_boot/bottom_boot
        # equation 13 in th paper
        diff_boot =(A3_boot-R_boot)/(1-A3_boot**2-2/kappa_boot*A3_boot)
        kappa_boot = kappa_boot - diff_boot
    bootstrapped_kappas.append(kappa_boot)
bootstrapped_kappas = np.array(bootstrapped_kappas)
bootstrapped_sigmas = np.degrees(1/np.sqrt(bootstrapped_kappas))
print('16th percentile bootstrapped kappa =',np.percentile(bootstrapped_kappas,16))
print('84th percentile bootstrapped kappa =',np.percentile(bootstrapped_kappas,84))
print('16th percentile bootstrapped sigma =',np.percentile(bootstrapped_sigmas,16))
print('84th percentile bootstrapped sigma =',np.percentile(bootstrapped_sigmas,84))
#%%
'''
Rayleigh parameter estimation from Wikipedia, which references
https://nvlpubs.nist.gov/nistpubs/jres/66D/jresv66Dn2p167_A1b.pdf
JOURNAL OF RESEARCH of the National Bureau of Standards-D. Radio Propagation
Vol. 66D, No.2, March- April 1962
Some Problems Connected With Rayleigh Distributions
M. M. Siddiqui
Contribution from Boulder Laboratories, National Bureau of Standards, Boulder, Colo.
(October 19, 1961)
'''
# paragraph between equations 3.10 and 3.11 in the reference
zvec = np.radians(vmf_relative_inclinations_degrees)**2 # with variance gamma/2, so std = sqrt(gamma)/sqrt(2)
# equation 4.7 in the reference
ccc = np.sum(zvec)/Nobj # maximum likelihood estimate of gamma
# equation 4.10 in the reference
gamma_bound_1 = 2*Nobj*ccc/stats.chi2.isf(0.84,df=2*Nobj)
gamma_bound_2 = 2*Nobj*ccc/stats.chi2.isf(0.16,df=2*Nobj)
std = np.sqrt(ccc)/np.sqrt(2)
std_bound_1 = np.sqrt(gamma_bound_1)/np.sqrt(2)
std_bound_2 = np.sqrt(gamma_bound_2)/np.sqrt(2)
print('siddiqui statistics',np.degrees([std_bound_2,std,std_bound_1]))
#%% compare vmf inclination and truncated Rayleigh inclination for appendix B
kappa_vec = np.array([30,15,4,1])
colors = ['green','brown','magenta','blue']
degree_vec = np.linspace(start=0,stop=180,num=1000,endpoint=False)
fig = plt.figure(figsize=(3,2))
ax1  = fig.add_subplot(111)
ax1.set_xlabel('Inclination (degrees)')
ax1.set_ylabel('pdf')
plt.rcParams['font.size'] = 6
axins = inset_axes(ax1,width="50%",height="50%",borderpad=1)
axins.set_xlabel('Inclination (degrees)')
axins.set_ylabel('Rayleigh - VMF')
mindiff = []
maxdiff = []
maxy = []
for ik in range(4):
    kappa_here = kappa_vec[ik]
    color = colors[ik]
    sigma_here = 1/np.sqrt(kappa_here)
    # equation 22 in paper
    curve_vec_22 = kappa_here/(np.exp(kappa_here)-np.exp(-kappa_here))*\
        np.exp(kappa_here*np.cos(np.radians(degree_vec)))*np.sin(np.radians(degree_vec))
    curve_vec_22 = np.radians(curve_vec_22)
    # equation 28 in paper
    C_R = 1/( 1 - np.exp(-np.pi**2/2/sigma_here**2) )
    # equation 27 in paper
    curve_vec_27 = C_R/sigma_here**2 * np.radians(degree_vec) * \
        np.exp(-np.radians(degree_vec)**2/2/sigma_here**2)
    curve_vec_27 = np.radians(curve_vec_27)
    h = ax1.plot(degree_vec[0:-1:10],curve_vec_22[0:-1:10],color=color,zorder=1,lw=0.5,linestyle='solid')
    h2 = ax1.plot(degree_vec[0:-1:10],curve_vec_27[0:-1:10],color=color,zorder=2,lw=0.5,linestyle='dotted')
    h3 = axins.plot(degree_vec[0:-1:10],curve_vec_27[0:-1:10]-curve_vec_22[0:-1:10],color=color,\
                    lw=0.5,linestyle='solid')
    mindiff.append(np.min(curve_vec_27-curve_vec_22))
    maxdiff.append(np.max(curve_vec_27-curve_vec_22))
    maxy.append(1.1*np.max([np.max(curve_vec_22),np.max(curve_vec_27)]))
ax1.set_xlim([0,180])
ax1.set_ylim([0,np.max(maxy)])
axins.set_xlim([0,180])
titlestr = 'figB1'
plt.savefig(titlestr + '.eps',dpi=300,bbox_inches='tight',pad_inches=0)
plt.show()
