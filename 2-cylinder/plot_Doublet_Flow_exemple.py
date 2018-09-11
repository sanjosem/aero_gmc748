import numpy as np # scientific computing module
import matplotlib.patches as pa # plotting patches
import matplotlib.pyplot as plt # plotting functions module

# Cartesion grid parameter 
nx = 1000
ny = 1000
x = np.linspace(-2, 2, nx)
y = np.linspace(-2, 2, ny)
X, Y = np.meshgrid(x, y)


# Angle of attack of the uniform flow
aoa = np.deg2rad(8.0)

# Doublet  point
Kappa=4.0
# 
THETA = np.arctan2(Y,X)-aoa
R = np.sqrt((X)**2+(Y)**2)

PSI1 = -Kappa/(2.*np.pi) * np.sin(THETA)/R
PHI1 = Kappa/(2.*np.pi) * np.cos(THETA) / R
# They diverge at the center point

# Uniform velocity
Vinf = 1.
PSI2 = Vinf * R * np.sin(THETA)
PHI2 = Vinf * R * np.cos(THETA)

# Nb contours to plot
N=27
# We define the contours to be plotted because of the divergence around 0
V1=np.linspace(PHI1[0,0],PHI1[-1,-1],N)
V2=np.linspace(PHI1[0,0]+PHI2[0,0],PHI1[-1,-1]+PHI2[-1,-1],N)

# Circle radius
circ=pa.Circle((0.,0.), radius=np.sqrt(Kappa/(2.*np.pi*Vinf)),ec='black',color='grey',alpha=0.5)

# Parameter for plot : use latex in labels, and set size for title, ticks and legend
# plt.rc('text', usetex=True)
# plt.rc('font', family='serif')

SMALL_SIZE = 24
MEDIUM_SIZE = 26
BIGGER_SIZE = 28

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


fig,ax = plt.subplots(1,1, sharex=True, sharey=False,figsize=(9,8))
plt.subplots_adjust(hspace=0.1,wspace = 0.0,left  = 0.17,right  = 0.97,bottom  = 0.1,top  = 0.93)

CSS = plt.contour(X,Y,PSI1+PSI2,V2)
CSP = plt.contour(X,Y,PHI1+PHI2,V2,linestyles='dashed')
ax.add_patch(circ)
plt.xlabel(r'x')
plt.ylabel(r'y')
plt.title(r'Superposition - Cylindre')
plt.axis('equal')
# ax.set_adjustable('box')
# plt.xlim(-1,1)
# plt.ylim(-1,1)
plt.clabel(CSS,inline=1,fontsize=0.7*SMALL_SIZE)
plt.clabel(CSP,inline=1,fontsize=0.7*SMALL_SIZE)
plt.savefig('plot_Cylindre.pdf')
plt.show()

a=np.sqrt(2./(np.pi * Vinf))
hide=R<=a
PSI=(PSI1+PSI2)
PHI=(PHI1+PHI2)
mask=np.zeros_like(PSI, dtype=bool)
mask[hide]=True
PSI=np.ma.array(PSI,mask=mask)
PHI=np.ma.array(PHI,mask=mask)

fig,ax = plt.subplots(1,1, sharex=True, sharey=False,figsize=(9,8))
plt.subplots_adjust(hspace=0.1,wspace = 0.0,left  = 0.17,right  = 0.97,bottom  = 0.1,top  = 0.93)

circ=pa.Circle((0.,0.), radius=np.sqrt(Kappa/(2.*np.pi*Vinf)),ec='black',color='grey',alpha=0.5)
CSS = plt.contour(X,Y,PSI,V2,corner_mask=True)
CSP = plt.contour(X,Y,PHI,V2,linestyles='dashed',corner_mask=True)
ax.add_patch(circ)
plt.xlabel(r'x')
plt.ylabel(r'y')
plt.title(r'Superposition - Cylindre Externe')
plt.axis('equal')
# ax.set_adjustable('box')
# plt.xlim(-1,1)
# plt.ylim(-1,1)
plt.clabel(CSS,inline=1,fontsize=0.7*SMALL_SIZE)
plt.clabel(CSP,inline=1,fontsize=0.7*SMALL_SIZE)
plt.savefig('plot_Cylindre_externe.pdf')
plt.show()


# Define velocity
ur=Vinf*np.cos(THETA)*(1.+a**2/R**2-2*a**2/R**2)
utheta=-Vinf*np.sin(THETA)*(1.-a**2/R**2+2*a**2/R**2)
u=ur*np.cos(THETA+aoa)-utheta*np.sin(THETA+aoa)
v=ur*np.sin(THETA+aoa)+utheta*np.cos(THETA+aoa)


# Seeding of particle for trajectory calculation
seedline=np.empty((N,2))
seedline[:,0]=-2
seedline[:,1]=np.linspace(-2,2,N)
# More particle close to axis
seedline2=np.empty((N,2))
seedline2[:,0]=-2
seedline2[:,1]=np.linspace(-0.5,0.5,N)

fig,ax = plt.subplots(1,1, sharex=True, sharey=False,figsize=(9,8))
plt.subplots_adjust(hspace=0.1,wspace = 0.0,left  = 0.17,right  = 0.97,bottom  = 0.1,top  = 0.93)
# Function to plot streamlines of a velocity vector field. The streamlines are emitted from the point locations provided in the start_point argument. 
plt.streamplot(x, y, u, v, density=1, linewidth=1, arrowsize=1, arrowstyle='->',start_points=seedline, minlength=0.1)
plt.streamplot(x, y, u, v, color='C0', density=5, linewidth=1, arrowsize=1, arrowstyle='->',start_points=seedline2, minlength=0.1)
circ=pa.Circle((0.,0.), radius=np.sqrt(Kappa/(2.*np.pi*Vinf)),ec='black',color='grey',alpha=0.5)
ax.add_patch(circ)
plt.contour(X,Y,PSI,[0.0,],colors='red',linewidths=3.)
plt.xlabel(r'x')
plt.ylabel(r'y')
plt.title(r'Superposition - Cylindre')
plt.axis('equal')
# ax.set_adjustable('box')
# plt.xlim(-1,1)
# plt.ylim(-1,1)
plt.savefig('plot_streamplot_Cylindre.pdf')
plt.show()
