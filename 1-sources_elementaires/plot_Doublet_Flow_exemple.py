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
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

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