import numpy as np # scientific computing module
import matplotlib.pyplot as plt # plotting functions module

# Cartesion grid parameter 
nx = 500
ny = 500
x = np.linspace(-1, 1, nx)
y = np.linspace(-1, 1, ny)
X, Y = np.meshgrid(x, y)

# NB: This field is discontinuous, we will see a cluster of lines due to that discontinuity in -pi/+pi
THETA = np.arctan2(Y,X)
R = np.sqrt(X**2+Y**2)
# To avoid the discontinuity one need to create a theta and r grid, and then compute coordinates
# See commented section below

# # Cylindrical grid parameters
# nr = 200
# nt = 360
# r = np.linspace(0, 2*np.sqrt(2), nr)
# t = np.linspace(0, 2*np.pi, nt, endpoint=False)
# R,THETA = np.meshgrid(r, t)
# X = R * np.cos(THETA)
# Y = R * np.sin(THETA)

# Point Source 
Lambda=1.0
PSI1 = Lambda/(2.*np.pi) * THETA
PHI1 = Lambda/(2.*np.pi) * np.log(R)

# Uniform velocity
aoa = np.deg2rad(5.0)
Vinf = 1.
PSI2 = Vinf * R * np.sin(THETA-aoa)
PHI2 = Vinf * R * np.cos(THETA-aoa)

# Nb contours to plot
N=20

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


# Create a figure (frame is adjust to let the text)
fig,ax = plt.subplots(1,1, sharex=True, sharey=False,figsize=(9,8))
plt.subplots_adjust(hspace=0.1,wspace = 0.0,left  = 0.17,right  = 0.97,bottom  = 0.1,top  = 0.93)

CSS = plt.contour(X,Y,PSI1,N)
CSP = plt.contour(X,Y,PHI1,N,linestyles='dashed')
plt.xlabel(r'x')
plt.ylabel(r'y')
plt.title(r'Source Ponctuelle')
plt.axis('equal')
# ax.set_adjustable('box')
# plt.xlim(-1,1)
# plt.ylim(-1,1)
plt.clabel(CSS,inline=1,fontsize=0.7*SMALL_SIZE)
plt.clabel(CSP,inline=1,fontsize=0.7*SMALL_SIZE)
plt.savefig('plot_source_ponctuelle.pdf')
plt.show()

fig,ax = plt.subplots(1,1, sharex=True, sharey=False,figsize=(9,8))
plt.subplots_adjust(hspace=0.1,wspace = 0.0,left  = 0.17,right  = 0.97,bottom  = 0.1,top  = 0.93)

CSS = plt.contour(X,Y,PSI2,N)
CSP = plt.contour(X,Y,PHI2,N,linestyles='dashed')
plt.xlabel(r'x')
plt.ylabel(r'y')
plt.title(r'\'{E}coulement Uniforme')
plt.axis('equal')
# ax.set_adjustable('box')
# plt.xlim(-1,1)
# plt.ylim(-1,1)
plt.clabel(CSS,inline=1,fontsize=0.7*SMALL_SIZE)
plt.clabel(CSP,inline=1,fontsize=0.7*SMALL_SIZE)
plt.savefig('plot_ecoulement_uniforme.pdf')
plt.show()

# To superimpose fields, we need to aligne the velocity with the reference frame of the source
# Cartesion grid parameter 
nx = 500
ny = 500
x = np.linspace(-1, 1, nx)
y = np.linspace(-1, 1, ny)
X, Y = np.meshgrid(x, y)

# Angle rotation by the angle of attack
THETA = np.arctan2(Y,X)-aoa
R = np.sqrt(X**2+Y**2)

# We recompute the fields for these coordinates (parameters as defined before only R, THETA change)
# Point Source 
PSI1 = Lambda/(2.*np.pi) * THETA
PHI1 = Lambda/(2.*np.pi) * np.log(R)

# Uniform velocity
PSI2 = Vinf * R * np.sin(THETA)
PHI2 = Vinf * R * np.cos(THETA)


fig,ax = plt.subplots(1,1, sharex=True, sharey=False,figsize=(9,8))
plt.subplots_adjust(hspace=0.1,wspace = 0.0,left  = 0.17,right  = 0.97,bottom  = 0.1,top  = 0.93)

CSS = plt.contour(X,Y,PSI1+PSI2,N)
CSP = plt.contour(X,Y,PHI1+PHI2,N,linestyles='dashed')
plt.xlabel(r'x')
plt.ylabel(r'y')
plt.title(r'Superposition')
plt.axis('equal')
# ax.set_adjustable('box')
# plt.xlim(-1,1)
# plt.ylim(-1,1)
plt.clabel(CSS,inline=1,fontsize=0.7*SMALL_SIZE)
plt.clabel(CSP,inline=1,fontsize=0.7*SMALL_SIZE)
plt.savefig('plot_superposition.pdf')
plt.show()