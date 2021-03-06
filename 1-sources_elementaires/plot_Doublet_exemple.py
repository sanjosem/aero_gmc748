import numpy as np # scientific computing module
import matplotlib.pyplot as plt # plotting functions module

# Cartesion grid parameter 
nx = 500
ny = 500
x = np.linspace(-2, 2, nx)
y = np.linspace(-2, 2, ny)
X, Y = np.meshgrid(x, y)


# Doublet  point
Kappa=0.1
# 
THETA = np.arctan2(Y,X)
R = np.sqrt((X)**2+(Y)**2)

PSI1 = -Kappa/(2.*np.pi) * np.sin(THETA)/R
PHI1 = Kappa/(2.*np.pi) * np.cos(THETA) / R
# They diverge at the center point

# Nb contours to plot
N=20
# We define the contours to be plotted because of the divergence around 0
Mphi1=PHI1[nx//2+1,ny//2+1]
mphi1=PHI1[-1,-1]
scale=3.5
V1=np.linspace(mphi1*scale,Mphi1/scale,N)
V1D=np.concatenate((-V1[::-1],V1))

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

CSS = plt.contour(X,Y,PSI1,V1D)
CSP = plt.contour(X,Y,PHI1,V1D,linestyles='dashed')
plt.xlabel(r'x')
plt.ylabel(r'y')
plt.title(r'Doublet')
plt.axis('equal')
# ax.set_adjustable('box')
# plt.xlim(-1,1)
# plt.ylim(-1,1)
plt.clabel(CSS,inline=1,fontsize=0.7*SMALL_SIZE)
plt.clabel(CSP,inline=1,fontsize=0.7*SMALL_SIZE)
plt.savefig('plot_Doublet.pdf')
plt.show()
