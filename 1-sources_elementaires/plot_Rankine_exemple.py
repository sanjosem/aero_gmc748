import numpy as np # scientific computing module
import matplotlib.pyplot as plt # plotting functions module

# Cartesion grid parameter 
nx = 500
ny = 500
x = np.linspace(-2, 2, nx)
y = np.linspace(-2, 2, ny)
X, Y = np.meshgrid(x, y)

THETA = np.arctan2(Y,X)
R = np.sqrt((X)**2+(Y)**2)

# Source  point
Lambda1=2.0
# center location
x1=-0.75
y1=0.
# 
THETA1 = np.arctan2(Y-y1,X-x1)
R1 = np.sqrt((X-x1)**2+(Y-y1)**2)

PSI1 = Lambda1/(2.*np.pi) * THETA1
PHI1 = Lambda1/(2.*np.pi) * np.log(R1)

# Sink point
Lambda2=-Lambda1
# center location
x2=-x1
y2=0.
# 
THETA2 = np.arctan2(Y-y2,X-x2)
R2 = np.sqrt((X-x2)**2+(Y-y2)**2)

PSI2 = Lambda2/(2.*np.pi) * THETA2
PHI2 = Lambda2/(2.*np.pi) * np.log(R2)

# Uniform velocity
Vinf = 1.
PSI3 = Vinf * R * np.sin(THETA)
PHI3 = Vinf * R * np.cos(THETA)

# Nb contours to plot
N=21

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
# plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title



fig,ax = plt.subplots(1,1, sharex=True, sharey=False,figsize=(9,8))
plt.subplots_adjust(hspace=0.1,wspace = 0.0,left  = 0.17,right  = 0.97,bottom  = 0.1,top  = 0.93)

CSS = plt.contour(X,Y,PSI1+PSI2+PSI3,N)
CSP = plt.contour(X,Y,PHI1+PHI2+PHI3,N,linestyles='dashed')
plt.contour(X,Y,PSI1+PSI2+PSI3,[0.0,],colors='red',linewidths=3.)
bbox_props = dict(boxstyle="round", fc="w", ec="0.5", alpha=0.9)
ax.annotate(r'puits', xy=(0.75, 0), xytext=(0.5, -0.75),
            arrowprops=dict(facecolor='black',edgecolor='black', shrink=0.05),
            fontsize=SMALL_SIZE,bbox=bbox_props,zorder=10)
ax.annotate(r'source', xy=(-0.75, 0), xytext=(-0.6, -0.75),
            arrowprops=dict(facecolor='black',edgecolor='black', shrink=0.05),
            fontsize=SMALL_SIZE,bbox=bbox_props,zorder=10)

a=np.sqrt(x2**2+Lambda1*x2/(np.pi*Vinf))
ax.annotate(r'$-a$', xy=(-a, 0), xytext=(-1.4, -1.3),
            arrowprops=dict(facecolor='red',edgecolor='red', shrink=0.05),
            fontsize=BIGGER_SIZE,color='red',bbox=bbox_props,zorder=10)
ax.annotate(r'$a$', xy=(a, 0), xytext=(1.3, -1.3),
            arrowprops=dict(facecolor='red',edgecolor='red', shrink=0.05),
            fontsize=BIGGER_SIZE,color='red',bbox=bbox_props,zorder=10)
plt.xlabel(r'x')
plt.ylabel(r'y')
plt.title(r'Superposition - Ovale Rankine')
plt.axis('equal')
# ax.set_adjustable('box')
# plt.xlim(-1,1)
# plt.ylim(-1,1)
plt.clabel(CSS,inline=1,fontsize=0.7*SMALL_SIZE)
plt.clabel(CSP,inline=1,fontsize=0.7*SMALL_SIZE)
plt.savefig('plot_Rankine.pdf')
plt.show()

# Usage of streamplot
u = Vinf + Lambda1/(2.*np.pi) * ( (X-x1)/((X-x1)**2+Y**2) - (X-x2)/((X-x2)**2+Y**2) )
v = Lambda1/(2.*np.pi) * ( Y/((X-x1)**2+Y**2) - Y/((X-x2)**2+Y**2) )

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
plt.contour(X,Y,PSI1+PSI2+PSI3,[0.0,],colors='red',linewidths=3.)
plt.xlabel(r'x')
plt.ylabel(r'y')
plt.title(r'Superposition - Ovale Rankine')
plt.axis('equal')
# ax.set_adjustable('box')
# plt.xlim(-1,1)
# plt.ylim(-1,1)
plt.savefig('plot_streamplot_Rankine.pdf')
plt.show()
