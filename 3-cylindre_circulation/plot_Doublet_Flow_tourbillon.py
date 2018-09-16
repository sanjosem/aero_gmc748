import numpy as np # scientific computing module
import matplotlib.patches as pa # plotting patches
import matplotlib.pyplot as plt # plotting functions module

# Cartesion grid parameter 
nx = 1000
ny = 1000
x = np.linspace(-6, 6, nx)
y = np.linspace(-6, 6, ny)
X, Y = np.meshgrid(x, y)


# Angle of attack of the uniform flow
aoa = np.deg2rad(5.0)

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

a=np.sqrt(2./(np.pi * Vinf))

# Compute velocity field for stream plot (here without circulation which added in the loop)
Vr=(1.-a**2/R**2)*Vinf*np.cos(THETA)
Vt=-(1.+a**2/R**2)*Vinf*np.sin(THETA) 

# Vortex
for Gamma in np.linspace(0.,15.,10):
	PSI3 = Gamma/(2. *np.pi) * np.log(R/a)
	PHI3 = - Gamma/(2. *np.pi) * THETA

	# Nb contours to plot
	N=31
	PSI=(PSI1+PSI2+PSI3)
	PHI=(PHI1+PHI2+PHI3)
	# We define the contours to be plotted because of the divergence around 0
	V2=np.linspace(np.amin(PHI)/50.,np.amax(PHI)/50.,N)


	# Seeding of particle for trajectory calculation
	seedline=np.empty((N,2))
	seedline[:,0]=-6
	seedline[:,1]=np.linspace(-6,6,N)


	# Circle radius
	circ=pa.Circle((0.,0.), radius=np.sqrt(Kappa/(2.*np.pi*Vinf)),ec='black',color='grey',alpha=0.5)

	# Parameter for plot : use latex in labels, and set size for title, ticks and legend
	# plt.rc('text', usetex=True)
	# plt.rc('font', family='serif')

	SMALL_SIZE = 22
	MEDIUM_SIZE = 24
	BIGGER_SIZE = 26

	plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
	plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
	plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
	plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
	plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
	plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
	# plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title



	hide=R<=a
	mask=np.zeros_like(PSI, dtype=bool)
	mask[hide]=True
	PSI=np.ma.array(PSI,mask=mask)
	PHI=np.ma.array(PHI,mask=mask)

	fig,ax = plt.subplots(1,1, sharex=True, sharey=False,figsize=(9,8))
	plt.subplots_adjust(hspace=0.1,wspace = 0.0,left  = 0.17,right  = 0.97,bottom  = 0.1,top  = 0.93)

	circ=pa.Circle((0.,0.), radius=np.sqrt(Kappa/(2.*np.pi*Vinf)),ec='black',color='grey',alpha=0.5)
	CSS = plt.contour(X,Y,PSI,N,corner_mask=True)
	CSS2 = plt.contour(X,Y,PSI,[-0.05,0.0,0.05],corner_mask=True)
	CSP = plt.contour(X,Y,PHI,N,linestyles='dashed',corner_mask=True)
	ax.add_patch(circ)
	plt.xlabel(r'x')
	plt.ylabel(r'y')
	plt.title(r'Superposition - Cylindre + Tourbillon $\Gamma=${0:.1f}'.format(Gamma))
	plt.axis('equal')
	# ax.set_adjustable('box')
	# plt.xlim(-1,1)
	# plt.ylim(-1,1)
	plt.clabel(CSS,inline=1,fontsize=0.7*SMALL_SIZE)
	plt.clabel(CSP,inline=1,fontsize=0.7*SMALL_SIZE)
	# plt.clabel(CSS2,inline=1,fontsize=0.7*SMALL_SIZE)
	plt.savefig('plot_Cylindre_externe_tourbillon_Gamma{0:s}.png'.format('{0:.1f}'.format(Gamma).replace('.','_')))
	plt.show()


	# Compute velocity field in the unrotated plane (+ aoa)
	Vt_with_circ = Vt - Gamma/(2*np.pi*R)
	U = Vr * np.cos(THETA+aoa) - Vt_with_circ * np.sin(THETA+aoa)
	V = Vr * np.sin(THETA+aoa) + Vt_with_circ * np.cos(THETA+aoa)
	Vmag = np.sqrt(U**2+V**2)
	Vmag[R<a] = 0.0 # mask Vmag for color

	fig,ax = plt.subplots(1,1, sharex=True, sharey=False,figsize=(9,9))
	plt.subplots_adjust(hspace=0.1,wspace = 0.0,left  = 0.17,right  = 0.97,bottom  = 0.1,top  = 0.93)

	circ=pa.Circle((0.,0.), radius=np.sqrt(Kappa/(2.*np.pi*Vinf)),ec='black',color='grey',alpha=0.5)
	plt.streamplot(x, y, U, V, color='black', density=3, linewidth=1.5, arrowsize=1, arrowstyle='->',start_points=seedline, minlength=0.3, cmap='jet')
	im=plt.contourf(X,Y,Vmag,cmap='hot')
	# plt.streamplot(x, y, U, V, color='C0', density=2, linewidth=1, arrowsize=1, arrowstyle='->',start_points=seedline2, minlength=0.1)
	ax.add_patch(circ)
	plt.xlabel(r'x')
	plt.ylabel(r'y')
	plt.title(r'Superposition - Cylindre + Tourbillon $\Gamma=${0:.1f}'.format(Gamma))
	# plt.xlim(-1,1)
	# plt.ylim(-1,1)
	plt.colorbar(im,shrink=0.5,ax=ax)
	plt.axis('equal')
	ax.set_adjustable('box')
	# plt.clabel(CSS2,inline=1,fontsize=0.7*SMALL_SIZE)
	plt.savefig('streamplot_Cylindre_externe_tourbillon_Gamma{0:s}.png'.format('{0:.1f}'.format(Gamma).replace('.','_')))
	plt.show()


