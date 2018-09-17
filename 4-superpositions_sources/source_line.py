import numpy as np # scientific computing module
import matplotlib.patches as pa # plotting patches
import matplotlib.pyplot as plt # plotting functions module
import elementary_flows as ef # bank of potential flow defined in the lessons

# Cartesion grid parameter 
nx = 1000
ny = 1000
xbounds=(-3,3)
ybounds=(-3,3)
x = np.linspace(xbounds[0],xbounds[1], nx)
y = np.linspace(ybounds[0],ybounds[1], ny)
X, Y = np.meshgrid(x, y)

# Generation of the elementary sources

# Uniform flow
Vinf = 1.0
uniform_flow = ef.Uniform(Vinf)
# Compute velocity 
uniform_flow.velocity(X,Y)
# Compute stream function 
uniform_flow.stream_function(X,Y)


# Array of point sources
N_sources = 5
total_strength = 5.0

# Create an array of five sources side by side of equal intensity
strength_source = (total_strength / float(N_sources)) * np.ones((N_sources,1),dtype=float)

xs = np.zeros((N_sources,1),dtype=float) # x center
ys = np.linspace(-1.0,1.0,N_sources) # y center
sources =  np.empty(N_sources,dtype=object)

for i in range(N_sources):
	sources[i] = ef.Source(strength_source[i],xs[i],ys[i])
	sources[i].velocity(X,Y)
	sources[i].stream_function(X,Y)

# Principle of superposition : usage of copy is required here to avoid issue when restarting script
U = (uniform_flow.u).copy()
V = (uniform_flow.v).copy()
PSI = (uniform_flow.psi).copy()

for si in sources:
	U += si.u
	V += si.v
	PSI += si.psi

# Find stagnation points: point with the minimum velocity magnitude
# That is closest  to zero since we are in discretized quantities
Vmag=np.sqrt(U**2+V**2)
# We know that there are 2 so we take the two first minimal locations
idx_stag=np.argsort(Vmag,axis=None)[0]
# That index is in a flatten matrix X[idx_stag], Y[idx_stag]
# as if all values were flatten in a vector of size NxN
# in fact this is how it is stored in memory
# If you want (istag,jstag), so that X[jstag,istag]=X.item(idx_stag)
jstag, istag = np.unravel_index(idx_stag, Vmag.shape)
# We note first index j because it is the index on which Y is varying
# Second index is the index on which X is varying
print('min velocity point: (x,y)=({0:.2f},{1:.2f})'.format(x[istag],y[jstag]))

# Plot parameters
Ncontours=27
seedline=np.empty((Ncontours,2))
seedline[:,0]=xbounds[0]
seedline[:,1]=np.linspace(ybounds[0],ybounds[1],Ncontours)
seedline2=np.empty((Ncontours,2))
seedline2[:,0]=xbounds[0]
seedline2[:,1]=np.linspace(-1.,1.,Ncontours)

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
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


# fig,ax = plt.subplots(1,1, sharex=True, sharey=False,figsize=(9,8))
# plt.subplots_adjust(hspace=0.1,wspace = 0.0,left  = 0.17,right  = 0.97,bottom  = 0.1,top  = 0.93)

# CSS = plt.contour(X,Y,PSI,Ncontours)
# plt.xlabel(r'x')
# plt.ylabel(r'y')
# plt.title(r'Superposition - Uniform flow + Source + Sink')
# plt.axis('equal')
# # ax.set_adjustable('box')
# # plt.xlim(-1,1)
# # plt.ylim(-1,1)
# plt.clabel(CSS,inline=1,fontsize=0.7*SMALL_SIZE)
# plt.savefig('plot_stream_unif_src_sink.pdf')
# plt.show()


fig,ax = plt.subplots(1,1, sharex=True, sharey=False,figsize=(9,8))
plt.subplots_adjust(hspace=0.1,wspace = 0.0,left  = 0.17,right  = 0.97,bottom  = 0.1,top  = 0.93)
# Function to plot streamlines of a velocity vector field. The streamlines are emitted from the point locations provided in the start_point argument. 
plt.streamplot(x, y, U, V, color=PSI, density=1, linewidth=1, arrowsize=1, arrowstyle='->',start_points=seedline, minlength=0.1)
# We plot two times to have a refined discretisation around axis
plt.streamplot(x, y, U, V, color=PSI, density=5, linewidth=1, arrowsize=1, arrowstyle='->',start_points=seedline2, minlength=0.2)
plt.plot(xs,ys,marker='o',color='red')
plt.plot(x[istag],y[jstag],marker='X',color='black')
plt.xlabel(r'x')
plt.ylabel(r'y')
plt.title(r'Superposition - Uniform flow + Sources')
plt.axis('equal')
# ax.set_adjustable('box')
# plt.xlim(-1,1)
# plt.ylim(-1,1)
plt.savefig('plot_streamplot_unif_line_sources.png')
plt.show()




