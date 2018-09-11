import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

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


fig,ax = plt.subplots(1,1, sharex=True, sharey=False,figsize=(9,5.5))
plt.subplots_adjust(hspace=0.1,wspace = 0.0,left  = 0.13,right  = 0.97,bottom  = 0.165,top  = 0.97)

theta=np.linspace(0,2*np.pi,180)
Cp=1-4*np.sin(theta)**2

plt.plot(theta,Cp,linewidth=2.)
plt.xlabel(r'$\theta$ [rad]')
plt.ylabel(r'$C_p$')
ax.set_xticks([0., .5*np.pi, np.pi, 1.5*np.pi, 2*np.pi])
ax.set_xticklabels(["$0$", r"$\frac{\pi}{2}$",
                     r"$\pi$", r"$\frac{3\pi}{2}$", r"$2\pi$"])
plt.grid()
plt.savefig('plot_Cp.png',dpi=300)
plt.show()