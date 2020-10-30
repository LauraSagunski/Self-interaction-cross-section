from __future__ import division

import numpy as np
import scipy.integrate as scint
from scipy.interpolate import interp2d
from scipy.interpolate import RectBivariateSpline
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import brentq
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
import matplotlib.ticker
from cross_sections import averagedsigma

fontsize=14
legendfontsize=14
font = {'size' : fontsize}
rc('font',**font)
rc('text', usetex=True)
rc('font', family='serif', serif='Computer Modern Roman')

mode = 'V'
sign = 'repulsive'

outputname = 'sigma'+mode+'_'+sign+'_'
outputname_data = 'sigma'+mode+'list_'+sign+'.txt'

InvGev3tocm2g = 2184e-7
kmsToSpeedOfLight = 3336e-9

beta0 = lambda mp, mx, a, v0: 2. * a * mp / (mx *  (v0 * kmsToSpeedOfLight)**2.)
kappa0 = lambda mp, mx, a, v0: mx * v0 * kmsToSpeedOfLight / (2. * mp)

def sigmaovermx(mp, mx, a, v0):
  return np.pi / mp**2. / mx * averagedsigma(kappa0(mp,mx,a,v0), beta0(mp,mx,a,v0), mode, sign) * InvGev3tocm2g

v0Cluster = 1900. / (4. / np.sqrt(np.pi))
somCluster = 0.2
v0Dwarf =  50. / (4. / np.sqrt(np.pi))

mxgrid = np.logspace(1,3, 101, endpoint=True)
mpgrid = np.logspace(-4,-1, 101, endpoint=True)

alpha = 0.5
sDwarfgrid = []
sClustergrid = []
kDwarfgrid = []

for mp in mpgrid:

  for mx in mxgrid:

    sDwarfgrid.append([mp, mx, sigmaovermx(mp,mx,alpha,v0Dwarf)])
    sClustergrid.append([mp, mx, sigmaovermx(mp,mx,alpha,v0Cluster)])
    kDwarfgrid.append([mp, mx, kappa0(mp,mx,alpha,v0Dwarf)])

fig, ax = plt.subplots(figsize=(6,4.5))

plt.xscale('log')
plt.yscale('log')
plt.xlim((10,1000))
plt.ylim((0.0005,0.02))

CS = ax.tricontour(np.array(sDwarfgrid)[:,1], np.array(sDwarfgrid)[:,0], np.array(sDwarfgrid)[:,2], [0.001,0.1,1,10,50],zorder=1,colors=('#1f77b4'))

fmt = matplotlib.ticker.LogFormatterSciNotation()

ax.clabel(CS, inline=1, fmt=fmt, fontsize=14)

CS2 = ax.tricontour(np.array(sClustergrid)[:,1], np.array(sClustergrid)[:,0], np.array(sClustergrid)[:,2], [1e-6,1e-5,1e-4,1e-3,1e-2,0.1,1],zorder=1,colors=('#ff7f0e'),linestyle='--')

ax.clabel(CS2, inline=1, fmt=fmt, fontsize=14)

ax.contourf(mxgrid, mpgrid, [[sigmaovermx(mp,mx,alpha,v0Dwarf) for mx in mxgrid] for mp in mpgrid], [50,10e10],zorder=3,colors=('#1f77b4'),alpha=0.3)

ax.contourf(mxgrid, mpgrid, [[sigmaovermx(mp,mx,alpha,v0Cluster) for mx in mxgrid] for mp in mpgrid], [1,10e10],zorder=3,colors=('#ff7f0e'),alpha=0.3)

ax.tricontour(np.array(kDwarfgrid)[:,1], np.array(kDwarfgrid)[:,0], np.array(kDwarfgrid)[:,2], [1],zorder=1,colors=('gray'))



#plt.plot(mxminDwarfgrid, alphagrid,  color='grey',zorder=3)
#plt.plot(mxmaxgrid, alphagrid,  color='grey',zorder=3)

#ax.fill_between(np.concatenate(([mxgrid[0]],mxminDwarfgrid)), np.concatenate(([alphagrid[0]],alphagrid)),1,  color='xkcd:silver',zorder=2)
#ax.fill_between(np.concatenate((mxmaxgrid,[mxgrid[-1]])), np.concatenate((alphagrid,[alphagrid[-1]])),0,  color='xkcd:silver',zorder=2)

if alpha >= 0.3:
  ax.text(15, 0.0015, r'$\sigma_\mathrm{Dwarf} / m_\chi > 50 \, \mathrm{cm^2/g}$', fontsize=14,  color='#1f77b4')
  ax.text(15, 0.004, r'$\sigma_\mathrm{Cluster} / m_\chi > 1 \, \mathrm{cm^2/g}$', fontsize=14,  color='#ff7f0e',rotation=-45)

ax.text(80, 0.004, r'$\kappa_\mathrm{Dwarf} = 1$', fontsize=14,  color='gray',rotation=45)

plt.xlabel(r'$m_\chi$', fontsize=16)
plt.ylabel(r'$m_\phi$', fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.title(r'$\alpha_\chi = '+str(alpha)+'$', fontsize=16)

plt.tight_layout()
plt.savefig(outputname+'alpha_'+str(alpha)+'.pdf')
plt.show()

exit()
