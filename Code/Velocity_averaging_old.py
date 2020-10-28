from __future__ import division

import numpy as np
import scipy.integrate as scint
from scipy.interpolate import interp2d
from scipy.interpolate import RectBivariateSpline
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import brentq
import matplotlib.pyplot as plt
import matplotlib.ticker
from cross_sections import sigma_combined as sigma

mode = 'T'
sign = 'attractive'

outputname_plot1 = 'sigma'+mode+'_'+sign+'.pdf'
outputname_plot2 = 'sigma'+mode+'_'+sign+'_Cluster.pdf'
outputname_plot3 = 'sigma'+mode+'_'+sign+'_Dwarf.pdf'
outputname_plot4 = 'sigma'+mode+'_'+sign+'_'
outputname_data = 'sigma'+mode+'list_'+sign+'.txt'

beta0grid = np.logspace(-4,4, 81, endpoint=True)
kappa0grid = np.logspace(-3,3, 61, endpoint=True)

plt.xscale('log')
plt.yscale('log')

def velocity_averaging(xfunction):

  if mode == 'T':
    weighting_function = lambda x: np.exp(-x**2./4.) * x**4./(32.*np.sqrt(2/np.pi))
  elif mode == 'V':
    weighting_function = lambda x: np.exp(-x**2./4.) * x**5./(48.)
  else:
    print('Mode not recognized')
    exit()

  return scint.quad(lambda x: weighting_function(x) * xfunction(x), 0.1, 10)[0]

def wrong_velocity_averaging(xfunction):

  weighting_function = lambda x: np.exp(-x**2./4.) * x**3./(4 * np.sqrt(np.pi))

  return scint.quad(lambda x: weighting_function(x) * xfunction(x), 0.1, 10)[0]

averagedsigmagrid = []
averagedsigmagrid_wrong = []

for kappa0 in kappa0grid:
  for beta0 in beta0grid:

    xfunction1 = lambda x: sigma(kappa0 * x, beta0 / x**2., mode, sign) 

    averagedsigmagrid.append([beta0, kappa0, velocity_averaging(xfunction1)])
    averagedsigmagrid_wrong.append([beta0, kappa0, wrong_velocity_averaging(xfunction1)])

averagedsigmaarray = np.array(averagedsigmagrid)[:,2].reshape((len(kappa0grid),len(beta0grid)))
averagedsigmaarray_wrong = np.array(averagedsigmagrid_wrong)[:,2].reshape((len(kappa0grid),len(beta0grid)))

averagedsigmainter = RectBivariateSpline(np.log10(kappa0grid), np.log10(beta0grid), np.log10(averagedsigmaarray))
averagedsigmainter_wrong = RectBivariateSpline(np.log10(kappa0grid), np.log10(beta0grid), np.log10(averagedsigmaarray_wrong))

averagedsigma = lambda x, y: 10**averagedsigmainter(np.log10(x),np.log10(y))[0,0]
averagedsigma_wrong = lambda x, y: 10**averagedsigmainter_wrong(np.log10(x),np.log10(y))[0,0]

fig, ax1 = plt.subplots()

ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlim((1e-3,1e4))
ax1.set_ylim((1e-6,1e2))

kappaplot = 10.

ax1.plot(beta0grid, np.array([averagedsigma(kappaplot,beta0) for beta0 in beta0grid]), label='$\overline{\sigma_\mathrm{'+mode+'}}$ ($\kappa_0 = 10$)', color='xkcd:blue',zorder=2)
ax1.plot(beta0grid, np.array([averagedsigma_wrong(kappaplot,beta0) for beta0 in beta0grid]), label=r'$\langle \sigma_\mathrm{'+mode+'} v_\mathrm{rel} \\rangle / \langle v_\mathrm{rel} \\rangle$', color='xkcd:blue',linestyle='--',zorder=2)
ax1.plot(beta0grid, [sigma(4./np.sqrt(np.pi)*kappaplot,beta0*np.pi/16., mode, sign) for beta0 in beta0grid], label=r'$\sigma_\mathrm{'+mode+'}(\langle v_\mathrm{rel} \\rangle)$', color='xkcd:blue',linestyle=':',zorder=2)

ax1.legend(fontsize=12,loc=2, frameon=False)
ax1.set_xlabel(r'$\beta_0$',fontsize=16)
ax1.set_ylabel(r'$\sigma \, m_\phi^2 / \pi$',fontsize=16)
ax1.tick_params(axis='x', labelsize=14)
ax1.tick_params(axis='y', labelsize=14)

ax2 = ax1.inset_axes([0.5, 0.08, 0.47, 0.47])

ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlim((5e-1,1e2))
ax2.set_ylim((1e-1,2e1))

ax2.plot(beta0grid, np.array([averagedsigma(kappaplot,beta0) for beta0 in beta0grid]), label='$\overline{\sigma_\mathrm{'+mode+'}}$ ($\kappa_0 = 10$)', color='xkcd:blue',zorder=2)
ax2.plot(beta0grid, np.array([averagedsigma_wrong(kappaplot,beta0) for beta0 in beta0grid]), label=r'$\langle \sigma_\mathrm{'+mode+'} v_\mathrm{rel} \\rangle / \langle v_\mathrm{rel} \\rangle$ ($\kappa_0 = 10$)', color='xkcd:blue',linestyle='--',zorder=2)
ax2.plot(beta0grid, [sigma(4./np.sqrt(np.pi)*kappaplot,beta0*np.pi/16., mode, sign) for beta0 in beta0grid], label=r'$\sigma_\mathrm{'+mode+'}(\langle v_\mathrm{rel} \\rangle)$', color='xkcd:blue',linestyle=':',zorder=2)

ax2.tick_params(axis='x', labelsize=10)
ax2.tick_params(axis='y', labelsize=10)

ax1.indicate_inset_zoom(ax2)

plt.tight_layout()
plt.savefig(outputname_plot1)
plt.show()

np.savetxt(outputname_data, averagedsigmagrid)

InvGev3tocm2g = 2184e-7
kmsToSpeedOfLight = 3336e-9

beta0 = lambda mp, mx, a, v0: 2. * a * mp / (mx *  (v0 * kmsToSpeedOfLight)**2.)
kappa0 = lambda mp, mx, a, v0: mx * v0 * kmsToSpeedOfLight / (2. * mp)

def sigmaovermx(mp, mx, a, v0):
  return np.pi / mp**2. / mx * averagedsigma(kappa0(mp,mx,a,v0), beta0(mp,mx,a,v0)) * InvGev3tocm2g

v0Cluster = 1900. / (4. / np.sqrt(np.pi))
somCluster = 0.2
v0Dwarf =  50. / (4. / np.sqrt(np.pi))

mxgrid = np.logspace(1,3, 101, endpoint=True)
mpgrid = np.logspace(-3,0, 101, endpoint=True)

alpha = 0.5
sDwarfgrid = []
sClustergrid = []
kDwarfgrid = []
bDwarfgrid = []

print(sigmaovermx(0.003,190,0.5,v0Dwarf))
print(sigmaovermx(0.003,190,0.5,v0Cluster))

for mp in mpgrid:

  for mx in mxgrid:

    sDwarfgrid.append([mp, mx, sigmaovermx(mp,mx,alpha,v0Dwarf)])
    sClustergrid.append([mp, mx, sigmaovermx(mp,mx,alpha,v0Cluster)])
    kDwarfgrid.append([mp, mx, kappa0(mp,mx,alpha,v0Dwarf)])
    bDwarfgrid.append([mp, mx, beta0(mp,mx,alpha,v0Cluster)])

fig, ax = plt.subplots()

plt.xscale('log')
plt.yscale('log')
plt.xlim((10,1000))
plt.ylim((0.001,0.02))

CS = ax.tricontour(np.array(sDwarfgrid)[:,1], np.array(sDwarfgrid)[:,0], np.array(sDwarfgrid)[:,2], [0.001,0.1,1,10,50],zorder=1,colors=('#1f77b4'))

fmt = matplotlib.ticker.LogFormatterSciNotation()

ax.clabel(CS, inline=1, fmt=fmt, fontsize=14)

CS2 = ax.tricontour(np.array(sClustergrid)[:,1], np.array(sClustergrid)[:,0], np.array(sClustergrid)[:,2], [1e-6,1e-5,1e-4,1e-3,1e-2,0.1,1],zorder=1,colors=('#ff7f0e'),linestyle='--')

ax.clabel(CS2, inline=1, fmt=fmt, fontsize=14)

ax.contourf(mxgrid, mpgrid, [[sigmaovermx(mp,mx,alpha,v0Dwarf) for mx in mxgrid] for mp in mpgrid], [50,10e10],zorder=3,colors=('#1f77b4'),alpha=0.3)

ax.contourf(mxgrid, mpgrid, [[sigmaovermx(mp,mx,alpha,v0Cluster) for mx in mxgrid] for mp in mpgrid], [1,10e10],zorder=3,colors=('#ff7f0e'),alpha=0.3)

ax.tricontour(np.array(kDwarfgrid)[:,1], np.array(kDwarfgrid)[:,0], np.array(kDwarfgrid)[:,2], [1],zorder=1,colors=('gray'))

ax.tricontour(np.array(bDwarfgrid)[:,1], np.array(bDwarfgrid)[:,0], np.array(bDwarfgrid)[:,2], [100],zorder=1,colors=('gray'))


#plt.plot(mxminDwarfgrid, alphagrid,  color='grey',zorder=3)
#plt.plot(mxmaxgrid, alphagrid,  color='grey',zorder=3)

#ax.fill_between(np.concatenate(([mxgrid[0]],mxminDwarfgrid)), np.concatenate(([alphagrid[0]],alphagrid)),1,  color='xkcd:silver',zorder=2)
#ax.fill_between(np.concatenate((mxmaxgrid,[mxgrid[-1]])), np.concatenate((alphagrid,[alphagrid[-1]])),0,  color='xkcd:silver',zorder=2)

if alpha >= 0.3:
  ax.text(15, 0.0015, r'$\sigma_\mathrm{Dwarf} / m_\chi > 50 \, \mathrm{cm^2/g}$', fontsize=14,  color='#1f77b4')
  ax.text(11, 0.0035, r'$\sigma_\mathrm{Cluster} / m_\chi > 1 \, \mathrm{cm^2/g}$', fontsize=14,  color='#ff7f0e',rotation=-35)

ax.text(100, 0.01, r'$\kappa_\mathrm{Dwarf} = 1$', fontsize=14,  color='gray',rotation=32)

plt.xlabel(r'$m_\chi$', fontsize=16)
plt.ylabel(r'$m_\phi$', fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.title(r'$\alpha_\chi = '+str(alpha)+'$', fontsize=16)
plt.tight_layout()

plt.savefig(outputname_plot4+str(alpha)+'.pdf')
plt.show()

exit()
