from __future__ import division

import numpy as np
import scipy.integrate as scint
from scipy.interpolate import interp2d
from scipy.interpolate import RectBivariateSpline
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import brentq
import matplotlib.pyplot as plt
from cross_sections import sigma

mode = 'V'
sign = 'attractive'

inputname1 = 'sigma'+mode+'list_'+sign+'_v2.npy'
inputname2 = 'sigma'+mode+'list_'+sign+'_smallkappa.npy'
outputname_plot1 = 'sigma'+mode+'_'+sign+'.pdf'
outputname_plotdiff = 'sigma'+mode+'_'+sign+'_difference.pdf'
outputname_plotrel = 'sigma'+mode+'_'+sign+'_reldifference.pdf'
outputname_plot2 = 'sigma'+mode+'_'+sign+'_Cluster.pdf'
outputname_plot3 = 'sigma'+mode+'_'+sign+'_Dwarf.pdf'
outputname_data = 'sigma'+mode+'list_'+sign+'.txt'

T1 = np.load(inputname1)
T2 = np.load(inputname2)
Tjoined = np.concatenate((T2,T1))

kappalist = np.log10(np.concatenate((np.array([0.5,2,4]),np.logspace(np.log10(5), 2, 10))))
betalist = np.log10(Tjoined[0][:,0])
sigmalist = np.log10(np.array([Tjoined[i][:,1] for i in range(len(kappalist))]))


f = RectBivariateSpline(kappalist, betalist, sigmalist, kx=1,ky=3)

def sigmainter(kappa, beta):
  if kappa < 0.5:
    return 0
  else:
    return 10**f(np.amin([np.log10(kappa),kappalist[-1]]),np.log10(beta))[0,0]

print(sigma(2,0.05,mode, sign),sigmainter(2,0.05))
print(sigma(2,0.25,mode, sign),sigmainter(2,0.25))
print(sigma(2,0.75,mode, sign),sigmainter(2,0.75))
print(sigma(2,2,mode, sign),sigmainter(2,2))
print(sigma(2,10,mode, sign),sigmainter(2,10))
print(sigma(2,100,mode, sign),sigmainter(2,100))

beta0grid = np.logspace(-3,4, 36, endpoint=True)
kappa0grid = np.logspace(0.4,2, 17, endpoint=True)

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
averagedsigmagrid_analytic = []

for kappa0 in kappa0grid:
  for beta0 in beta0grid:

    xfunction1 = lambda x: sigmainter(kappa0 * x, beta0 / x**2.) 
    xfunction2 = lambda x: sigma(kappa0 * x, beta0 / x**2., mode, sign) 

    averagedsigmagrid.append([beta0, kappa0, velocity_averaging(xfunction1)])
    averagedsigmagrid_wrong.append([beta0, kappa0, wrong_velocity_averaging(xfunction1)])
    averagedsigmagrid_analytic.append([beta0, kappa0, velocity_averaging(xfunction2)])

averagedsigmaarray = np.array(averagedsigmagrid)[:,2].reshape((len(kappa0grid),len(beta0grid)))
averagedsigmaarray_wrong = np.array(averagedsigmagrid_wrong)[:,2].reshape((len(kappa0grid),len(beta0grid)))
averagedsigmaarray_analytic = np.array(averagedsigmagrid_analytic)[:,2].reshape((len(kappa0grid),len(beta0grid)))

print(averagedsigmaarray_analytic)

averagedsigmainter = RectBivariateSpline(np.log10(kappa0grid), np.log10(beta0grid), np.log10(averagedsigmaarray))
averagedsigmainter_wrong = RectBivariateSpline(np.log10(kappa0grid), np.log10(beta0grid), np.log10(averagedsigmaarray_wrong))
averagedsigmainter_analytic = RectBivariateSpline(np.log10(kappa0grid), np.log10(beta0grid), np.log10(averagedsigmaarray_analytic))

averagedsigma = lambda x, y: 10**averagedsigmainter(np.log10(x),np.log10(y))[0,0]
averagedsigma_wrong = lambda x, y: 10**averagedsigmainter_wrong(np.log10(x),np.log10(y))[0,0]
averagedsigma_analytic = lambda x, y: 10**averagedsigmainter_analytic(np.log10(x),np.log10(y))[0,0]

fig, ax1 = plt.subplots()

ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlim((1e-3,1e4))
ax1.set_ylim((1e-6,1e2))

ax1.plot(beta0grid, np.array([averagedsigma(10.,beta0) for beta0 in beta0grid]), label='$\overline{\sigma_\mathrm{'+mode+'}}$ ($\kappa_0 = 10$)', color='xkcd:blue',zorder=2)
ax1.plot(beta0grid, np.array([averagedsigma_wrong(10.,beta0) for beta0 in beta0grid]), label=r'$\langle \sigma_\mathrm{'+mode+'} v_\mathrm{rel} \\rangle / \langle v_\mathrm{rel} \\rangle$', color='xkcd:blue',linestyle='--',zorder=2)
ax1.plot(beta0grid, [sigmainter(4./np.sqrt(np.pi)*10.,beta0*np.pi/16.) for beta0 in beta0grid], label=r'$\sigma_\mathrm{'+mode+'}(\langle v_\mathrm{rel} \\rangle)$', color='xkcd:blue',linestyle=':',zorder=2)
ax1.plot(beta0grid, np.array([averagedsigma(2.,beta0) for beta0 in beta0grid]), label='$\overline{\sigma_\mathrm{'+mode+'}}$ ($\kappa_0 = 2$)', color='xkcd:green',zorder=1)

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

ax2.plot(beta0grid, np.array([averagedsigma(10.,beta0) for beta0 in beta0grid]), label='$\overline{\sigma_\mathrm{'+mode+'}}$ ($\kappa_0 = 10$)', color='xkcd:blue',zorder=2)
ax2.plot(beta0grid, np.array([averagedsigma_wrong(10.,beta0) for beta0 in beta0grid]), label=r'$\langle \sigma_\mathrm{'+mode+'} v_\mathrm{rel} \\rangle / \langle v_\mathrm{rel} \\rangle$ ($\kappa_0 = 10$)', color='xkcd:blue',linestyle='--',zorder=2)
ax2.plot(beta0grid, [sigmainter(4./np.sqrt(np.pi)*10.,beta0*np.pi/16.) for beta0 in beta0grid], label=r'$\sigma_\mathrm{'+mode+'}(\langle v_\mathrm{rel} \\rangle)$ ($\kappa_0 = 10$)', color='xkcd:blue',linestyle=':',zorder=2)
ax2.plot(beta0grid, np.array([averagedsigma(2.,beta0) for beta0 in beta0grid]), label='$\overline{\sigma_\mathrm{'+mode+'}}$ ($\kappa_0 = 2$)', color='xkcd:green',zorder=1)

ax2.tick_params(axis='x', labelsize=10)
ax2.tick_params(axis='y', labelsize=10)

ax1.indicate_inset_zoom(ax2)

plt.tight_layout()
plt.savefig(outputname_plot1)
plt.show()

np.savetxt(outputname_data, averagedsigmagrid)

fig, ax1 = plt.subplots()

ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlim((1e-3,1e4))
ax1.set_ylim((1e-6,1e2))

ax1.plot(beta0grid, np.array([averagedsigma(10.,beta0) for beta0 in beta0grid]), label='Numeric ($\kappa_0 = 10$)', color='xkcd:blue',zorder=2)
ax1.plot(beta0grid, np.array([averagedsigma(2.,beta0) for beta0 in beta0grid]), label='Numeric ($\kappa_0 = 2$)', color='xkcd:green',zorder=1)
ax1.plot(beta0grid, np.array([averagedsigma_analytic(10.,beta0) for beta0 in beta0grid]), label='Analytic ($\kappa_0 = 10$)',linestyle='--', color='xkcd:blue',zorder=2)
ax1.plot(beta0grid, np.array([averagedsigma_analytic(2.,beta0) for beta0 in beta0grid]), label='Analytic ($\kappa_0 = 2$)',linestyle='--', color='xkcd:green',zorder=1)

ax1.legend(fontsize=10,loc=2, frameon=False)
ax1.set_xlabel(r'$\beta_0$',fontsize=16)
ax1.set_ylabel(r'$\sigma \, m_\phi^2 / \pi$',fontsize=16)
ax1.tick_params(axis='x', labelsize=14)
ax1.tick_params(axis='y', labelsize=14)

ax2 = ax1.inset_axes([0.5, 0.08, 0.47, 0.47])

ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlim((5e-1,1e2))
ax2.set_ylim((1e-1,2e1))

ax2.plot(beta0grid, np.array([averagedsigma(10.,beta0) for beta0 in beta0grid]), label='Numeric ($\kappa_0 = 10$)', color='xkcd:blue',zorder=2)
ax2.plot(beta0grid, np.array([averagedsigma(2.,beta0) for beta0 in beta0grid]), label='Numeric ($\kappa_0 = 5$)', color='xkcd:green',zorder=1)
ax2.plot(beta0grid, np.array([averagedsigma_analytic(10.,beta0) for beta0 in beta0grid]), label='Analytic ($\kappa_0 = 10$)',linestyle='--', color='xkcd:blue',zorder=2)
ax2.plot(beta0grid, np.array([averagedsigma_analytic(2.,beta0) for beta0 in beta0grid]), label='Analytic ($\kappa_0 = 2$)',linestyle='--', color='xkcd:green',zorder=1)

ax2.tick_params(axis='x', labelsize=10)
ax2.tick_params(axis='y', labelsize=10)

ax1.indicate_inset_zoom(ax2)

plt.tight_layout()
plt.savefig(outputname_plotdiff)
plt.show()

np.savetxt(outputname_data, averagedsigmagrid)


plt.xscale('log')

plt.plot(beta0grid, np.array([(averagedsigma_analytic(5.,beta0) - averagedsigma(5.,beta0))/(averagedsigma(5.,beta0)) for beta0 in beta0grid]), label=r'$\kappa_0 = 5$', color='xkcd:green')
plt.plot(beta0grid, np.array([(averagedsigma_analytic(20.,beta0) - averagedsigma(20.,beta0))/(averagedsigma(20.,beta0)) for beta0 in beta0grid]), label=r'$\kappa_0 = 20$', color='xkcd:blue')
plt.plot(beta0grid, np.array([(averagedsigma_analytic(50.,beta0) - averagedsigma(50.,beta0))/(averagedsigma(50.,beta0)) for beta0 in beta0grid]), label=r'$\kappa_0 = 50$', color='xkcd:purple')

plt.legend(fontsize=14)
plt.xlabel(r'$\beta_0$',fontsize=16)
plt.ylabel('Relative difference',fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

plt.tight_layout()
plt.savefig(outputname_plotrel)
plt.show()

InvGev3tocm2g = 2184e-7
kmsToSpeedOfLight = 3336e-9

beta0 = lambda mp, mx, a, v0: 2. * a * mp / (mx *  (v0 * kmsToSpeedOfLight)**2.)
kappa0 = lambda mp, mx, a, v0: mx * v0 * kmsToSpeedOfLight / (2. * mp)

def sigmaovermx(mp, mx, a, v0):
  return np.pi / mp**2. / mx * averagedsigma(np.clip(kappa0(mp,mx,a,v0),1.,100.), np.clip(beta0(mp,mx,a,v0),1.e-3,1.e4)) * InvGev3tocm2g

v0Cluster = 1900. / (4. / np.sqrt(np.pi))
somCluster = 0.2
v0Dwarf =  50. / (4. / np.sqrt(np.pi))

def FixmpCluster(mx, a):
  sofmp = lambda mp: sigmaovermx(mp, mx, a, v0Cluster) - somCluster
  
  mpmax = mx * v0Cluster * kmsToSpeedOfLight / 2.
  mpmin =   1e-3 * mx *  (v0Cluster * kmsToSpeedOfLight)**2. / (2. * a)

  if mpmin > mpmax:
    return -1
  elif sofmp(mpmin) < 0:
    return 1e-10
  elif sofmp(mpmax) > 0:
    return 1
  else:

    return brentq(sofmp,mpmin,mpmax) 

def mxmin(a):
  mpmax = lambda mx: mx * v0Cluster * kmsToSpeedOfLight / 2.
  sofmp = lambda mx: sigmaovermx(mpmax(mx), mx, a, v0Cluster) - somCluster
  return brentq(sofmp,10.,1000.) 

def mxmax(a):
  mpmin = lambda mx: 1e-3 * mx *  (v0Cluster * kmsToSpeedOfLight)**2. / (2. * a)
  sofmp = lambda mx: sigmaovermx(mpmin(mx), mx, a, v0Cluster) - somCluster
  return brentq(sofmp,10.,1000.)

def mxminDwarf(a):
  rootfunc = lambda mx: (2. * FixmpCluster(mx, a)) / (v0Dwarf * kmsToSpeedOfLight) - mx
  return brentq(rootfunc, mxmin(a), mxmax(a))

print(FixmpCluster(130., 0.3))

mxgrid = np.logspace(1,3, 101, endpoint=True)
alphagrid = np.logspace(-2,0, 101, endpoint=True)

mxminClustergrid = [mxmin(alpha) for alpha in alphagrid]
mxminDwarfgrid = [ mxminDwarf(alpha) for alpha in alphagrid]
mxmaxgrid = [mxmax(alpha) for alpha in alphagrid]

mpgrid = []
sDwarfgrid = []
for alpha in alphagrid:
  mxminD = mxminDwarf(alpha)
  mxminC = mxmin(alpha)
  mxmaxC = mxmax(alpha)

  for mx in mxgrid:

    if mx < mxmaxC:
      if mx > mxminC:  
        mpgrid.append([alpha, mx, FixmpCluster(mx,alpha)])
      if mx > mxminD:
        sDwarfgrid.append([alpha, mx, sigmaovermx(FixmpCluster(mx,alpha),mx,alpha,v0Dwarf)])

fig, ax = plt.subplots()

plt.xscale('log')
plt.yscale('log')
plt.xlim((10,1000))
plt.ylim((0.01,1))

#CS = ax.contour(mxgrid, alphagrid, mparray, [2e-3,5e-3,1e-2,1.5e-2,2.5e-2])

CS = ax.tricontour(np.array(mpgrid)[:,1], np.array(mpgrid)[:,0], np.array(mpgrid)[:,2], [2e-3,5e-3,1e-2,1.5e-2,2.5e-2], zorder=1,colors=('#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'))

ax.clabel(CS, inline=1, fontsize=14)

plt.plot(mxminClustergrid, alphagrid,  color='grey',zorder=3)
plt.plot(mxmaxgrid, alphagrid,  color='grey',zorder=3)

ax.fill_between(np.concatenate(([mxgrid[0]],mxminClustergrid)), np.concatenate(([alphagrid[0]],alphagrid)),1,  color='xkcd:silver',zorder=2)
ax.fill_between(np.concatenate((mxmaxgrid,[mxgrid[-1]])), np.concatenate((alphagrid,[alphagrid[-1]])),0,  color='xkcd:silver',zorder=2)

ax.text(11, 0.6, r'$\kappa_\mathrm{Cluster} < 1$', fontsize=16,  color='grey')
ax.text(200, 0.04, r'$\beta_\mathrm{Cluster} < 10^{-3}$', fontsize=16,  color='grey')

plt.xlabel(r'$m_\chi$', fontsize=16)
plt.ylabel(r'$\alpha_\chi$', fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

plt.tight_layout()

plt.savefig(outputname_plot2)
plt.show()

fig, ax = plt.subplots()

plt.xscale('log')
plt.yscale('log')
plt.xlim((10,1000))
plt.ylim((0.01,1))

CS = ax.tricontour(np.array(sDwarfgrid)[:,1], np.array(sDwarfgrid)[:,0], np.array(sDwarfgrid)[:,2], [10,20,50,200],zorder=1,colors=('#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'))

ax.clabel(CS, inline=1, fmt='%1.0f', fontsize=14)

plt.plot(mxminDwarfgrid, alphagrid,  color='grey',zorder=3)
plt.plot(mxmaxgrid, alphagrid,  color='grey',zorder=3)

ax.fill_between(np.concatenate(([mxgrid[0]],mxminDwarfgrid)), np.concatenate(([alphagrid[0]],alphagrid)),1,  color='xkcd:silver',zorder=2)
ax.fill_between(np.concatenate((mxmaxgrid,[mxgrid[-1]])), np.concatenate((alphagrid,[alphagrid[-1]])),0,  color='xkcd:silver',zorder=2)

ax.text(15, 0.5, r'$\kappa_\mathrm{Dwarf} < 1$', fontsize=16,  color='grey')
ax.text(200, 0.04, r'$\beta_\mathrm{Cluster} < 10^{-3}$', fontsize=16,  color='grey')

plt.xlabel(r'$m_\chi$', fontsize=16)
plt.ylabel(r'$\alpha_\chi$', fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

plt.tight_layout()

plt.savefig(outputname_plot3)
plt.show()

exit()
