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
from cross_sections import sigma_combined as sigma

fontsize=14
legendfontsize=14
font = {'size' : fontsize}
rc('font',**font)
rc('text', usetex=True)
rc('font', family='serif', serif='Computer Modern Roman')

fig, axs = plt.subplots(2, 2, figsize=(12,12/1.5))

modes = ['T','V']
signs = ['attractive','repulsive']
color_list = ['green', 'cyan', 'blue']
kappas = [2,5,20]

lx = [1,10]
ly = [1e-20,1e-20]
axs[0,1].plot(lx,ly,color=color_list[0])
axs[0,1].plot(lx,ly,color=color_list[1])
axs[0,1].plot(lx,ly,color=color_list[2])
axs[0,1].legend([r'$\kappa_0=2$', r'$\kappa_0=5$', r'$\kappa_0=20$'],frameon=False)

axs[0,0].plot(lx,ly,color=color_list[2])
axs[0,0].plot(lx,ly,linestyle='--',color=color_list[2])
axs[0,0].plot(lx,ly,linestyle=':',color=color_list[2])
axs[0,0].legend([r'$\overline{\sigma_\mathrm{T}}$', r'$\langle \sigma_\mathrm{T} v_\mathrm{rel}\rangle / \langle v_\mathrm{rel} \rangle$', r'$\sigma_\mathrm{T}(\langle v_\mathrm{rel} \rangle)$'],frameon=False)

beta0grid = np.logspace(-4,4, 81, endpoint=True)
kappa0grid = np.logspace(-3,3, 61, endpoint=True)

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

for i in [0,1]:
  mode = modes[i]
  for j in [0,1]:
    sign = signs[j]

    averagedsigmagrid = []
    averagedsigmagrid_wrong = []

    for kappa0 in kappa0grid:
      for beta0 in beta0grid:

        xfunction1 = lambda x: sigma(kappa0 * x, beta0 / x**2., mode, sign) 

        averagedsigmagrid.append([beta0, kappa0, velocity_averaging(xfunction1)])
        averagedsigmagrid_wrong.append([beta0, kappa0, wrong_velocity_averaging(xfunction1)])

    outputname_data = 'sigma'+mode+'list_'+sign+'.txt'
    np.savetxt(outputname_data, averagedsigmagrid)

    averagedsigmaarray = np.array(averagedsigmagrid)[:,2].reshape((len(kappa0grid),len(beta0grid)))
    averagedsigmaarray_wrong = np.array(averagedsigmagrid_wrong)[:,2].reshape((len(kappa0grid),len(beta0grid)))

    averagedsigmainter = RectBivariateSpline(np.log10(kappa0grid), np.log10(beta0grid), np.log10(averagedsigmaarray))
    averagedsigmainter_wrong = RectBivariateSpline(np.log10(kappa0grid), np.log10(beta0grid), np.log10(averagedsigmaarray_wrong))

    averagedsigma = lambda x, y: 10**averagedsigmainter(np.log10(x),np.log10(y))[0,0]
    averagedsigma_wrong = lambda x, y: 10**averagedsigmainter_wrong(np.log10(x),np.log10(y))[0,0]

    kappaplot = 20

    axs[i,j].plot(beta0grid, np.array([np.pi*averagedsigma_wrong(kappaplot,beta0) for beta0 in beta0grid]), color=color_list[2],linestyle='--',zorder=2)
    axs[i,j].plot(beta0grid, [np.pi*sigma(4./np.sqrt(np.pi)*kappaplot,beta0*np.pi/16., mode, sign) for beta0 in beta0grid],  color=color_list[2],linestyle=':',zorder=2)

    for k in [0,1,2]:

      kappaplot = kappas[k]
    
      axs[i,j].plot(beta0grid, np.array([np.pi*averagedsigma(kappaplot,beta0) for beta0 in beta0grid]), color=color_list[k] ,zorder=2)

axs[0, 0].set_title('Attractive')
axs[0, 1].set_title('Repulsive')

for ax in axs:
    for a in ax:
        a.set_xscale('log')
        a.set_yscale('log')
        a.set_xlim([1e-3,1e3])
        a.set_ylim([1e-5,1e3])

axs[1,0].set_xlabel(r'$\beta_0$', fontsize=18)
axs[1,1].set_xlabel(r'$\beta_0$', fontsize=18)

axs[0,0].set_ylabel(r'$\sigma_\mathrm{T} \, m_\phi^2$', fontsize=18)
axs[1,0].set_ylabel(r'$\sigma_\mathrm{V} \, m_\phi^2$', fontsize=18)

plt.savefig('velocity_averaging.pdf', bbox_inches='tight', pad_inches=0)
plt.show()

exit() 
