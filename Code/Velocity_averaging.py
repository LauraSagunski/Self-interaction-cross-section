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

modes = ['T','V','even','odd','scalar','fermion','vector']
signs = ['attractive','repulsive']

beta0grid = np.logspace(-5,5, 101, endpoint=True)
kappa0grid = np.logspace(-3,3, 61, endpoint=True)

def velocity_averaging(xfunction):

  if mode == 'T':
    weighting_function = lambda x: np.exp(-x**2./4.) * x**4./(32.*np.sqrt(2/np.pi))
  elif mode in modes:
    weighting_function = lambda x: np.exp(-x**2./4.) * x**5./(48.)
  else:
    print('Mode not recognized')
    exit()

  return scint.quad(lambda x: weighting_function(x) * xfunction(x), 0.1, 10)[0]

def wrong_velocity_averaging(xfunction):

  weighting_function = lambda x: np.exp(-x**2./4.) * x**3./(4 * np.sqrt(np.pi))

  return scint.quad(lambda x: weighting_function(x) * xfunction(x), 0.1, 10)[0]

for i, mode in enumerate(modes):
  for j, sign in enumerate(signs):

    averagedsigmagrid = []
    averagedsigmagrid_wrong = []

    for kappa0 in kappa0grid:
      for beta0 in beta0grid:

        xfunction1 = lambda x: sigma(kappa0 * x, beta0 / x**2., mode, sign) 

        averagedsigmagrid.append([beta0, kappa0, velocity_averaging(xfunction1)])
        averagedsigmagrid_wrong.append([beta0, kappa0, wrong_velocity_averaging(xfunction1)])

    outputname_data = 'sigma'+mode+'list_'+sign+'.txt'
    np.savetxt(outputname_data, averagedsigmagrid)

exit() 
