from __future__ import division

import numpy as np
from numpy import sqrt, heaviside, pi, sin, cos, log, exp, euler_gamma
from scipy.special import kn, erfi, lambertw
import matplotlib.pyplot as plt

approximate_eta = False

lmin = lambda beta, kappa: max(1./2.,beta*kappa)
lminp = lambda beta, kappa: max(1.,2.*beta*kappa)
turn = lambda beta, betalow, a: exp(-(max(beta, betalow) - betalow)*a)

if approximate_eta:
    eta = lambda x: -2.*log(x/2.)-1-2.*euler_gamma+(1-euler_gamma-log(x/2.))*x**2.
else:
    eta = lambda x: x**2 * (- kn(1,x)**2 + kn(2,x)*kn(0,x))

zeta = lambda kappa, beta, lmin: (max(lmin, beta*kappa)**2 - lmin**2)/(2*kappa**2*beta**2) + eta(max(lmin, beta*kappa)/kappa)

lambdaT = (1.+cos(2.)+2*sin(2.))/2.
lambdaV = (9.-cos(4.)-4.*sin(4.))/16.

sigmaT_smallbeta = lambda beta, kappa: 2. * beta**2. * zeta(kappa, beta, 0.5)

sigmaV_smallbeta = lambda beta, kappa, lmin: 4. * beta**2. * zeta(kappa, 2.*beta, lmin)

def sigmaTatt(beta, kappa):
  if beta < 1: return sigmaT_smallbeta(beta,kappa)*turn(beta,0.2,-0.63854)
  elif beta > 50: return 2. * log(beta) * (log(log(beta)) + 1)
  else: return 4.70847*log(beta + 0.822474)

def sigmaTrep(beta, kappa):
  if beta <1: return sigmaT_smallbeta(beta,kappa)*turn(beta,0.2,0.532971)
  elif beta > 50: return lambdaT * (log(2.*beta)-log(log(2.*beta)))**2.
  else: return 2.89927*log(beta + 0.46495)

def sigmaVatt(beta, kappa, lmin):
  if beta < 0.5: return sigmaV_smallbeta(beta,kappa,lmin)*turn(beta,0.1,-0.671439)
  elif beta > 25: return (1 + log(beta)- 1/(2.*log(beta)))**2/2.
  else: return 2.53258*log(beta + 1.04944)

def sigmaVrep(beta, kappa, lmin):
  if beta < 0.5: return sigmaV_smallbeta(beta,kappa,lmin)*turn(beta,0.1,0.370562)
  elif beta > 25: return  log(2. * beta) * (lambdaV * log(2. * beta) - (2.*lambdaV - 1) * log(log(2.*beta)))
  else: return 2.77162*log(beta + 0.801796)

def sigma(kappa, beta, mode = 'T', sign = 'attractive'):
  if kappa < 0.5: return 0.
  if mode == 'T':
    if sign == 'attractive': return sigmaTatt(beta, kappa)
    elif sign == 'repulsive': return sigmaTrep(beta, kappa)
    else: print('Sign not recognized') 
  elif mode == 'V':
    if sign == 'attractive': return sigmaVatt(beta, kappa, 1.)
    elif sign == 'repulsive': return sigmaVrep(beta, kappa, 1.)
    else: print('Sign not recognized') 
  elif mode == 'even':
    if sign == 'attractive': return sigmaVatt(beta, kappa, 0.5)
    elif sign == 'repulsive': return sigmaVrep(beta, kappa, 0.5)
  elif mode == 'odd':
    if sign == 'attractive': return sigmaVatt(beta, kappa, 1.5)
    elif sign == 'repulsive': return sigmaVrep(beta, kappa, 1.5)
  elif mode == 'scalar':
    return sigma(kappa, beta, mode = 'even', sign = sign)
  elif mode == 'fermion':
    return 0.75*sigma(kappa, beta, mode = 'odd', sign = sign) + 0.25*sigma(kappa, beta, mode = 'even', sign = sign)
  elif mode == 'vector':
    return 1/3.*sigma(kappa, beta, mode = 'odd', sign = sign) + 2/3.*sigma(kappa, beta, mode = 'even', sign = sign)
  else:
    print('Mode not recognized')
    exit()

beta0grid = np.logspace(-2,1, 61, endpoint=True)
plt.xscale('log')
plt.yscale('log')
plt.plot(beta0grid, np.array([2.*sigma(2.,beta0,mode='fermion') for beta0 in beta0grid]),linestyle='--')
plt.plot(beta0grid, np.array([2.*sigma(5.,beta0,mode='fermion') for beta0 in beta0grid]),linestyle='--')
plt.plot(beta0grid, np.array([2.*sigma(20.,beta0,mode='fermion') for beta0 in beta0grid]),linestyle='--')

#plt.plot(beta0grid, np.array([sigma(2.,beta0,mode='scalar') for beta0 in beta0grid]))
#plt.plot(beta0grid, np.array([sigma(2.,beta0,mode='fermion') for beta0 in beta0grid]))
#plt.plot(beta0grid, np.array([sigma(2.,beta0,mode='vector') for beta0 in beta0grid]))
plt.show()
