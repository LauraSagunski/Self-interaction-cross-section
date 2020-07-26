from __future__ import division

import numpy as np
from numpy import sqrt, heaviside, pi, sin, cos, log, exp
from scipy.special import kn, erfi, lambertw

eta = lambda x: x**2.*(kn(0,x)*kn(2,x) - kn(1,x)**2.)
print(eta(0.01))

gamma = lambda beta: sqrt(pi/beta) * erfi(sqrt(lambertw(2.*beta)/2.))-1/lambertw(2.*beta)

lmin = lambda beta, kappa: max(1.,beta*kappa)
lminp = lambda beta, kappa: max(3./4.,beta*kappa)

lambdaT = (1.+cos(2.)+2*sin(2.))/2.
lambdaV = (9.-cos(4.)-4.*sin(4.))/16.

sigmaT_smallbeta = lambda beta, kappa: 2. * beta**2. * ((lmin(beta,kappa)**2-1)/(2.* beta**2 * kappa**2) +eta(lmin(beta,kappa)/kappa))

sigmaV_smallbeta = lambda beta, kappa: 4. * beta**2. * ((lminp(beta,kappa)**2-9./16.)/(2. * beta**2 * kappa**2)+eta(2.*lminp(beta,kappa)/kappa))

def sigmaTatt(beta, kappa):
  if beta < 0.2: return sigmaT_smallbeta(beta,kappa)
  elif beta > 50: return (1. + log(beta)- 1./(2.*log(beta)))**2. - 1./gamma(beta)**2.
  else: return 2.778 * exp(-0.0115 * beta) * lambertw(2*beta)**2.

def sigmaTrep(beta, kappa):
  if beta < 0.2: return sigmaT_smallbeta(beta,kappa)
  elif beta > 10: return lambdaT * lambertw(2.*beta)**2.
  else: return 1.669 * lambertw(2*beta)**1.582

def sigmaVatt(beta, kappa):
  if beta < 0.2: return sigmaV_smallbeta(beta,kappa)
  elif beta > 20: return (1 + log(beta)- 1/(2*log(beta)))**2/2.
  else: return 1.7444 * lambertw(2*beta)**1.447

def sigmaVrep(beta, kappa):
  if beta < 0.2: return sigmaV_smallbeta(beta,kappa)
  elif beta > 5: return lambdaV * lambertw(2*beta)**2 + (lambertw(4*pi*beta**2)**2/4.-lambertw(2*beta)**2)/2.
  else: return 1.520 * lambertw(2*beta)**1.334

def sigma(kappa, beta, mode = 'T', sign = 'attractive'):
  if kappa < 0.5: return 0.
  if mode == 'T':
    if sign == 'attractive': return sigmaTatt(beta, kappa)
    elif sign == 'repulsive': return sigmaTrep(beta, kappa)
    else: print('Sign not recognized') 
  elif mode == 'V':
    if sign == 'attractive': return sigmaVatt(beta, kappa)
    elif sign == 'repulsive': return sigmaVrep(beta, kappa)
    else: print('Sign not recognized') 
  else:
    print('Mode not recognized')
    exit()

