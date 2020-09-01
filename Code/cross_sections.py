from __future__ import division

import numpy as np
from numpy import sqrt, heaviside, pi, sin, cos, log, exp, euler_gamma
from scipy.special import kn, erfi, lambertw

lmin = lambda beta, kappa: max(1.,beta*kappa)
lminp = lambda beta, kappa: max(3./2.,2.*beta*kappa)
eta = lambda x: 2.*log(2.*x)-1-2.*euler_gamma+(1-euler_gamma+log(2.*x))/x**2.
turn = lambda beta, betalow, a: exp(-(max(beta, betalow) - betalow)*a)

lambdaT = (1.+cos(2.)+2*sin(2.))/2.
lambdaV = (9.-cos(4.)-4.*sin(4.))/16.

sigmaT_smallbeta = lambda beta, kappa: 2. * beta**2. * ((lmin(beta,kappa)**2-1)/(2.* beta**2 * kappa**2) + eta(kappa/lmin(beta,kappa)))

sigmaV_smallbeta = lambda beta, kappa,: 4. * beta**2. * ((lminp(beta,kappa)**2-9./4.)/(8. * beta**2 * kappa**2)+eta(kappa/lminp(beta,kappa)))

def sigmaTatt(beta, kappa):
  if beta < 1: return sigmaT_smallbeta(beta,kappa)*turn(beta,0.2,-0.63854)
  elif beta > 50: return 2. * log(beta) * (log(log(beta)) + 1)
  else: return 4.70847*log(beta + 0.822474)

def sigmaTrep(beta, kappa):
  if beta <1: return sigmaT_smallbeta(beta,kappa)*turn(beta,0.2,0.532971)
  elif beta > 50: return lambdaT * (log(2.*beta)-log(log(2.*beta)))**2.
  else: return 2.89927*log(beta + 0.46495)

def sigmaVatt(beta, kappa):
  if beta < 0.5: return sigmaV_smallbeta(beta,kappa)*turn(beta,0.1,-0.671439)
  elif beta > 25: return (1 + log(beta)- 1/(2.*log(beta)))**2/2.
  else: return 2.53258*log(beta + 1.04944)

def sigmaVrep(beta, kappa):
  if beta < 0.5: return sigmaV_smallbeta(beta,kappa)*turn(beta,0.1,0.370562)
  elif beta > 25: return  log(2. * beta) * (lambdaV * log(2. * beta) - (2.*lambdaV - 1) * log(log(2.*beta)))
  else: return 2.77162*log(beta + 0.801796)

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

