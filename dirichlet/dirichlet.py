# Copyright (C) 2012 Eric J. Suh
#
# This file is subject to the terms and conditions defined in file
# 'LICENSE.txt', which is part of this source code package.

'''Dirichlet.py

Maximum likelihood estimation and likelihood ratio tests of Dirichlet
distribution models of data.

Most of this package is a port of Thomas P. Minka's wonderful Fastfit MATLAB
code. Much thanks to him for that and his clear paper "Estimating a Dirichlet
distribution". See the following URL for more information:

    http://research.microsoft.com/en-us/um/people/minka/'''

import sys
import scipy as sp
import scipy.stats as stats
from scipy.special import (psi, polygamma, gammaln)
from numpy import (array, asanyarray, ones, arange, log, diag, vstack, exp,
        asarray, ndarray, zeros, isscalar, var)
from numpy.linalg import norm
import numpy as np
from . import simplex

try:
    # python 2
    MAXINT = sys.maxint
except AttributeError:
    # python 3
    MAXINT = sys.maxsize

try:
    # python 2
    xrange
except NameError:
    # python 3
    xrange = range

__all__ = [
    'pdf',
    'test_same_distribution',
    'test_uniform',
    'plot',
    'mle',
    'meanprecision',
    'loglikelihood',
]

euler = -1*psi(1) # Euler-Mascheroni constant

def test_same_distribution(D1, D2, method='meanprecision', maxiter=None):
    '''Test for statistical difference between observed proportions.

    Parameters
    ----------
    D1 : array
    D2 : array
        Both ``D1`` and ``D2`` must have the same number of columns, which are
        the different levels or categorical possibilities. Each row of the
        matrices must add up to 1.
    method : string
        One of ``'fixedpoint'`` and ``'meanprecision'``, designates method by
        which to find MLE Dirichlet distribution. Default is
        ``'meanprecision'``, which is faster.
    maxiter : int
        Maximum number of iterations to take calculations. Default is
        ``sys.maxint``.

    Returns
    -------
    D : float
        Test statistic, which is ``-2 * log`` of likelihood ratios.
    p : float
        p-value of test.
    a0 : array
    a1 : array
    a2 : array
        MLE parameters for the Dirichlet distributions fit to 
        ``D1`` and ``D2`` together, ``D1``, and ``D2``, respectively.'''

    N1, K1 = np.atleast_2d(D1).shape
    N2, K2 = np.atleast_2d(D2).shape
    if K1 != K2:
        raise Exception("D1 and D2 must have the same number of columns")

    D0 = vstack((D1, D2))
    a0 = mle(D0, method=method, maxiter=maxiter)
    a1 = mle(D1, method=method, maxiter=maxiter)
    a2 = mle(D2, method=method, maxiter=maxiter)

    D = 2 * (loglikelihood(D1, a1) + loglikelihood(D2, a2)
         - loglikelihood(D0, a0))
    return (D, stats.chi2.sf(D, K1), a0, a1, a2)

def test_uniform(D, do_MWM_correction=True, tol=1e-7, method='meanprecision', maxiter=None):
    '''
    Likelihood-ratio test where the null hypothesis is that all means are equal. We implement a Bartlett-type correction for small sample in Morris Water Maze data where K=4.
    
    Parameters
    ----------
    D : array of size (N,K)
        N is the number of samples and K is the number of categories. Each row of the matrices must add up to 1.
    do_MWM_correction : bool
        If True and if K=4, implement the approximate Bartlett-type correction for the likelihood-ratio statistic.
    method : string
        One of ``'fixedpoint'`` and ``'meanprecision'``, designates method by
        which to find MLE Dirichlet distribution. Default is
        ``'meanprecision'``, which is faster.
    maxiter : int
        Maximum number of iterations to take calculations. Default is
        ``sys.maxint``.

    Returns
    -------
    lr : float
        Test statistic, which is ``-2 * log`` of likelihood ratios, corrected if K=4
        and do_MWM_correction is True.
    p : float
        p-value of test.
    a0 : array
    a1 : array
        MLE parameters for the Dirichlet distributions fit to the data under the null
        hypothesis H0 (all parameters of the Dirichlet distribution are equal) and
        the alternative hypothesis (no constraint).
    '''
    
        
    N, K = D.shape
    logp = log(D).mean(axis=0)
    
    a0 = ((K-1.)/(K**2 * var(D))-1.)/K
    a0 = a0 * ones(K)
    
    # print "Initial guess ", a0
    a0 = _fit_s(D, a0, logp, tol=tol)
    
    a1 = mle(D, method=method, maxiter=maxiter)
    
    lr = 2 * (loglikelihood(D, a1) - loglikelihood(D, a0))
    
    if do_MWM_correction :
        if not K == 4:
            raise Exception("The Bartlett correction is only implemented for K=4.")
        lr *= 3./(3.+5.9*N**(-1.4))
    
    return (lr, stats.chi2.sf(lr, K-1), a1, a0)

def plot(data, label=None, do_test_uniform=True, do_MWM_correction=True, save_figure=None):
    '''
    Plot constant-sum data as stacked bars.
    
    Parameters
    ----------
    data : array of size (N,K)
        N is the number of samples and K is the number of categories. Each row of the matrices must add up to 1.
    label : string
        Name of the group.
    do_test_uniform : bool
        Whether to perform the unifornity test and show its result (p-value and parameters with approximate error bars).
    do_MWM_correction : bool
        Whether to correct the uniformity test statistic for small sample with MWM data.
    save_figure : string
        If not None, path and name of the figure file to be saved.
    
    '''
    
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import scipy
    
    def p2star(p):
        bins = [0., 0.001, 0.01, 0.05, 1.]
        i = np.digitize(p, bins[::-1])-1
        s = ''
        if i==0:
            return '(-)'
        for _ in range(i):
            s += '*'
        return '('+s+')'
    
    n, k = data.shape
    
    # Create the figure
    colors = mpl.cm.YlOrRd_r(np.linspace(0.1,.8,k))
    fig = plt.figure(figsize=(2,3))
    fig.subplots_adjust(left=.1, right=.8, bottom=0.05, top=1.)
    
    pval_cor = []
    
    # Perform the uniformity test
    if do_test_uniform:
        tu = test_uniform(data, do_MWM_correction=do_MWM_correction)
        amle = tu[2]
        pval = tu[1]
        # Lower bounds on error bars from inverse Fisher information
        vararray = 1./((_trigamma(amle)-_trigamma(amle.sum()))*n)

    ind = np.linspace(-.5,+.5,n)
    w = (ind[1]-ind[0])*0.9
    plt.bar(ind, data[:,0], w, color=colors[0])
    for j in range(1,k):
        plt.bar(ind, data[:,j], w, bottom=np.sum(data[:,:j], axis=1), color=colors[j])
        if do_test_uniform:
            plt.hlines(np.sum(amle[:j]/np.sum(amle)), ind[0]-w, ind[-1]+w, color='k', linestyles='--')
            plt.errorbar(0, np.sum(amle[:j]/np.sum(amle)), yerr=np.sqrt(vararray[j])/np.sum(amle), color='k', capsize=3)

    plt.ylim(-1e-2,1.2)
    plt.xlim(-.5-w,.5+w)
    
    # Add grid in background
    for y in np.linspace(0,1,k+1):
        plt.axhline(y=y, c='0.8', ls='--', lw=1, zorder=-32)
        plt.text(ind[-1]+1.2*w, y, '%.2g'%y, color='0.8', va='center')
    
    if label is not None:
        plt.text(0, 1.12, label, va='center', ha='center', fontsize=12)

    if do_test_uniform:
        plt.text(0, 1.05, ('$p_{\\rm Dirichlet}=%.2g$ '%pval)+p2star(pval), va='center', ha='center', fontsize=8)
        
    plt.axis('off')
        
    if save_figure is not None:
        plt.savefig(save_figure, dpi=300)

def pdf(alphas):
    '''Returns a Dirichlet PDF function'''
    alphap = alphas - 1
    c = np.exp(gammaln(alphas.sum()) - gammaln(alphas).sum())
    def dirichlet(xs):
        '''N x K array'''
        return c * (xs**alphap).prod(axis=1)
    return dirichlet

def meanprecision(a):
    '''Mean and precision of Dirichlet distribution.

    Parameters
    ----------
    a : array
        Parameters of Dirichlet distribution.

    Returns
    -------
    mean : array
        Numbers [0,1] of the means of the Dirichlet distribution.
    precision : float
        Precision or concentration parameter of the Dirichlet distribution.'''

    s = a.sum()
    m = a / s
    return (m,s)

def loglikelihood(D, a):
    '''Compute log likelihood of Dirichlet distribution, i.e. log p(D|a).

    Parameters
    ----------
    D : 2D array
        where ``N`` is the number of observations, ``K`` is the number of
        parameters for the Dirichlet distribution.
    a : array
        Parameters for the Dirichlet distribution.

    Returns
    -------
    logl : float
        The log likelihood of the Dirichlet distribution'''
    N, K = D.shape
    logp = log(D).mean(axis=0)
    return N*(gammaln(a.sum()) - gammaln(a).sum() + ((a - 1)*logp).sum())

def mle(D, tol=1e-7, method='meanprecision', maxiter=None):
    '''Iteratively computes maximum likelihood Dirichlet distribution
    for an observed data set, i.e. a for which log p(D|a) is maximum.

    Parameters
    ----------
    D : 2D array
        ``N x K`` array of numbers from [0,1] where ``N`` is the number of
        observations, ``K`` is the number of parameters for the Dirichlet
        distribution.
    tol : float
        If Euclidean distance between successive parameter arrays is less than
        ``tol``, calculation is taken to have converged.
    method : string
        One of ``'fixedpoint'`` and ``'meanprecision'``, designates method by
        which to find MLE Dirichlet distribution. Default is
        ``'meanprecision'``, which is faster.
    maxiter : int
        Maximum number of iterations to take calculations. Default is
        ``sys.maxint``.

    Returns
    -------
    a : array
        Maximum likelihood parameters for Dirichlet distribution.'''

    if method == 'meanprecision':
        return _meanprecision(D, tol=tol, maxiter=maxiter)
    else:
        return _fixedpoint(D, tol=tol, maxiter=maxiter)

def _fixedpoint(D, tol=1e-7, maxiter=None):
    '''Simple fixed point iteration method for MLE of Dirichlet distribution'''
    N, K = D.shape
    logp = log(D).mean(axis=0)
    a0 = _init_a(D)

    # Start updating
    if maxiter is None:
        maxiter = MAXINT
    for i in xrange(maxiter):
        a1 = _ipsi(psi(a0.sum()) + logp)
        # if norm(a1-a0) < tol:
        if abs(loglikelihood(D, a1)-loglikelihood(D, a0)) < tol: # much faster
            return a1
        a0 = a1
    raise Exception('Failed to converge after {} iterations, values are {}.'
                    .format(maxiter, a1))

def _meanprecision(D, tol=1e-7, maxiter=None):
    '''Mean and precision alternating method for MLE of Dirichlet
    distribution'''
    N, K = D.shape
    logp = log(D).mean(axis=0)
    a0 = _init_a(D)
    s0 = a0.sum()
    if s0 < 0:
        a0 = a0/s0
        s0 = 1
    elif s0 == 0:
        a0 = ones(a.shape) / len(a)
        s0 = 1
    m0 = a0/s0

    # Start updating
    if maxiter is None:
        maxiter = MAXINT
    for i in xrange(maxiter):
        a1 = _fit_s(D, a0, logp, tol=tol)
        s1 = sum(a1)
        a1 = _fit_m(D, a1, logp, tol=tol)
        m = a1/s1
        # if norm(a1-a0) < tol:
        if abs(loglikelihood(D, a1)-loglikelihood(D, a0)) < tol: # much faster
            return a1
        a0 = a1
    raise Exception('Failed to converge after {} iterations, values are {}.'
                    .format(maxiter, a1))

def _fit_s(D, a0, logp, tol=1e-7, maxiter=1000):
    '''Assuming a fixed mean for Dirichlet distribution, maximize likelihood
    for preicision a.k.a. s'''
    N, K = D.shape
    s1 = a0.sum()
    m = a0 / s1
    mlogp = (m*logp).sum()
    for i in xrange(maxiter):
        s0 = s1
        g = psi(s1) - (m*psi(s1*m)).sum() + mlogp
        h = _trigamma(s1) - ((m**2)*_trigamma(s1*m)).sum()

        if g + s1 * h < 0:
            s1 = 1/(1/s0 + g/h/(s0**2))
        if s1 <= 0:
            s1 = s0 * exp(-g/(s0*h + g)) # Newton on log s
        if s1 <= 0:
            s1 = 1/(1/s0 + g/((s0**2)*h + 2*s0*g)) # Newton on 1/s
        if s1 <= 0:
            s1 = s0 - g/h # Newton
        if s1 <= 0:
            raise Exception('Unable to update s from {}'.format(s0))

        a = s1 * m
        if abs(s1 - s0) < tol:
            return a

    raise Exception('Failed to converge after {} iterations, s is {}'
            .format(maxiter, s1))

def _fit_m(D, a0, logp, tol=1e-7, maxiter=1000):
    '''With fixed precision s, maximize mean m'''
    N,K = D.shape
    s = a0.sum()

    for i in xrange(maxiter):
        m = a0 / s
        a1 = _ipsi(logp + (m*(psi(a0) - logp)).sum())
        a1 = a1/a1.sum() * s

        if norm(a1 - a0) < tol:
            return a1
        a0 = a1

    raise Exception('Failed to converge after {} iterations, s is {}'
            .format(maxiter, s))

def _piecewise(x, condlist, funclist, *args, **kw):
    '''Fixed version of numpy.piecewise for 0-d arrays'''
    x = asanyarray(x)
    n2 = len(funclist)
    if isscalar(condlist) or \
            (isinstance(condlist, np.ndarray) and condlist.ndim == 0) or \
            (x.ndim > 0 and condlist[0].ndim == 0):
        condlist = [condlist]
    condlist = [asarray(c, dtype=bool) for c in condlist]
    n = len(condlist)

    zerod = False
    # This is a hack to work around problems with NumPy's
    #  handling of 0-d arrays and boolean indexing with
    #  numpy.bool_ scalars
    if x.ndim == 0:
        x = x[None]
        zerod = True
        newcondlist = []
        for k in range(n):
            if condlist[k].ndim == 0:
                condition = condlist[k][None]
            else:
                condition = condlist[k]
            newcondlist.append(condition)
        condlist = newcondlist

    if n == n2-1:  # compute the "otherwise" condition.
        totlist = condlist[0]
        for k in range(1, n):
            totlist |= condlist[k]
        condlist.append(~totlist)
        n += 1
    if (n != n2):
        raise ValueError(
                "function list and condition list must be the same")

    y = zeros(x.shape, x.dtype)
    for k in range(n):
        item = funclist[k]
        if not callable(item):
            y[condlist[k]] = item
        else:
            vals = x[condlist[k]]
            if vals.size > 0:
                y[condlist[k]] = item(vals, *args, **kw)
    if zerod:
        y = y.squeeze()
    return y

def _init_a(D):
    '''Initial guess for Dirichlet alpha parameters given data D'''
    E = D.mean(axis=0)
    E2 = (D**2).mean(axis=0)
    return ((E[0] - E2[0])/(E2[0]-E[0]**2)) * E

def _ipsi(y, tol=1.48e-9, maxiter=10):
    '''Inverse of psi (digamma) using Newton's method. For the purposes
    of Dirichlet MLE, since the parameters a[i] must always
    satisfy a > 0, we define ipsi :: R -> (0,inf).'''
    y = asanyarray(y, dtype='float')
    x0 = _piecewise(y, [y >= -2.22, y < -2.22],
            [(lambda x: exp(x) + 0.5), (lambda x: -1/(x+euler))])
    for i in xrange(maxiter):
        x1 = x0 - (psi(x0) - y)/_trigamma(x0)
        if norm(x1 - x0) < tol:
            return x1
        x0 = x1
    raise Exception(
        'Unable to converge in {} iterations, value is {}'.format(maxiter, x1))

def _trigamma(x):
    return polygamma(1, x)
