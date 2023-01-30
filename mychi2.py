import configparser
import numpy as np
from scipy.interpolate import interp1d

from mydata import get_index

def fwd_subst(R, b):
    '''Forward substitution.'''
    n = b.size
    x = np.zeros(n)
    for i in range(n):
        x[i] = b[i]
        for j in range(i):
            x[i] -= R[j, i] * x[j]
        x[i] /= R[i, i]
    return x

def bwd_subst(R, b):
    '''Backward substitution.'''
    n = b.size
    x = np.zeros(n)
    for i in range(n-1, -1, -1):
        x[i] = b[i]
        for j in range(i+1, n):
            x[i] = x[i] - R[i][j] * x[j]
        x[i] /= R[i][i]
    return x

def init_lstsq(sd, Rcov, npoly, imin, imax):
    '''Initialize the least square fitting for nuisance parameters.
    Refs:
      Numerical Recipes 3rd edition, section 15.4.
      (But the matrix manipulations are simplified here.)
    Arguments:
      sd: the s array for data;
      Rcov: the upper triangular matrix from QR decomp of cov;
      npoly: the number of nuisance parameters;
      imin: the minimum index for the fitting range;
      imax: the maximum index for the fitting range.
    Return: [basis function, design matrix, RHS matrix for LS].'''
    # The basis function and design maxtrix
    nd = sd.size
    basis = np.zeros([npoly, nd])
    A = np.zeros([npoly, nd])
    for i in range(npoly):
        basis[i, imin:imax] = sd[imin:imax]**(i - 2)
        A[i, :] = fwd_subst(Rcov, basis[i, :])

    # Construct U from the QR decomposition of A
    U = np.linalg.qr(np.transpose(A), mode='r')
    # Construct M = U^{-T} . A^T
    M = np.zeros([npoly, nd])
    col = np.zeros(npoly)
    for i in range(nd):
        col = fwd_subst(U, A[:, i])
        M[:, i] = col
    # Construct M . R^{-T}
    for i in range(npoly):
        M[i, :] = bwd_subst(Rcov, M[i, :])

    return [basis, U, M]

class Chi2Class():
    def __init__(self, xim_var, xid_var, covmat_var):

        self.covmat = covmat_var
        self.xim = xim_var
        self.xid = xid_var
        
        self.list_params = self.xim.list_params
        self.chi2_norm = (self.covmat.nmock - self.xid.nidx - 2) / (self.covmat.nmock - 1)
      
        print('STATUS: Initialize matrices for least square fitting of the nuisance params.')
        self.basis, self.A, self.M = init_lstsq(self.xid.sd, self.covmat.Rcov, self.xim.npoly, self.xid.imin, self.xid.imax)

    def chi2_func(self, alpha, params):
        '''Define the (log) likelihood.'''
        B = params[0]
        a0 = params[2]
        a1 = params[3]
        a2 = params[4]
        

        imin = self.xid.imin
        imax = self.xid.imax
        # Compute the model 2PCF with a given alpha
        fxim = interp1d(self.xim.sm, self.xim.xi_model(params), kind='cubic')
        xi = fxim(self.xid.sd[imin:imax] * alpha)

        # Least square fitting of nuisance parameters
        dxi = self.xid.xid[imin:imax] - (xi * (B**2) + (a0 * np.ones(len(self.xid.sd[imin:imax]))) + (a1  / self.xid.sd[imin:imax]) + (a2 / (self.xid.sd[imin:imax]**2)))
        
        
        # Compute chi-squared
        diff = np.zeros(self.xid.ndbin)
        diff[imin:imax] = dxi
        
        matrix = self.covmat.icov
        chisq = np.dot(diff.T, np.dot(matrix, diff))
        return chisq * self.chi2_norm * self.covmat.rescale_chi2_by

    def best_fit(self, alpha, params):
        '''Compute the best-fit theoretical Xi curve.'''
        B = params[0]
        imin = self.xid.imin
        imax = self.xid.imax
        sm = self.xim.sm

        # Compute the model 2PCF with a given alpha
        fxim = interp1d(sm, self.xim.xi_model(params), kind='cubic')
        fxim_nw = interp1d(sm, self.xim.xi_model_nw(params), kind='cubic')
        xi = fxim(self.xid.sd[imin:imax] * alpha)
        xi_nw = fxim_nw(self.xid.sd[imin:imax] * alpha)

        # Least square fitting of nuisance parameters
        dxi = self.xid.xid[imin:imax] - xi * (B**2)
        poly = np.dot(self.M[:, imin:imax], dxi)
        a_poly = bwd_subst(self.A, poly)

        # Compute best-fit
        if alpha >= 1:
            imin0, imax0, nidx0 = get_index(sm, sm[1], sm[-2] / alpha)
        else:
            imin0, imax0, nidx0 = get_index(sm, sm[1] / alpha, sm[-2])
        
        best = B**2 * fxim(sm[imin0:imax0] * alpha)
        best_nw = B**2 * fxim_nw(sm[imin0:imax0] * alpha)
        best_broadband = np.zeros(nidx0)
        for i in range(self.xim.npoly):
            best += a_poly[i] * sm[imin0:imax0]**(i - 2)
            best_nw += a_poly[i] * sm[imin0:imax0]**(i - 2)
            best_broadband += a_poly[i] * sm[imin0:imax0]**(i - 2)
        return [sm[imin0:imax0], best, best_nw, best_broadband]
