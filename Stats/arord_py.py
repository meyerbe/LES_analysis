
import numpy

'''
%ARORD	Evaluates criteria for selecting the order of an AR model.
%
%  [SBC,FPE]=ARORD(R,m,mcor,ne,pmin,pmax) returns approximate values
%  of the order selection criteria SBC and FPE for models of order
%  pmin:pmax. The input matrix R is the upper triangular factor in the
%  QR factorization of the AR model; m is the dimension of the state
%  vectors; the flag mcor indicates whether or not an intercept vector
%  is being fitted; and ne is the number of block equations of size m
%  used in the estimation. The returned values of the order selection
%  criteria are approximate in that in evaluating a selection
%  criterion for an AR model of order p < pmax, pmax-p initial values
%  of the given time series are ignored.
%
%  ARORD is called by ARFIT. 
%	
%  See also ARFIT, ARQR.

%  For testing purposes, ARORD also returns the vectors logdp and np,
%  containing the logarithms of the determinants of the (scaled)
%  covariance matrix estimates and the number of parameter vectors at
%  each order pmin:pmax.

%  Modified 17-Dec-99
%  Author: Tapio Schneider
%          tapio@gps.caltech.edu
'''

# function [sbc, fpe, logdp, np] = arord(R, m, mcor, ne, pmin, pmax)
def arord(R, m, mcor, ne, pmin, pmax):
    print('calling arord')
    imax 	  = pmax-pmin        # maximum index of output vectors

    # initialize output vectors
    sbc     = numpy.zeros((imax+1),dtype=int)        # Schwarz's Bayesian Criterion
    fpe     = numpy.zeros((imax+1),dtype=int)        # log of Akaike's Final Prediction Error
    logdp   = numpy.zeros((imax+1),dtype=int)        # determinant of (scaled) covariance matrix
    mp      = numpy.zeros((imax+1),dtype=int)        # number of parameter vectors of length m
    mp[imax]= numpy.int(m*pmax+mcor)

    # Get lower right triangle R22 of R:
    #   | R11  R12 |
    # R=|          |
    #   | 0    R22 |

    print(mp[imax], m)
    print('R',R.shape, ne)
    R22     = R[ mp[imax] : mp[imax]+m, mp[imax] : mp[imax]+m]
    print('')
    print('')
    # print(R)
    print(R[-1,-1])
    print('R22',R22)
    print('')
    print('')

    # From R22, get inverse of residual cross-product matrix for model of order pmax
    invR22  = numpy.linalg.inv(R22)
    Mp      = invR22*invR22.T

    # For order selection, get determinant of residual cross-product matrix
    # logdp = log det(residual cross-product matrix)
    logdp[imax] = 2.*numpy.log(numpy.abs(numpy.prod(numpy.diag(R22))));

    # Compute approximate order selection criteria for models of order pmin:pmax
    # i = imax
    # for p = pmax:-1:pmin
    for p in range(pmax,pmin-1,-1):
        print('p',p, pmax, pmin)
        # mp[imax] = m*p + mcor	# number of parameter vectors of length m
        if p < pmax:
        # Downdate determinant of residual cross-product matrix
        # Rp: Part of R to be added to Cholesky factor of covariance matrix
        #     Rp = R(np(i)+1:np(i)+m, np(imax)+1:np(imax)+m);
            Rp = R[mp[imax]+1:mp[imax] + m, mp[imax] + 1:mp[imax] + m]
            print(Rp.shape)

            # Get Mp, the downdated inverse of the residual cross-product matrix, using the Woodbury formula
            # L = chol(eye(m) + Rp*Mp*Rp.T).T
            L = numpy.linalg.cholesky(numpy.identity(m) + Rp * Mp * Rp.T).T
            # N = L \ Rp * Mp
            N = numpy.divide(L,Rp) * Mp
            Mp = Mp - N.T*N

            # Get downdated logarithm of determinant
            # logdp(i) = logdp(i+1) + 2.* log(abs(prod(diag(L))))
            logdp[i] = logdp[i + 1] + 2. * numpy.log(numpy.abs(numpy.prod(numpy.diag(L))))
        #     # end
        #
        # # Schwarz's Bayesian Criterion
        # sbc(i) = logdp(i)/m - log(ne) * (ne-np(i))/ne
        #
        # #  logarithm of Akaike's Final Prediction Error
        # fpe(i) = logdp(i)/m - log(ne*(ne-np(i))/(ne+np(i)))
        #
        # # Modified Schwarz criterion (MSC):
        # # msc(i) = logdp(i)/m - (log(ne) - 2.5) * (1 - 2.5*np(i)/(ne-np(i)));
        #
        # i      = i-1;                # go to next lower order
        # # end

    return

