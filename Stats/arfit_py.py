import sys
from arqr_py import arqr
from arord_py import arord

'''
%ARFIT	Stepwise least squares estimation of multivariate AR model.
%
%  [w,A,C,SBC,FPE,th]=ARFIT(v,pmin,pmax) produces estimates of the
%  parameters of a multivariate AR model of order p,
%
%      v(k,:)' = w' + A1*v(k-1,:)' +...+ Ap*v(k-p,:)' + noise(C),
%
%  where p lies between pmin and pmax and is chosen as the optimizer
%  of Schwarz's Bayesian Criterion. The input matrix v must contain
%  the time series data, with columns of v representing variables
%  and rows of v representing observations.  ARFIT returns least
%  squares estimates of the intercept vector w, of the coefficient
%  matrices A1,...,Ap (as A=[A1 ... Ap]), and of the noise covariance
%  matrix C.
%
%  As order selection criteria, ARFIT computes approximations to
%  Schwarz's Bayesian Criterion and to the logarithm of Akaike's Final
%  Prediction Error. The order selection criteria for models of order
%  pmin:pmax are returned as the vectors SBC and FPE.
%
%  The matrix th contains information needed for the computation of
%  confidence intervals. ARMODE and ARCONF require th as input
%  arguments.
%
%  If the optional argument SELECTOR is included in the function call,
%  as in ARFIT(v,pmin,pmax,SELECTOR), SELECTOR is used as the order
%  selection criterion in determining the optimum model order. The
%  three letter string SELECTOR must have one of the two values 'sbc'
%  or 'fpe'. (By default, Schwarz's criterion SBC is used.) If the
%  bounds pmin and pmax coincide, the order of the estimated model
%  is p=pmin=pmax.
%
%  If the function call contains the optional argument 'zero' as the
%  fourth or fifth argument, a model of the form
%
%         v(k,:)' = A1*v(k-1,:)' +...+ Ap*v(k-p,:)' + noise(C)
%
%  is fitted to the time series data. That is, the intercept vector w
%  is taken to be zero, which amounts to assuming that the AR(p)
%  process has zero mean.

%  Modified 14-Oct-00
%  Authors: Tapio Schneider
%           tapio@gps.caltech.edu
%
%           Arnold Neumaier
%           neum@cma.univie.ac.at
'''

# def function [w, A, C, sbc, fpe, th]=arfit(v, pmin, pmax, selector, no_const):
def arfit(v, pmin, pmax, selector='sbc',no_const=1):
    # n: number of observations; m: dimension of state vectors
    [n,m]   = v.shape
    print('n,m', n,m)

    if (pmin != round(pmin) or pmax != round(pmax)):
        print('Order must be integer.')
        sys.exit()
    if (pmax < pmin):
        print('PMAX must be greater than or equal to PMIN.')

    print('test succesful')

    # set defaults and check for optional arguments
    if selector == 'zero':
        mcor = 0            # no intercept vector to be fitted
        selector = 'sbc'
    elif selector == 'sbc':
        mcor = 1            # fit intercept vector
    else:
        mcor = 1            # fit intercept vector
    print('selector: ' + selector + ', mcor: ' +str(mcor))

    # # !!!
    # if no_const == 0:
    #     mcor = 0            # no intercept vector to be fitted
    # else:
    #     print('Bad argument. Usage: ')
    #     sys.exit()
    # # !!!

    ne  	= n - pmax;               # number of block equations of size m
    mpmax	= m*pmax + mcor;          # maximum number of parameter vectors of length m
    if (ne <= mpmax):
        print('Time series too short.')
        sys.exit()

    # (1) compute QR factorization for model of order pmax
    [R, scale] = arqr(v, pmax, mcor)
    # print('R', R)

    # (2) compute approximate order selection criteria for models
    # of order pmin:pmax
    # [sbc, fpe]   = arord(R, m, mcor, ne, pmin, pmax);
    arord(R, m, mcor, ne, pmin, pmax)

    # get index iopt of order that minimizes the order selection
    # criterion specified by the variable selector
    # [val, iopt]  = min(eval(selector));

  # % select order of model
  # popt         = pmin + iopt-1; % estimated optimum order
  # np           = m*popt + mcor; % number of parameter vectors of length m
  #
  # % decompose R for the optimal model order popt according to
  # %
  # %   | R11  R12 |
  # % R=|          |
  # %   | 0    R22 |
  # %
  # R11   = R(1:np, 1:np);
  # R12   = R(1:np, mpmax+1:mpmax+m);
  # R22   = R(np+1:mpmax+m, mpmax+1:mpmax+m);
  #
  # % get augmented parameter matrix Aaug=[w A] if mcor=1 and Aaug=A if mcor=0
  # if (np > 0)
  #   if (mcor == 1)
  #     % improve condition of R11 by re-scaling first column
  #     con 	= max(scale(2:mpmax+m)) / scale(1);
  #     R11(:,1)	= R11(:,1)*con;
  #   end;
  #   Aaug = (R11\R12)';
  #
  #   %  return coefficient matrix A and intercept vector w separately
  #   if (mcor == 1)
  #     % intercept vector w is first column of Aaug, rest of Aaug is
  #     % coefficient matrix A
  #     w = Aaug(:,1)*con;        % undo condition-improving scaling
  #     A = Aaug(:,2:np);
  #   else
  #     % return an intercept vector of zeros
  #     w = zeros(m,1);
  #     A = Aaug;
  #   end
  # else
  #   % no parameters have been estimated
  #   % => return only covariance matrix estimate and order selection
  #   % criteria for ``zeroth order model''
  #   w   = zeros(m,1);
  #   A   = [];
  # end
  #
  # % return covariance matrix
  # dof   = ne-np;                % number of block degrees of freedom
  # C     = R22'*R22./dof;        % bias-corrected estimate of covariance matrix
  #
  # % for later computation of confidence intervals return in th:
  # % (i)  the inverse of U=R11'*R11, which appears in the asymptotic
  # %      covariance matrix of the least squares estimator
  # % (ii) the number of degrees of freedom of the residual covariance matrix
  # invR11 = inv(R11);
  # if (mcor == 1)
  #   % undo condition improving scaling
  #   invR11(1, :) = invR11(1, :) * con;
  # end
  # Uinv   = invR11*invR11';
  # th     = [dof zeros(1,size(Uinv,2)-1); Uinv];



  # return [w, A, C, sbc, fpe, th]
    return




# ____________
