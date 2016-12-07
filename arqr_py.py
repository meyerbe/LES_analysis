import numpy

# function [R, scale]=arqr(v, p, mcor)
'''
%ARQR	QR factorization for least squares estimation of AR model.
%
%  [R, SCALE]=ARQR(v,p,mcor) computes the QR factorization needed in
%  the least squares estimation of parameters of an AR(p) model. If
%  the input flag mcor equals one, a vector of intercept terms is
%  being fitted. If mcor equals zero, the process v is assumed to have
%  mean zero. The output argument R is the upper triangular matrix
%  appearing in the QR factorization of the AR model, and SCALE is a
%  vector of scaling factors used to regularize the QR factorization.
%
%  ARQR is called by ARFIT. 
%
%  See also ARFIT.

%  Modified 29-Dec-99
%  Author: Tapio Schneider
%          tapio@gps.caltech.edu
'''
def arqr(v, p, mcor):
    # v_i: state vector at given time step (one observation)
    # m: dimension of state vectors v_i
    # n: number of time steps / observations
    # v: matrix with state vectors v_i (i=1..n)
    [n,m] = v.shape
    # u_i: predictor vector (i=1..n)
    # np: dimension of predictor vectors u_i
    # ne: number of block equations
    mp    = m*p + mcor
    # ne    = n-p                     # number of block equations of size m

  # % If the intercept vector w is to be fitted, least squares (LS)
  # % estimation proceeds by solving the normal equations for the linear
  # % regression model
  # %
  # %                  v(k,:)' = Aaug*u(k,:)' + noise(C)
  # %
  # % with Aaug=[w A] and `predictors'
  # %
  # %              u(k,:) = [1 v(k-1,:) ...  v(k-p,:)].
  # %
  # % If the process mean is taken to be zero, the augmented coefficient
  # % matrix is Aaug=A, and the regression model
  # %
  # %                u(k,:) = [v(k-1,:) ...  v(k-p,:)]
  # %
  # % is fitted.
  # % The number np is the dimension of the `predictors' u(k).



    K = numpy.zeros(shape=(n, mp+m))
    if (mcor == 1):
        K[:,0] = numpy.ones((n,))

    print('mcor', mcor)
    print(K)

    print(numpy.mod(5,3))
    for i in range(n):
        for j in range(mcor,mp):
            b = numpy.mod((j-mcor),mp)
            a = ((j-mcor)-b)/mp
            print(i,j,'a',a,b,n,m,mp)

    for i in range(n):
        for j in range(1,mp):
            for ia in range(p-1):
                ib = j - (1+ia*m)
                print(i, j, ia, ib, m)
                if ib > 0:
                    print(i, j, ia, ib, m)
                    K[i,j] = v[i-ia,ib]
        for j in range(m):
            K[i,j+mp] = v[i,j]


    # for a in range(p):
    #     for b in range(m):
    #         K[i,j] = v[]
    #         K[i,mp:(m+mp)] = v[i-a,]

  #
  #   # Assemble the data matrix K (of which a QR factorization will be computed)
  #   K = numpy.zeros(shape=(ne,np+m))               # initialize K
  #   if (mcor == 1):
  #   # first column of K consists of ones for estimation of intercept vector w
  #       K[:,1] = numpy.ones(ne,1)
  #   print('...',v.shape,K.shape)
  #
  #   # Assemble `predictors' u in K
  #   for j in range(1,p):
  #       # K(:, mcor+m*(j-1)+1:mcor+m*j) = [v(p-j+1:n-j, :)];
  #       K[:,mcor+m * (j-1) + 1:mcor + m * j] = v[p - j + 1:n - j,:]
  #
  #   # Add `observations' v (left hand side of regression model) to K
  #   K(:,np+1:np+m) = [v(p+1:n, :)];
  #
  # # % Compute regularized QR factorization of K: The regularization
  # # % parameter delta is chosen according to Higham's (1996) Theorem
  # # % 10.7 on the stability of a Cholesky factorization. Replace the
  # # % regularization parameter delta below by a parameter that depends
  # # % on the observational error if the observational error dominates
  # # % the rounding error (cf. Neumaier, A. and T. Schneider, 2001:
  # # % "Estimation of parameters and eigenmodes of multivariate
  # # % autoregressive models", ACM Trans. Math. Softw., 27, 27--57.).
  #   q     = np + m             # number of columns of K
  #   delta = (q^2 + q + 1)*eps  # Higham's choice for a Cholesky factorization
  #   scale = sqrt(delta)*sqrt(sum(K.^2))
  #   R     = triu(qr([K; diag(scale)]))
  #
  # # return [R, scale]
    return