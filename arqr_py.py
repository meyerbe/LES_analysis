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
    # v_i: state vector at given time step (one observation), dim(v_i) = m
    # m: dimension of state vectors v_i
    # n: number of time steps / observations
    # v: matrix with state vectors v_i (i=1..n), dim(v)=[n,m]
    [n,m] = v.shape
    # u_i: predictor vector (i=1..n)
    # np: dimension of predictor vectors u_i
    # ne: number of block equations
    mp = m*p + mcor
    ne    = n-p                     # number of block equations of size m

    # K: matrix that has to be QR decomposed, dim(K) = [n,mp+m] = [n,m*(p+1)+mcor]
    K = numpy.zeros(shape=(ne, mp + m))

    # If the intercept vector w is to be fitted, least squares (LS)
    # estimation proceeds by solving the normal equations for the linear
    # regression model
    #                  v(k,:)' = Aaug*u(k,:)' + noise(C)
    #
    # with Aaug=[w A] and `predictors'
    #              u(k,:) = [1 v(k-1,:) ...  v(k-p,:)].
    #
    # If the process mean is taken to be zero, the augmented coefficient
    # matrix is Aaug=A, and the regression model
    #                u(k,:) = [v(k-1,:) ...  v(k-p,:)]
    # is fitted.
    # The number np is the dimension of the `predictors' u(k).


    if (mcor == 1):
        K[:,0] = numpy.ones((ne,))

    print('mcor', mcor)
    print('p', p)
    # print(K)
    print('K', K.shape, n, ne)
    print('v', v.shape)

    for j in range(p):
        for k in range(n - p):
        # for k in range(1,n-p+1):
            for l in range(m):
            # for l in range(1,m+1):
                # Delta in K = m        --> dim(v_i) = m
                # Delta in v = n - p    --> dim(K) = [n-p,mp]
                K[k, mcor + m * j + (l + 1)] = v[p - j + k, l]
                # K[k,mcor+m*j+(l+1)] = v[p-(j+1)+(k+1),l]
        # K[:, mcor + m*j + 1 : mcor + m*(j+1)] = v[p - (j+1) + 1:n - (j+1),:]
        # K(:, mcor + m * (j - 1) + 1:mcor + m * j) = v(p - j + 1:n - j,:)


        # Add 'observations' v (left hand side of regression model) to K
    for k in range(n-p):
        for l in range(m):
            K[k, mp + l] = v[p + k, l]
    # K[:, mp + 1:mp + m] = v[p + 1:n,:]
    # K(:, mp + 1:mp + m) = v(p + 1:n,:)


    # % Compute regularized QR factorization of K: The regularization
    # % parameter delta is chosen according to Higham's (1996) Theorem
    # % 10.7 on the stability of a Cholesky factorization. Replace the
    # % regularization parameter delta below by a parameter that depends
    # % on the observational error if the observational error dominates
    # % the rounding error (cf. Neumaier, A. and T. Schneider, 2001:
    # % "Estimation of parameters and eigenmodes of multivariate
    # % autoregressive models", ACM Trans. Math. Softw., 27, 27--57.).
    q     = mp + m             # number of columns of K
    # print(numpy.spacing)
    # delta = (q^2 + q + 1)*eps  # Higham's choice for a Cholesky factorization
    delta = (q**2 + q + 1) * 1e-5 # Higham's choice for a Cholesky factorization
    # scale = sqrt(delta)*sqrt(sum(K.^2))
    scale = numpy.sqrt(delta) * numpy.sqrt(sum(K**2))
    # R     = triu(qr([K; diag(scale)]))
    # R = numpy.triu(qr([K;diag(scale)]))
    Q, R = numpy.linalg.qr(K)
    print('')
    print('')
    # print(Q)
    print(R)
    print('')
    print('')

    return [R, scale]
    # return