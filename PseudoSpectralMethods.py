import math
import numpy as np
from scipy.linalg import toeplitz
from scipy.fftpack import ifft
from numpy.linalg import inv
from numpy.linalg import pinv
import Tests

def chebdiff(N, M):
    
    """
    The function computes the differentiation
    matrices D1, D2, ..., DM on Chebyshev nodes.
    
    The code implements two strategies for enhanced accuracy suggested by
    W. Don and S. Solomonoff in SIAM J. Sci. Comp. Vol. 6, pp. 1253--1268
    (1994).
    
    The two strategies are:
      (a) the use of trigonometric identities to avoid the computation of
          differences x(k)-x(j)
    
      (b) the use of the "flipping trick" which is necessary since
          sin t can be computed to high relative precision when t is small
          whereas sin (pi-t) cannot.
    
    Note added May 2003: It may, in fact, be slightly better not to
    implement the strategies (a) and (b).
    
    Please consult the following paper for details:
    "Spectral Differencing with a Twist", by R. Baltensperger and
    M.R. Trummer, to appear in SIAM J. Sci. Comp.
    
    J.A.C. Weideman, S.C. Reddy 1998.  Help notes modified by
    JACW, May 2003.
      
      
    INPUTS:
             N:  size of required differentation matrix
             M:  highest derivative needed (integer).
                 Note:  0 < M <= N-1
    
    
    OUTPUT:
             x:  Interior Chebyshev points (length N) [Chebyshev nodes]
            DM:  differentiation matrix of size (M x N x N), containing 
                 differentiation matrices:
                     DM(0:N-1, 0:N-1, l) contains l-th derivative
                     matrix, l = 1..M
    
    """


    # Identity matrix with logical matrix
    I = np.identity(N)
    #I_logical = np.matrix(np.eye(N, dtype=bool))


    # n1 n2 are indices used for the flipping trick
    n1 = math.floor(N / 2.0)
    n2 = math.ceil(N / 2.0)

    # Get the theta vector
    k = np.arange(N)[:, np.newaxis]
    theta = k * math.pi / (N - 1.0)

    # Chebyshev points of the second kind or you can think of them as
    # extreme points on [-1,1] of T_{N-1}, the Chebyshev poolynomial of 
    # degree N-1
    x = np.sin(math.pi * np.arange(N-1, -N, -2)[:, np.newaxis] / (2.0 * (N-1)))

    T = np.tile(theta / 2.0, (1.0, N))

    # Trigonometric identity
    DX = 2.0 * np.sin(T.T + T) * np.sin(T.T - T)

    # Flipping trick
    DX = np.concatenate((DX[:n1], -1.0 * np.flipud(np.fliplr(DX[:n2]))), axis = 0)

    # Put 1s in the main diagonal of DX
    np.fill_diagonal(DX, 1.0)

    # C is the matrix with entries c(k)/c(j)
    C = toeplitz((-1.0)**k)
    C[0:1] = C[0:1] * 2.0       # multiply the first row by 2
    C[N-1:N] = C[N-1:N] * 2.0   # multiply the  last row by 2
    C[:,0] = C[:,0] / 2.0       #   divide the first col by 2
    C[:,N-1] = C[:,N-1] / 2.0   #   divide the  last col by 2

    # Z contains entries 1/(x(k)-x(j))
    Z = 1.0 / DX
    
    # with zeros on the diagonal.
    np.fill_diagonal(Z, 0.0)

    D = np.matrix(np.identity(N), copy=False)

    DM = np.empty((M, N, N))
    
    for l in range(0,M):
        A = np.tile(np.diag(D), N).reshape(N,N).T
        B = np.multiply(C, A) - D
        G = np.multiply(Z, B)
        D = (l + 1) * G
        #D1 = (l + 1) * np.multiply(Z, np.multiply(C, np.tile(np.diag(D), N).reshape(N,N).T) - D)

        #sum up the rows of D
        row_sums = -1.0 * D.sum(axis=1)
        # then assign the diagonals the value of the sum of the row.
        np.fill_diagonal(D, row_sums)
        DM[l,:,:] = D # DM[l::] = D
    
    return x, DM


def cheb4c(N, debug):
    """
      The function [x, D4] =  cheb4c(N) computes the fourth 
      derivative matrix on Chebyshev interior points, incorporating 
      the clamped boundary conditions u(1)=u'(1)=u(-1)=u'(-1)=0.
    
      The code implements two strategies for enhanced accuracy suggested by 
      W. Don and S. Solomonoff in SIAM J. Sci. Comp. Vol. 6, pp. 1253--1268 
      (1994).
      
      The two strategies are: 
        (a) the use of trigonometric 
            identities to avoid the computation of differences 
            x(k)-x(j) and 
    
        (b) the use of the "flipping trick"
            which is necessary since sin t can be computed to high
            relative precision when t is small whereas sin (pi-t) cannot.
      
      J.A.C. Weideman, S.C. Reddy 1998.
    
    INPUTS:
             N:  N-2 = Order of differentiation matrix.  
                 (The interpolant has degree N+1.)
    
    OUTPUTS:
             x:  Interior Chebyshev points (vector of length N-2)
            D4:  Fourth derivative matrix  (size (N-2)x(N-2))
    
    """ 
        
    # Identity matrix with logical matrix
    I = np.matrix(np.identity(N - 2.0), copy=False)
    I_logical = np.matrix(np.eye(N - 2.0, dtype=bool))
    
    # n1 n2 are indices used for the flipping trick
    n1 = math.floor((N / 2.0) - 1)
    n2 = math.ceil((N / 2.0) - 1)    
    
    # Get the theta vector
    k = np.arange(1, N - 1.0)[:, np.newaxis]
    theta = k * math.pi / (N - 1.0)
    
    # Chebyshev points
    x = np.sin(math.pi * np.arange(N - 3.0, 1.0 - N, -2.0)[:, np.newaxis] / (2.0 * (N - 1.0)))
    
    # s = sin(theta)
    s = np.concatenate((np.sin(theta[:n1]), np.flipud(np.sin(theta[:n2]))), axis = 0)
    
    # Compute weight function and its derivatives
    alpha = s**4.0
    beta1 = -4.0 * (s**2.0) * (x / alpha)
    beta2 = 4.0 * (3.0 * x**2.0 - 1.0) / alpha
    beta3 = 24.0 * x / alpha
    beta4 = 24.0 / alpha
    B = np.array([beta1.T, beta2.T, beta3.T, beta4.T])

    T = np.tile(theta / 2.0, (1.0, N - 2.0))

    # Trigonometric identity
    DX = 2.0 * np.sin(T.T + T) * np.sin(T.T - T)
    
    # Flipping trick
    DX = np.concatenate((DX[:n1], -1.0 * np.flipud(np.fliplr(DX[:n2]))), axis = 0)
    
    # Put 1s in the main diagonal of DX
    np.fill_diagonal(DX, 1.0)    
    
    # Compute the matrix with entries c(k)/c(j)
    ss = (s**2.0) * ((-1.0)**k)
    S = np.tile(ss, N - 2.0)
    C = S / S.T
    
    # Z contains entries 1/(x(k)-x(j))
    Z = 1.0 / DX
    
    # with zeros on the diagonal.
    np.fill_diagonal(Z, 0.0)
    
    # X is same as Z, but with diagonal entries removed
    X = Z.T
    
    # INSTEAD OF THIS, think of the arrays as vertical lists!!!, 
    # then when you get rid of the zeroes, 
    # you can concatenate the lists together side by side.
    Y = np.zeros((N-3, N-2.0))
    for i in range(0, int(N)-2):
        Y[:,i] = [a for a in X[:,i] if a != 0.0]
    
    X = Y

    # Initialize Y and D vectors.
    # Y contains matrix of cumulative sums, D scaled differentiation matrices.
    Y = np.ones((N-3, N-2.0))
    
    D = I
    A = np.ones((N - 2.0, N - 2.0))
    DM = np.empty((4, N - 2.0, N - 2.0))
    
    for l in range(0, 4):
        
#        test_01 = B[l, :]
#        test_02 = (l + 1) * Y[0:N-3, :] * X
#        test_03 = np.concatenate((B[l, :], (l + 1) * Y[0:N-3, :] * X), axis=0)
        
        # Recursion for diagonals
        Y = np.cumsum(np.concatenate((B[l, :], (l + 1) * Y[0:N-3, :] * X), axis=0), axis=0)
        
        # Off-diagonal
        a = np.tile(np.diag(D), N - 2.0).reshape(N - 2.0, N - 2.0).T
        b = np.multiply(C, a) - D
        g = np.multiply(Z, b)
        D = (l + 1) * g
        
        # Correct the diagonal, put the last row of Y into the empty diagonal elements in D
        np.fill_diagonal(D, Y[-1, :])
        
        DM[l,:,:] = D


    D4 = DM[3,:,:]
      
    return x, D4


def clencurt(n1):
    """
    CLENCURT
    
        nodes x (Chebyshev points) and weights w for Clenshaw-Curtis
        quadrature from Spectral Methods in MATLAB
    """
    # Trefethen's N is different from W&R's
    if n1 is 1:
        x = 0
        w = 2
    else:
        n = n1 - 1
        C = np.zeros((n1, 2))
        k = 2 * (1 + np.arange(np.floor(n / 2)))
        C[::2, 0] = 2 / np.hstack((1, 1 - k * k))
        C[1, 1] = -n
        V = np.vstack((C, np.flipud(C[1:n, :])))
        F = np.real(ifft(V, n=None, axis=0))
        x = -1.0 * F[0:n1, 1]  # I added the -1 at the beginning of this operation
        w = np.hstack((F[0, 0], 2 * F[1:n, 0], F[n, 0]))

    return x, w


def chebint(fk, x, debug):
    #  The function p = chebint(fk, x) computes the polynomial interpolant
    #  of the data (xk, fk), where xk are the Chebyshev nodes.  
    #  Two or more data points are assumed.
    #
    #  Input:
    #   - fk:  Vector of y-coordinates of data, at Chebyshev points 
    #          x(k) = cos((k-1)*pi/(N-1)), k = 1...N.
    #   - x:   Vector of x-values where polynomial interpolant is to be evaluated.
    #
    #  Output:
    #   - p:    Vector of interpolated values.
    #
    #
    #  The code implements the barycentric formula; see page 252 in
    #  P. Henrici, Essentials of Numerical Analysis, Wiley, 1982.
    #  (Note that if some fk > 1/eps, with eps the machine epsilon,
    #  the value of eps in the code may have to be reduced.)
    #
    #  J.A.C. Weideman, S.C. Reddy 1998     
    
    # Write some code to make sure that thte data passed in is column vecotors
    fk = np.expand_dims(fk, axis=1)
    x = np.expand_dims(x, axis=1)

    N = len(fk)
    M = len(x)

    # Compute the chebyshev points
    xk = np.sin(math.pi * np.arange(N - 1, -1.0 * N, -2.0)[:, np.newaxis] / (2.0 * (N - 1)))
    
    powers = np.arange(N)[:, np.newaxis]
    w = np.ones((N,1)) * ((-1.0)**(powers))
    
    w[0] = w[0] / 2.0
    w[-1] = w[-1] / 2.0
    
    # Compute quantities x-x(k)
    D = np.tile(x, N) - np.tile(xk, M).T

    # Compute the reciprocals
    D = 1./(D + np.spacing(1)*(D == 0.0))

    # Evaluate interpolant as matrix-vector products
    p = np.dot(D, w * fk) / (np.dot(D,w))
    
    return p


def my_ifft(u_tilde, alpha, beta, ffg):

    u = np.zeros((ffg.Nx, 3*ffg.modes, ffg.Nz), dtype=np.complex128)
    u_tilde = np.asarray(u_tilde)
    u_tilde = u_tilde[:, 0]

    for iz in range(0, ffg.Nz):
        for ix in range(0, ffg.Nx):
            u[ix, :, iz] = u_tilde * np.exp(1j * ((alpha * ffg.x[ix]) + (beta * ffg.z[iz]) - (alpha * ffg.c * ffg.t)))

    return u


def get_state_vectors(alpha, beta, Re, N, omega, bf, vel_profile):

    #================================================================
    #### Calculate the derivative matrices in the wall-normal direction
    #================================================================
    # y_cheb_full are the interpolated y co-ordinates, 
    # i.e. Chebyshev interior points (nodes), with endpoints included.
    y_cheb_full, DM = chebdiff(N, 2)
    D1 = DM[0, 1:-1, 1:-1]              # First  derivative matrix (partial_d_dy)
    D2 = DM[1, 1:-1, 1:-1]              # Second derivative matrix (partial_d_dyy)
    y_cheb, D4 = cheb4c(N, False)    # Fourth derivative matrix (partial_d_dyyyy)
    # y_cheb has endpoints removed, i.e. y_cheb_full[1:-1].
    # Hence the number of modes are given as
    modes = N-2 # because the endpoints are removed.

    #================================================================
    #### Assign mean flow
    #================================================================
    # Check to see if a turbulent mean is given
    if len(vel_profile) == 0: # No turbulent mean given.
        # Check which baseflow to use
        if bf == 'lam': # Laminar Base flow 
            U = np.identity(modes)
            np.fill_diagonal(U, 1.0 - y_cheb**2.0) # 1 at centreline
            dU_dy  = np.identity(modes)
            np.fill_diagonal(dU_dy, -2.0*y_cheb)
            d2U_dy2 = -2.0
        elif bf == 'cou': # Couette Base flow
            U = np.identity(modes)
            np.fill_diagonal(U, y_cheb)
            dU_dy  = np.identity(modes)
            np.fill_diagonal(dU_dy, 1.0)
            d2U_dy2 = 0.0
    else: # Use turbulent mean
        U = np.identity(len(vel_profile))
        np.fill_diagonal(U, vel_profile)
        dU_dy = np.identity(len(vel_profile))
#        np.fill_diagonal(dU_dy, 0.0)
#        d2U_dy2 = 0.0
        np.fill_diagonal(dU_dy, -2.0*y_cheb) # laminar 
        d2U_dy2 = -2.0                       # laminar 


    #================================================================
    #### Calculate Laplacian for constructing operators
    #================================================================
    I = np.identity(modes)                              # identity matrix
    Z = np.zeros(shape=(modes, modes))                  # matrix of zeroes
    kappa = (alpha*alpha) + (beta*beta)                 # variable calculated from (del)dot(del)
    del_hat_2 = D2 - kappa*I                            # Laplacian
    del_hat_4 = D4 - 2.0*D2*kappa + kappa*kappa*I       # Laplacian_2

    #================================================================
    #### Construct operators
    #================================================================
    # Calculations from pg60 Schmid Henningson eqns 3.29 and 3.30
    # Squire operator
    O_sq = ((del_hat_2/Re) - (1.0j * alpha * U))
    # Center operator
    O_c = -1.0j*beta*dU_dy
    # Orr-Sommerfeld operator
    a0=(del_hat_4 / Re)
    a1=( 1.0j * alpha * d2U_dy2 * I)
    a2=(-1.0j * alpha * np.asmatrix(U) * np.asmatrix(del_hat_2))
    O_os = a0 + a1 + a2
    x0 = np.linalg.solve(del_hat_2, O_os)

    #================================================================
    #### Calculate state operator
    #================================================================
    # Equation 2.7
    # (Moarref, Model-based scaling of the streamwise energy density in 
    # high-Reynolds number turbulent channels, 2013)
    #
    # State operator
    # A = | x0   Z  |
    #     |  C   SQ |
    A = np.vstack((np.hstack((x0,Z)), np.hstack((O_c, O_sq))))

    #================================================================
    #### Calculate unweighted resolvent/transfer function
    #================================================================
    I = np.eye(A.shape[0]) # identity matrix (length of state operator)
    L = -1.0j*omega*I - A
    R = inv(L) # resolvent of A

    #================================================================
    #### Calculate conversion matrix
    #================================================================
    # C maps the resolvent from being a function of wall-normal velocity (v)
    # and vorticity (eta) to primitives (u, v, w).
    I = np.identity(modes) # identity matrix (length of modes)
    C = np.vstack((np.hstack(((1.0j/kappa) * (alpha*D1), (-1.0j/kappa) * (beta*I))), 
                   np.hstack((                        I,                   Z)), 
                   np.hstack(( (1.0j/kappa) * (beta*D1),  ( 1.0j/kappa) * (alpha*I)))))
    C = np.asmatrix(C)

    #================================================================
    #### Calculate Clenshaw–Curtis quadrature
    #================================================================
    tmp, clencurt_quadrature = clencurt(N)
    clencurt_quadrature = np.diag(np.sqrt(clencurt_quadrature[1:-1]))
    # Stack it 3 times for each velocity compoenent (u, v, w):
    clencurt_quadrature = np.vstack((np.hstack((clencurt_quadrature,Z,Z)),
                                     np.hstack((Z,clencurt_quadrature,Z)),
                                     np.hstack((Z,Z,clencurt_quadrature))))
    clencurt_quadrature = np.asmatrix(clencurt_quadrature)

    #================================================================
    #### Grid weighted conversion matrix
    #================================================================
    C_w = clencurt_quadrature * C
    # Adjoint of C maps the forcing vector to resolvent
#    w_inv_C_adj = C.H * np.linalg.inv(clencurt_quadrature)
    invw_C_adj = pinv(C_w) # this is the calculation used in Ati's code. Why?

    #================================================================
    #### Calculate weighted resolvent/transfer function
    #================================================================
    H = C_w * R * invw_C_adj

    #================================================================
    #### Store necassary variables
    #================================================================
    state_vecs = {}
    state_vecs['H'] = H
    state_vecs['cq'] = clencurt_quadrature
    state_vecs['y_cheb'] = y_cheb
    state_vecs['y_cheb_full'] = y_cheb_full
    state_vecs['D1'] = D1
    state_vecs['U'] = np.diag(U)
    
    return state_vecs


## Debugging
#vecs = get_state_vectors(0.5, 1.0, [], 4000.0, 75, np.pi/12.0, "lam")
#U, S, V = svd(vecs['H'])
#S = np.diag(S)