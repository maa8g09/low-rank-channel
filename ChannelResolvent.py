import numpy as np
import PseudoSpectralMethods as ps
import Tests


def resolvent_formulation(ffg):
    x = []
    tmp = np.zeros((ffg.Nx, ffg.Nd*ffg.modes, ffg.Nz), dtype=np.complex128)

    # Loop through wavenumber triplets
    for i in range(0, len(ffg.kx)):
        # Fundamental wavenumbers from wavenumber triplets
        fund_alpha = ffg.kx[i]
        fund_beta  = ffg.kz[i]

        # The stationary wave modes being used to calculate spectral flow field
        modes_streamwise = fund_alpha * ffg.Mx
        modes_spanwise   = fund_beta * ffg.Mz

        text01='alpha:'+ str(fund_alpha)+ '  beta:'+ str(fund_beta)+ '    amplitude:'+ str(ffg.chi_tilde[i])
        print(text01)
        print('kx = mx * alpha        kz = mz * beta')

        # Loop through the stationary modes
        for ia in range(0, len(modes_streamwise)):
            for ib in range(0, len(modes_spanwise)):
                    # Wavenumbers
                    alpha = modes_streamwise[ia]
                    beta = modes_spanwise[ib]

                    if alpha == 0 or beta == 0:
                        continue

                    text02='(mx)kx: ('+str(ffg.Mx[ia])+') '+ str(alpha)+'    (mz)kz: ('+str(ffg.Mz[ib])+') '+ str(beta)
                    print(text02)

                    state_vecs = get_state_vectors(ffg, alpha, beta, x)
                    state_vecs['y_cheb_full'] = np.squeeze(np.asarray(state_vecs['y_cheb_full']))
                    ffg.set_y(state_vecs['y_cheb_full'])

                    vel_modes, singular_values, forcing_modes = np.linalg.svd(state_vecs['H'])
                    Tests.SVDNorm(vel_modes, singular_values, forcing_modes, state_vecs['H'])
                    Tests.orthogonality(vel_modes)
                    Tests.orthogonality(forcing_modes)

                    # Non-weighted resolvent modes (physical velocity modes)
                    resolvent_modes = np.linalg.solve(state_vecs['cq'], vel_modes)
                    Tests.divergence(resolvent_modes, alpha, beta, ffg.modes, state_vecs['D1'])
                    resolvent_modes = np.asmatrix(resolvent_modes)

                    # u_tilde = chi_tilde * Psi
                    u_tilde = resolvent_modes[:, 0] * ffg.chi_tilde[i] # Rank 1
                    u_tilde = np.asmatrix(u_tilde)

                    # Inverse fourier transform
                    physical_ff = np.zeros((ffg.Nx, 3*ffg.modes, ffg.Nz), dtype=np.complex128)
                    physical_ff = ps.my_ifft(u_tilde[:,0], alpha, beta, ffg)
                    tmp += physical_ff

        print('\n\n')

#    # Rearranging using numpy reshape
#    tmp2 = tmp.reshape((ffg.Nd, ffg.Nx, ffg.modes, ffg.Nz))
#    generated_ff2 = np.zeros((ffg.Nd, ffg.Nx, ffg.Ny, ffg.Nz), dtype=np.float64)
#    
#    if ffg.baseflow == 'lam':
#        generated_ff2[:, :,         1:ffg.Ny-1, :]   = tmp2[:, :,          0:ffg.modes, :].real
#        generated_ff2[:, :,  ffg.Ny+1:ffg.Ny*2-1, :] = tmp2[:, :,  ffg.modes:ffg.modes*2, :].real
#        generated_ff2[:, :,2*ffg.Ny+1:ffg.Ny*3-1, :] = tmp2[:, :,2*ffg.modes:ffg.modes*3, :].real
#        
#    elif ffg.baseflow == 'cou':
#        generated_ff_pos = np.ones((ffg.Nd, ffg.Nx, 0.5*ffg.Ny, ffg.Nz))
#        generated_ff_neg = -1.0*generated_ff_pos
#        
#        generated_ff = np.concatenate(((generated_ff_pos), (generated_ff_neg)), axis=2)
#        
#        generated_ff[:, :,         1:ffg.Ny-1, :]   = tmp2[:, :,          0:ffg.modes, :].real
#        generated_ff[:, :,  ffg.Ny+1:ffg.Ny*2-1, :] = tmp2[:, :,  ffg.modes:ffg.modes*2, :].real
#        generated_ff[:, :,2*ffg.Ny+1:ffg.Ny*3-1, :] = tmp2[:, :,2*ffg.modes:ffg.modes*3, :].real


    # My way of re-arranging...
    if ffg.baseflow == "lam":
        top_boundary = np.zeros((ffg.Nx, 1, ffg.Nz)) # no-slip boundary condition
        bot_boundary = top_boundary # no-slip boundary condition
    elif ffg.baseflow == "cou":
        top_boundary = np.ones((ffg.Nx, 1, ffg.Nz))
        bot_boundary = -1.0*top_boundary


    U_u = tmp.real[:,           0:ffg.modes  , :]
    U_v = tmp.real[:,   ffg.modes:ffg.modes*2, :]
    U_w = tmp.real[:, 2*ffg.modes:ffg.modes*3, :]

    generated_ff = np.zeros((ffg.Nd, ffg.Nx, ffg.Ny, ffg.Nz), dtype=np.float64)

    U_u = np.concatenate((top_boundary[:,:,:],
                          U_u[:,:,:],
                          bot_boundary[:,:,:]),
                          axis=1)

    U_v = np.concatenate((top_boundary[:,:,:],
                          U_v[:,:,:],
                          bot_boundary[:,:,:]),
                          axis=1)

    U_w = np.concatenate((top_boundary[:,:,:],
                          U_w[:,:,:],
                          bot_boundary[:,:,:]),
                          axis=1)

    for i in range(0, ffg.Nd):
        for nx in range(0, ffg.Nx):
            for ny in range(0, ffg.Ny):
                for nz in range(0, ffg.Nz):
                    if i == 0: # u direction
                        generated_ff[i, nx, ny, nz] = U_u[nx, ny, nz]
                    elif i == 1: # v direction
                        generated_ff[i, nx, ny, nz] = U_v[nx, ny, nz]
                    elif i == 2: # w direction
                        generated_ff[i, nx, ny, nz] = U_w[nx, ny, nz]

    return generated_ff



def resolvent_approximation(ffcf, rank, turb_profile, ffmean):

    kx = ffcf.Mx * ffcf.alpha      # list of wavenumbers to use (modes multiplied by fundamental alpha)
    kz = ffcf.Mz * ffcf.beta       # list of wavenumbers to use (modes multiplied by fundamental beta)
    
    # Ny = num of modes when approximating    
    Ny_orig = ffcf.Ny
    ffcf.set_Ny(Ny_orig+2)

    rank = min(rank, 3*ffcf.modes)
    ffcf.set_rank(rank)

    u_hat_approx  = np.zeros((len(kx), 3*Ny_orig, len(kz)), dtype=np.complex128)

    for mx in range(0, len(ffcf.Mx)):
        alpha = kx[mx]


        for mz in range(0, len(ffcf.Mz)):
            beta  = kz[mz]

            if alpha == 0 or beta == 0:
                # Make sure that the mean is copied over to the approximation
                u_hat_approx[mx, :, mz] = ffmean.velocityField[mx, :, mz]
                continue

            state_vecs = get_state_vectors(ffcf, alpha, beta, turb_profile)
            state_vecs['y_cheb_full'] = np.squeeze(np.asarray(state_vecs['y_cheb_full']))
            
            vel_modes, singular_values, forcing_modes = np.linalg.svd(state_vecs['H'])
            Tests.SVDNorm(vel_modes, singular_values, forcing_modes, state_vecs['H'])
            
            # Non-weighted resolvent modes
            resolvent_modes = np.linalg.solve(state_vecs['cq'], vel_modes)
            Tests.divergence(resolvent_modes, alpha, beta, ffcf.modes, state_vecs['D1'])
            Tests.orthogonality(vel_modes)

            # chi  = singular_values * eta
            chi = get_scalars(ffcf.ff[mx, :, mz], resolvent_modes, state_vecs['cq'], rank)
            chi = np.asarray(chi)
            chi = np.asmatrix(chi)
            print('chi shape')
            print(chi.shape)
            
            u_tilde_approx = vel_modes[:,:rank] * chi

            result = np.asmatrix(ffcf.ff[mx, :, mz]).T - u_tilde_approx
            result = np.linalg.norm(result)

#            text03='|| u_original - u_approx || = ' + str(result)
#            print(text03)

            u_hat_approx[mx, :, mz] = np.squeeze(np.asarray(u_tilde_approx))


    diff = np.abs(ffcf.ff - u_hat_approx)
    diff = np.linalg.norm(diff)
#    text01 = '\n|| u_original || =' + format(np.linalg.norm(ffcf.ff), '.8f')
#    text02 = '\n || u_approx ||  =' + format(np.linalg.norm(u_hat_approx), '.8f')
    text03 = '\n|| u_original - u_approx || = '+str(diff)+'\n'
#    print(text01)
#    print(text02)
    print(text03)
    diff_orig = np.linalg.norm(ffcf.ff)
    diff_aprx = np.linalg.norm(u_hat_approx)
    diff = diff_orig - diff_aprx
    text02='\n|| u_original || - || u_approx || = '+str(diff)+'\n'
    print(text02)

    ffcf.set_Ny(Ny_orig) # take away 2 from the value of Ny set at the start of this method.
    # this was done so that all of the wall-normal units are used.
    
#    # Rearranging using numpy reshape
#    u_hat_approx = u_hat_approx.reshape((ffcf.Nd, len(ffcf.Mx), ffcf.Ny, len(ffcf.Mz)))

    # My way of re-arranging...
    U_hat = np.zeros((ffcf.Nd, ffcf.Nx, ffcf.Ny, ffcf.Nz), dtype=np.complex128)
    U_hat_u = u_hat_approx[:,         0:ffcf.Ny  , :]
    U_hat_v = u_hat_approx[:,   ffcf.Ny:2*ffcf.Ny, :]
    U_hat_w = u_hat_approx[:, 2*ffcf.Ny:3*ffcf.Ny, :]

    for i in range(0, ffcf.Nd):
        for mx in range(0, len(ffcf.Mx)):
            for ny in range(0, ffcf.Ny):
                for mz in range(0, len(ffcf.Mz)):
                    if i == 0: # u direction
                        U_hat[i, mx, ny, mz] = U_hat_u[mx, ny, mz]
                    elif i == 1: # v direction
                        U_hat[i, mx, ny, mz] = U_hat_v[mx, ny, mz]
                    elif i == 2: # w direction
                        U_hat[i, mx, ny, mz] = U_hat_w[mx, ny, mz]

    ffcf.set_ff(U_hat)

    return ffcf



def get_state_vectors(ffg, alpha, beta, vel_profile):
    """
    We are calculating the state vectors in this function. The methodology
    followed here is given in the following reference in the "Formulation" 
    section, 
        2. Low-rank approximation to channel flow,
        
        (Moarref, Model-based scaling of the streamwise energy 
        density in high-Reynolds number turbulent channels, 2013)
    
    
    INPUTS:
            alpha:  streamwise wavenumber already in 2pialpha/Lx state
            beta:  spanwise wavenumber already in 2pibeta/Lz state
            
     
    OUTPUTS:
             C:  this operator maps the state vector onto the velocity vector
         C_adj:  adjoint of C, maps the forcing vector to the state vector
             A:  state operator
             W:  ClenCurt matrix (Clenshaw-Curtis quadrature)
             y:  grid-points in the y-axis
    """
    
    ## use the fundamental wavenumbers to multiply by the iterative tuples you use
    # i.e. if .geom file says:alpha = 1.14
    # then each alpha is a multiple of alpha
    
    
    
    # Calculate the differentiation matrix, DM, for the resolution in the 
    # y-axis, N. 
    # y_cheb are the interpolated y co-ordinates, i.e. Chebyshev interior points.
    y_cheb_full, DM = ps.chebdiff(ffg.Ny, 2)
    
    # First derivative matrix
    D1 = DM[0, 1:-1, 1:-1]
    
    # Second derivative matrix
    D2 = DM[1, 1:-1, 1:-1]


    # Fourth derivative matrix and clamped boundary conditions
    y_cheb, D4 = ps.cheb4c(ffg.Ny, False)
    # tmp is the same as y_cheb without endpoints, i.e. no [-1,1]
    
    
    # For the Orr-Sommerfeld equations we need to calculate the derivates
    #            D:  partial_dy
    #          Lap:  D**2.0 - K**2.0         (where K**2.0 = alpha**2.0 + beta**2.0)
    #         dUdy:  first derivative of Uo(y)
    #        dU2dy:  second derivative of Uo(y)
    #            f:  time derivative
    
    
    # Number of modes
    I = np.identity(ffg.modes)
    Z = np.zeros(shape=(ffg.modes, ffg.modes))
    K2 = (alpha**2.0) + (beta**2.0)
    Lap = D2 - K2*I #Laplacian
    del_hat_4 = D4 - 2.0*D2*K2 + K2*K2*I
    
    #print("\nVelocity profile length: " + str(len(vel_profile)))
    #print("")
    #print(vel_profile)
    if len(vel_profile) == 0:
        if ffg.baseflow == 'lam':
            # Laminar Base flow 
            U = np.identity(ffg.modes)
            np.fill_diagonal(U, 1.0 - y_cheb**2.0) # 1 at centreline
            
            dU_dy  = np.identity(ffg.modes)
            np.fill_diagonal(dU_dy, -2.0*y_cheb)
            
            d2U_dy2 = -2.0
    
        elif ffg.baseflow == 'cou':
            # Couette Base flow
            U = np.identity(ffg.modes)
            np.fill_diagonal(U, y_cheb)
            
            dU_dy  = np.identity(ffg.modes)
            np.fill_diagonal(dU_dy, 1.0)
            
            d2U_dy2 = 0.0
        
    else:
        # Use Turbulent Mean
        U = np.identity(len(vel_profile))
        np.fill_diagonal(U, vel_profile)
        
        dU_dy = np.identity(len(vel_profile))
        np.fill_diagonal(dU_dy, 0.0)
    
        d2U_dy2 = 0.0

    # pg 60 Schmid Henningson eqns 3.29 and 3.30
    SQ_operator = ((Lap/ffg.Re) - (1.0j * alpha * U))
    C_operator  = -1.0j*beta*dU_dy
    
    a0=(del_hat_4 / ffg.Re)
    a1=( 1.0j * alpha * d2U_dy2 * I)
    a2=(-1.0j * alpha * np.asmatrix(U) * np.asmatrix(Lap))
    
    OS_operator = a0 + a1 + a2
    x0 = np.linalg.solve(Lap, OS_operator)
    
    # Equation 2.7
    # (Moarref, Model-based scaling of the streamwise energy density in 
    # high-Reynolds number turbulent channels, 2013)
    #
    # State operator
    # A = | x0   Z  |
    #     |  C   SQ |
    A = np.vstack((np.hstack((x0,Z)), np.hstack((C_operator, SQ_operator))))

 
    # C maps state vector to the velocity vector
    C = np.vstack((np.hstack(((1.0j/K2) * (alpha*D1), (-1.0j/K2) * (beta*I))), 
                   np.hstack((                  I,                   Z)), 
                   np.hstack(((1.0j/K2) * (beta*D1), ( 1.0j/K2) * (alpha*I)))))
    
    C = np.asmatrix(C)
    
    
    tmp, clencurt_quadrature = ps.clencurt(ffg.Ny)
    clencurt_quadrature = np.diag(np.sqrt(clencurt_quadrature[1:-1]))
    clencurt_quadrature = np.vstack((np.hstack((clencurt_quadrature,Z,Z)),
                                     np.hstack((Z,clencurt_quadrature,Z)),
                                     np.hstack((Z,Z,clencurt_quadrature))))

    clencurt_quadrature = np.asmatrix(clencurt_quadrature)
    
    C_w = clencurt_quadrature * C
    
    # Adjoint of C maps the forcing vector to the state vector
    w_inv_C_adj = C.conjugate().T * np.linalg.inv(clencurt_quadrature)

    # Calculating the resolvent operator
    I = np.eye(A.shape[0])
    omega = alpha * ffg.c
    L = 1.0j*omega*I + A
    Linv = np.linalg.inv(L)
    R_A = -1.0 * Linv # resolvent of A

    # Calculating the transfer function
    H = C_w * R_A * w_inv_C_adj # transfer function
    
    state_vecs = {}
    state_vecs['H'] = H
    state_vecs['cq'] = clencurt_quadrature
    state_vecs['y_cheb'] = y_cheb
    state_vecs['y_cheb_full'] = y_cheb_full
    state_vecs['D1'] = D1
    
    return state_vecs



def get_scalars(u_hat, resolvent_modes, cq, rank):










    #========================================================================
    # Projecting with the required amount of column vectors==================
    #========================================================================
    vel_modes = np.asmatrix(vel_modes) # vel_modes have already been mutltiplied by the clencurt quadrature
    psi = vel_modes[: , :rank] # column vectors

    # Get the complex conjugate of the modes.
    psi_star = psi.H

    # Initialize the scalars vector (shape = Nd*Ny, long vector for u, v, w)
    chi = np.zeros((rank, 1), dtype=np.complex128)

    # Convert from array to matrix for multiplication later
    u_hat = np.asmatrix(u_hat).T

    chi = psi_star * u_hat

    return chi
    
