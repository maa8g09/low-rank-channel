import numpy as np
import sys

import PseudoSpectralMethods as ps
import Tests

from scipy.sparse.linalg import svds
from numpy.linalg import svd


def resolvent_formulation(ffg):

    tmp = np.zeros((ffg.Nx, ffg.Nd*ffg.modes, ffg.Nz), dtype=np.complex128)

    # Loop through wavenumber triplets
    for i in range(0, len(ffg.kx)):
        # Fundamental wavenumbers from wavenumber triplets
        fund_alpha = ffg.kx[i]
        fund_beta  = ffg.kz[i]

        # The stationary wave modes being used to calculate spectral flow field
        modes_streamwise = fund_alpha * ffg.Mx
        modes_spanwise   = fund_beta * ffg.Mz

        text01='alpha:'+ str(fund_alpha)+ '\t\tbeta:'+ str(fund_beta)+ '\t\tamplitude:'+ str(ffg.chi_tilde[i])
        
        print(text01)
        print('kx = mx * alpha\tkz = mz * beta')

        # Loop through the stationary modes
        for ia in range(0, len(modes_streamwise)):
            for ib in range(0, len(modes_spanwise)):
                    # Wavenumbers
                    alpha = modes_streamwise[ia]
                    beta = modes_spanwise[ib]

                    if alpha == 0 or beta == 0:
                        continue

                    text02 ='(mx)kx: ('+str(ffg.Mx[ia])+') '+ str(alpha) + '\t'
                    text02+='(mz)kz: ('+str(ffg.Mz[ib])+') '+ str(beta)
                    print(text02)

                    state_vecs = ps.get_state_vectors()
                    state_vecs['y_cheb_full'] = np.squeeze(np.asarray(state_vecs['y_cheb_full']))
                    y_cheb = state_vecs['y_cheb_full']

                    vel_modes, singular_values, forcing_modes_h = np.linalg.svd(state_vecs['H'])
                    
                    Tests.SVDNorm(vel_modes, singular_values, forcing_modes_h, state_vecs['H'])
                    Tests.orthogonality(vel_modes)
                    Tests.orthogonality(forcing_modes_h)

                    # Non-weighted resolvent modes (physical velocity modes)
                    resolvent_modes = np.linalg.solve(state_vecs['cq'], vel_modes)
                    # Non-weighted forcing modes (physical forcing modes)
                    unweighted_forcing_modes = np.linalg.solve(state_vecs['cq'], forcing_modes_h.conjugate().T)
                    
                    # Fix phase of first non-zero point based on critical layer
                    phase_shift = np.zeros((resolvent_modes[0,:].shape[1], resolvent_modes[0,:].shape[1]), dtype=np.complex128)
                    if ffg.c < 1.0:
                        inds = Tests.indices(state_vecs['U'], lambda x: x > ffg.c)
                        if len(inds) > 0:
                            ind0 = inds[0]
                            # ind1 = inds[-1]
                            # Don't need this point, because using the 
                            # first point multiplies both sides of the channel in 
                            # the wall-normal direction.
                        elif len(inds) == 0: # Use centreline
                            ind0 = ffg.Ny / 2 + 1
                    else:
                        ind0 = np.floor(ffg.modes/2)

                    phase_shift_tmp = np.exp(-1j * np.angle(resolvent_modes[ind0,:]))
                    np.fill_diagonal(phase_shift, phase_shift_tmp)
                    resolvent_modes *= phase_shift

# Alternative:
#                     # Non-weighted forcing modes (physical forcing modes)
#                    unweighted_forcing_modes = np.linalg.solve(state_vecs['cq'], forcing_modes_h.conjugate().T)
#                    unweighted_forcing_modes *= phase_shift
#                    unweighted_H = np.linalg.inv(state_vecs['cq']) * state_vecs['H'] * state_vecs['cq']
#                    unweighted_vel_modes = unweighted_H * unweighted_forcing_modes * np.diag(1.0/singular_values)
#                    
#                    delta_r = unweighted_vel_modes.real - resolvent_modes.real
#                    delta_i = unweighted_vel_modes.imag - resolvent_modes.imag

                    # Test the divergence
                    Tests.divergence(resolvent_modes, alpha, beta, ffg.modes, state_vecs['D1'])
                    resolvent_modes = np.asmatrix(resolvent_modes)

                    # u_tilde = chi_tilde * Psi
                    u_tilde = resolvent_modes[:, 0] * ffg.chi_tilde[i] # Rank 1
                    u_tilde = np.asmatrix(u_tilde)

                    # Inverse fourier transform
                    physical_ff = np.zeros((ffg.Nx, ffg.Nd*ffg.modes, ffg.Nz), dtype=np.complex128)
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
#    y_cheb = np.asarray(y_cheb)
    y_cheb = np.squeeze(np.asarray(y_cheb))
    
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

    generated_ff = np.zeros((ffg.Nd, ffg.Nx, ffg.Ny, ffg.Nz), dtype=np.float128)

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

    return generated_ff, y_cheb




def resolvent_approximation(original_ff_spectral,
                            mean_ff_spectral,
                            kx_array,
                            kz_array,
                            Nm,
                            c,
                            Re,
                            baseflow,
                            rank,
                            mean_profile,
                            sparse):

    #================================================================
    #### Initialize empty 3D array to store approximated velocity field
    #================================================================
    # Approximated velocity field is stacked in the wall-normal direction
    approximated_ff_spectral  = np.zeros((len(kx_array), 3*Nm, len(kz_array)), dtype=np.complex128)


    #================================================================
    #### Store the resolvent modes and amplitude coefficients 
    #    at each  Fourier mode pair
    #================================================================
    alpha_beta_resolvent_modes = np.zeros((len(kx_array), len(kz_array), 3*Nm, rank), dtype=np.complex128)
    alpha_beta_chi = np.zeros((len(kx_array), len(kz_array), rank), dtype=np.complex128)


    #================================================================
    #### Loop through wavenumbers 
    #================================================================
    for mx in range(0, len(kx_array)):
        kx = kx_array[mx]
        # print
        print('\n\nkx:'+ str(kx))

        for mz in range(0, len(kz_array)):
            kz  = kz_array[mz]
            # print kz
            sys.stdout.write(" "+ str(kz))
            sys.stdout.flush()

            if kx == 0 or kz == 0: # Zero Fourier modes
                #------------------------------------------------
                #### Set the zero Fourier modes to equal the mean flow
                #------------------------------------------------
                approximated_ff_spectral[mx, :, mz] = mean_ff_spectral[mx, :, mz]

                # [uvw] of approximation at the zero Fourier modes 
                # equals the [uvw] of mean at the zero Fourier modes.
                
                # Therefore the projection will be of the total flow,
                # i.e. original + mean flow

                continue # Start the loop again

            #------------------------------------------------
            #### Calculate the state vectors
            #------------------------------------------------
            omega = kx * c
            state_vecs = ps.get_state_vectors(kx, kz, Re, Nm+2, omega, baseflow, mean_profile)


            #------------------------------------------------
            #### Perform SVD
            #------------------------------------------------
            if sparse:
                #vel_modes, singular_values, forcing_modes = svds(state_vecs['H'], rank)
                vel_modes, singular_values, forcing_modes = svd(state_vecs['H'], full_matrices=False)
            else:
                vel_modes, singular_values, forcing_modes = svd(state_vecs['H'])

            vel_modes = np.asmatrix(vel_modes)


            #------------------------------------------------
            #### Check SVD
            #------------------------------------------------
            Tests.SVD(vel_modes, singular_values, forcing_modes, state_vecs['H'], sparse)


            #------------------------------------------------
            #### Retrieve non-grid-weighted resolvent modes (physical modes)
            #------------------------------------------------
            resolvent_modes = np.linalg.solve(state_vecs['cq'], vel_modes)


            #------------------------------------------------
            #### Check that the continuity criteria is satisfied
            #------------------------------------------------
            if not sparse:
                Tests.continuity(resolvent_modes, kx, kz, Nm, state_vecs['D1'])
#            elif sparse:
#                Tests.continuity(resolvent_modes[: , :rank], kx, kz, Nm, state_vecs['D1'])

            #------------------------------------------------
            #### Check that the weighted resovlent and forcing modes are orthogonal.
            #------------------------------------------------
            if not sparse:
                Tests.orthogonality(vel_modes)
                Tests.orthogonality(forcing_modes.T)


            #------------------------------------------------
            #### Fix phase of resolvent modes based on critical layer or centreline
            #------------------------------------------------
            phase_shift = np.zeros((resolvent_modes.shape[1], resolvent_modes.shape[1]), dtype=np.complex128)
            ind0 = Nm/2 + 1 # Use centreline, unless
            if c < 1.0:
                inds = Tests.indices(state_vecs['U'], lambda x: x > c)
                if len(inds) > 0:
                    ind0 = inds[0]

            np.fill_diagonal(phase_shift, np.exp(-1j * np.angle(resolvent_modes[ind0,:])))
            resolvent_modes *= phase_shift


            #------------------------------------------------
            #### Project resolvent modes to get amplitude coefficients
            #------------------------------------------------
            # denoted chi, defined as
            # chi  = singular_values * eta

            # Initialize the scalars vector
            chi = np.zeros((rank, 1), dtype=np.complex128)    

            # Projection
            chi = resolvent_modes[: , :rank].H * state_vecs['cq'].H * state_vecs['cq'] * np.asmatrix(original_ff_spectral[mx, :, mz]).T

            # Store the resolvent modes and amplitudes 
            # for reconstruction at a later date
            alpha_beta_resolvent_modes[mx, mz, :, :] = resolvent_modes[: , :rank]
            alpha_beta_chi[mx, mz, :] = np.squeeze(np.asarray(chi))


            #------------------------------------------------
            #### Calculate approximate flow field
            #------------------------------------------------
            # w * u_hat = w * psi * chi
            #     u_hat = [inv(w) * w ]* psi * chi
            #     u_hat = psi * chi
            tmp =  np.linalg.inv(state_vecs['cq']) * state_vecs['cq'] * resolvent_modes[:,:rank] * chi
            orig = np.asmatrix(original_ff_spectral[mx, :, mz]).T
            diff = np.linalg.norm(tmp - orig)
            approximated_ff_spectral[mx, :, mz] = np.squeeze(np.asarray(tmp))


    projected_field = {}
    projected_field['approximated_ff_spectral'] = approximated_ff_spectral
    projected_field['alpha_beta_resolvent_modes'] = alpha_beta_resolvent_modes
    projected_field['alpha_beta_chi'] = alpha_beta_chi

    return approximated_ff_spectral, alpha_beta_chi



