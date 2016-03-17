import numpy as np

import PseudoSpectralMethods as ps
import Tests
import Utils as ut

from scipy.sparse.linalg import svds
from numpy.linalg import svd

import sys

def resolvent_formulation(ffg):
    
    x = []
    y_cheb = np.zeros((ffg.Ny))
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

                    state_vecs = get_state_vectors(ffg, alpha, beta, x)
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

    return generated_ff, y_cheb



def resolvent_approximation(ffcf, rank, turb_mean_profile, ffmean):

    #================================================================
    # Create arrays of Fourier modes to use
    # (modes multiplied with fundamental wavenumbers)
    #================================================================
    kx = ffcf.Mx * ffcf.alpha
    kz = ffcf.Mz * ffcf.beta

    #================================================================
    # Add 2 to the amount of wall-normal points to be able to use 
    # chebyshev nodes. Since the wall-points are omitted in the calculation
    # of the transfer function, H.
    #
    # This is done so that all of the wall-normal nodes are used. (including at teh wall)
    #
    # Ny = num of modes when approximating (2 more than original)
    #================================================================
    Ny_orig = ffcf.Ny
    ffcf.set_Ny(Ny_orig+2)

    #================================================================
    # Ensure valid rank is specified
    #================================================================
    rank = min(rank, 3*ffcf.modes)
    ffcf.set_rank(rank)

    #================================================================
    # Initialize empty 4D array to store approximated velocity field
    #================================================================
    u_hat_approx  = np.zeros((len(kx), 3*Ny_orig, len(kz)), dtype=np.complex128)

    #================================================================
    # Loop through wavenumbers 
    #================================================================
    for mx in range(0, len(ffcf.Mx)):
        alpha = kx[mx]

        for mz in range(0, len(ffcf.Mz)):
            beta  = kz[mz]

            if alpha == 0 or beta == 0:
                #------------------------------------------------
                # Set the zero wavenumbers to equal the mean profile
                #------------------------------------------------
                u_hat_approx[mx, :, mz] = ffmean.velocityField[mx, :, mz]
                
                # Start the loop again
                continue

            text02 ='(mx)kx: ('+str(ffcf.Mx[mx])+') '+ str(alpha) + '\t'
            text02+='(mz)kz: ('+str(ffcf.Mz[mz])+') '+ str(beta)
            print(text02)

            #------------------------------------------------
            # Calculate the resolvent (R_A) and 
            # transfer function (H)
            #------------------------------------------------
            state_vecs = get_state_vectors(ffcf, alpha, beta, turb_mean_profile)

            #------------------------------------------------
            # Perform SVD (Singular Value Decomposition)
            #------------------------------------------------
            vel_modes, singular_values, forcing_modes = svds(state_vecs['H'], rank)
            vel_modes = np.asmatrix(vel_modes)
            #------------------------------------------------
            # Check it all went smoothly...
            #------------------------------------------------
#            Tests.SVDNorm(vel_modes, singular_values, forcing_modes, state_vecs['H'])

            #------------------------------------------------
            # Retrieve physical resolvent modes (non-grid-weighted)
            #------------------------------------------------
            resolvent_modes = np.linalg.solve(state_vecs['cq'], vel_modes)

            #------------------------------------------------
            # Check that the divergence criteria is satisfied
            #------------------------------------------------
#            Tests.divergence(resolvent_modes, alpha, beta, ffcf.modes, state_vecs['D1'])

            #------------------------------------------------
            # Check that the weighted resovlent 
            # and forcing modes are orthogonal.
            #------------------------------------------------
#            Tests.orthogonality(vel_modes)
#            Tests.orthogonality(forcing_modes.T)

            #------------------------------------------------
            # Fix phase of resolvent modes based on critical layer
            #------------------------------------------------
            phase_shift = np.zeros((resolvent_modes.shape[1], resolvent_modes.shape[1]), dtype=np.complex128)
            if ffcf.c < 1.0:
                inds = Tests.indices(state_vecs['U'], lambda x: x > ffcf.c)
                if len(inds) > 0:
                    ind0 = inds[0]
                    # ind1 = inds[-1]
                    # Don't need this point, because using the 
                    # first point multiplies both sides of the channel in 
                    # the wall-normal direction.
                elif len(inds) == 0: # Use centreline
                    ind0 = ffcf.Ny / 2 + 1
            else:
                ind0 = np.floor(ffcf.modes/2)

            phase_shift_tmp = np.exp(-1j * np.angle(resolvent_modes[ind0,:]))
            np.fill_diagonal(phase_shift, phase_shift_tmp)
            resolvent_modes *= phase_shift

            #------------------------------------------------
            # Project resolvent modes to get amplitude 
            # coefficients, chi, defined as
            # chi  = singular_values * eta
            #------------------------------------------------
            chi = get_scalars(ffcf.velocityField[mx, :, mz], resolvent_modes, state_vecs['cq'], rank)
#            chi = np.asarray(chi)
            chi = np.asmatrix(np.asarray(chi))
       
            #------------------------------------------------
            # Calculate approximate flow field.
            #------------------------------------------------
            u_hat_approx[mx, :, mz] = np.squeeze(np.asarray(vel_modes[:,:rank] * chi))


    #================================================================
    # Print difference in full and approximated velocity fields
    #================================================================
    diff = np.abs(ffcf.velocityField - u_hat_approx)
    diff = np.linalg.norm(diff)
    text03 = '\n|| u_original - u_approx || = '+str(diff)+'\n'
    print(text03)

    diff_orig = np.linalg.norm(ffcf.velocityField)
    diff_aprx = np.linalg.norm(u_hat_approx)
    diff = diff_orig - diff_aprx
    text02='\n|| u_original || - || u_approx || = '+str(diff)+'\n'
    print(text02)


    #================================================================
    # Reset Ny
    #================================================================
    # Take away 2 from the value of Ny set at the start of this method.    
    ffcf.set_Ny(Ny_orig) 


    #================================================================
    # Rearranging using numpy reshape
    #================================================================
    # u_hat_approx = u_hat_approx.reshape((ffcf.Nd, len(ffcf.Mx), ffcf.Ny, len(ffcf.Mz)))

    #================================================================
    # My way of re-arranging
    #================================================================
    U_hat = np.zeros((ffcf.Nd, len(ffcf.Mx), ffcf.Ny, len(ffcf.Mz)), dtype=np.complex128)
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

    ffcf.set_ff(U_hat, "sp")

    return ffcf




def resolvent_approximation2(ffcf, rank, turb_mean_profile, ffmean, sparse):
    print("\nUsing the new approximation method...")


    #================================================================
    #### Store the original flowfields in physical form for use at the end.
    #================================================================
    original_U = ffcf.velocityField[:, :, :, :]
    original_u = original_U[0, :, :, :].real
    original_v = original_U[1, :, :, :].real
    original_w = original_U[2, :, :, :].real

    # FFT test
    ffcf.make_xz_spectral()
    ffcf.make_xz_physical()

    diff = np.linalg.norm(original_U.real - ffcf.velocityField.real)
    if diff >= 1e-12:
        ut.error("FFT and IFFT went wrong...")

    #================================================================
    #### Remove the wall boundaries
    #================================================================
    print("\nRemoving the xz-plane at y=1 and y=-1,\nso that the chebyshev nodes can be used to construct the transfer function.")
    modified_U = np.zeros((ffcf.Nd, ffcf.Nx, ffcf.modes, ffcf.Nz))
    modified_U = ffcf.velocityField[:, :, 1:-1, :]
    ffcf.set_ff(modified_U.real, "pp")
    diff = np.linalg.norm(modified_U - ffcf.velocityField)
    if diff >= 1e-12:
        ut.error("The modified flow field has not been set correctly in the flow field class")

    modified_mean = np.zeros((ffcf.Nd, ffcf.Nx, ffcf.modes, ffcf.Nz))
    modified_mean = ffmean.velocityField[:, :, 1:-1, :]
    ffmean.set_ff(modified_mean.real, "pp")
    diff = np.linalg.norm(modified_mean - ffmean.velocityField)
    if diff >= 1e-12:
        ut.error("The modified mean flow field has not been set correctly in the flow field class")

    #================================================================
    #### Fourier transform the flow fields in xz directions
    #================================================================
    print("\nFFT instantaneous and turbu_mean flow fields")
    ffcf.make_xz_spectral()
    ffmean.make_xz_spectral()

    #================================================================
    #### Concatenate the flowfields so that they are stacked uvw in the wall-normal direction
    #================================================================
    spectral_U = np.concatenate((ffcf.velocityField[0, :, :, :],
                                 ffcf.velocityField[1, :, :, :],
                                 ffcf.velocityField[2, :, :, :]),
                                 axis=1)

    spectral_mean = np.concatenate((ffmean.velocityField[0, :, :, :],
                                    ffmean.velocityField[1, :, :, :],
                                    ffmean.velocityField[2, :, :, :]),
                                    axis=1)

    #================================================================
    #### Create arrays of Fourier modes to use
    # (modes multiplied with fundamental wavenumbers)
    #================================================================
    kx = ffcf.Mx * ffcf.alpha
    kz = ffcf.Mz * ffcf.beta

    #================================================================
    #### Ensure valid rank is specified
    #================================================================
    rank = min(rank, 3*ffcf.modes)
    ffcf.set_rank(rank)

    #================================================================
    #### Initialize empty 3D array to store approximated velocity field
    #================================================================
    u_hat_approx  = np.zeros((len(kx), 3*ffcf.modes, len(kz)), dtype=np.complex128)

    #================================================================
    #### Store the amplitude coefficients per streamwise Fourier mode
    #================================================================
    alpha_beta_chi = np.zeros((len(kx), len(kz), rank))

    #================================================================
    #### Loop through wavenumbers 
    #================================================================
    for mx in range(0, len(ffcf.Mx)):
        alpha = kx[mx]
        text02 ='\n(mx)kx: ('+str(ffcf.Mx[mx])+') '+ str(alpha)
        print(text02)

        for mz in range(0, len(ffcf.Mz)):
            beta  = kz[mz]

            sys.stdout.write(" "+ str(beta))
            sys.stdout.flush()

            if alpha == 0 or beta == 0:
                #------------------------------------------------
                #### Set the zero Fourier modes to equal the mean flow
                #------------------------------------------------
                if len(turb_mean_profile) == 0:
#                    print("No mean given")
#                    print(len(u_hat_approx[mx, :, mz]))
#                    print(len(spectral_U[mx, :, mz]))
                    u_hat_approx[mx, :, mz] = spectral_U[mx, :, mz]
                elif len(turb_mean_profile) != 0:
#                    print("Mean given")
#                    print(len(u_hat_approx[mx, :, mz]))
#                    print(len(spectral_U[mx, :, mz]))
                    u_hat_approx[mx, :, mz] = spectral_mean[mx, :, mz]

                # uvw of approximation at the zero Fourier modes equals the 
                # uvw of mean at the zero Fourier modes
                
                # Start the loop again
                continue

            #------------------------------------------------
            #### Calculate the resolvent (R_A) and transfer function (H)
            #------------------------------------------------
            state_vecs = get_state_vectors(ffcf, alpha, beta, turb_mean_profile[1:-1])
#            print("\n\nTransfer function shape:")
#            print(state_vecs['H'].shape)
#            print("")

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
            if not sparse:
                Tests.SVDNorm(vel_modes, singular_values, forcing_modes, state_vecs['H'])

            #------------------------------------------------
            #### Retrieve non-grid-weighted resolvent modes (physical modes)
            #------------------------------------------------
            resolvent_modes = np.linalg.solve(state_vecs['cq'], vel_modes)

            #------------------------------------------------
            #### Check that the divergence criteria is satisfied
            #------------------------------------------------
            if not sparse:
                Tests.divergence(resolvent_modes, alpha, beta, ffcf.modes, state_vecs['D1'])

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
            if ffcf.c < 1.0:
                inds = Tests.indices(state_vecs['U'], lambda x: x > ffcf.c)
                if len(inds) > 0:
                    ind0 = inds[0]
                elif len(inds) == 0: # Use centreline
                    ind0 = ffcf.Ny / 2 + 1
            else:
                ind0 = np.floor(ffcf.modes/2)

            phase_shift_tmp = np.exp(-1j * np.angle(resolvent_modes[ind0,:]))
            np.fill_diagonal(phase_shift, phase_shift_tmp)
            resolvent_modes *= phase_shift

            #------------------------------------------------
            #### Project resolvent modes to get amplitude coefficients
            # denoted chi, defined as
            # chi  = singular_values * eta
            #------------------------------------------------
            # Initialize the scalars vector
            chi = np.zeros((rank, 1), dtype=np.complex128)      
            # Projection
            chi = resolvent_modes[: , :rank].H * state_vecs['cq'].H * state_vecs['cq'] * np.asmatrix(spectral_U[mx, :, mz]).T 
            # Store the absolute value of the coefficients
            alpha_beta_chi[mx, mz, :] = np.squeeze(np.asarray(np.linalg.norm(chi)))

            #------------------------------------------------
            #### Calculate approximate flow field
            #------------------------------------------------
            # w * u_hat = w * psi * chi
            #     u_hat = [inv(w) * w ]* psi * chi
            #     u_hat = psi * chi
            tmp =  np.linalg.inv(state_vecs['cq']) * state_vecs['cq'] * resolvent_modes[:,:rank] * chi
            orig = np.asmatrix(spectral_U[mx, :, mz]).T
            diff = np.linalg.norm(tmp - orig)
            u_hat_approx[mx, :, mz] = np.squeeze(np.asarray(tmp))



    print("\nFinished approximating.\n")


    #================================================================
    # Print difference in full and approximated velocity fields
    #================================================================
#    diff = spectral_U - u_hat_approx
    diff = np.linalg.norm(spectral_U - u_hat_approx)
    text03 = '\n|| u_original - u_approx || = '+str(diff)
    print(text03)

    diff_orig = np.linalg.norm(spectral_U)
    diff_aprx = np.linalg.norm(u_hat_approx)
    diff = diff_orig - diff_aprx
    text02='\n|| u_original || - || u_approx || = '+str(diff)+'\n'
    print(text02)


    #================================================================
    #### Inverse Fourier transform the approximated flow field in xz directions
    #================================================================
    u_hat_approx = u_hat_approx.reshape((ffcf.Nd, len(ffcf.Mx), ffcf.modes, len(ffcf.Mz)))    
    print("Inverse FFT of approximated flow field")
    u_hat_approx_pp = np.fft.ifft(u_hat_approx, axis=1)
    u_hat_approx_pp = np.fft.ifft(u_hat_approx_pp, axis=3)


    u_hat_approx_pp_u = u_hat_approx_pp[0, :, :, :].real
    u_hat_approx_pp_v = u_hat_approx_pp[1, :, :, :].real
    u_hat_approx_pp_w = u_hat_approx_pp[2, :, :, :].real

    # IFFT the mean
    ffmean.make_xz_physical()
    
    #================================================================
    #### Add wall boundaries
    #================================================================
    approx_u = np.concatenate((original_U[0,:,:1,:].real,
                               u_hat_approx_pp_u.real,
                               original_U[0,:,-1:,:].real),
                               axis=1)
    diff = np.linalg.norm(approx_u.real - original_u.real)


    approx_v = np.concatenate((original_U[1,:,:1,:].real,
                               u_hat_approx_pp_v.real,
                               original_U[1,:,-1:,:].real),
                               axis=1)
    diff = np.linalg.norm(approx_v.real - original_v.real)


    approx_w = np.concatenate((original_U[2,:,:1,:].real,
                               u_hat_approx_pp_w.real,
                               original_U[2,:,-1:,:].real),
                               axis=1)

    diff = np.linalg.norm(approx_w.real - original_w.real)


    full_approx_U = np.zeros((ffcf.Nd, ffcf.Nx, ffcf.Ny, ffcf.Nz), dtype=np.float128)
    full_approx_U[0, :, :, :] = approx_u.real
    full_approx_U[1, :, :, :] = approx_v.real
    full_approx_U[2, :, :, :] = approx_w.real


    return full_approx_U, alpha_beta_chi





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
#    print("vp shape:\t", len(vel_profile))
#    print("U shape:\t", U.shape)
#    print("Lap len:\t", Lap.shape)
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
    state_vecs['U'] = np.diag(U)
    
    return state_vecs

