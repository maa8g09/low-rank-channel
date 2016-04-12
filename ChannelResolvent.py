import numpy as np
import sys
import time
import PseudoSpectralMethods as ps
import Tests
from datetime import datetime
from scipy.sparse.linalg import svds
from numpy.linalg import svd
from numpy.linalg import inv


def resolvent_formulation(ffg):
    #================================================================
    #### Initialize empty 3D array to store approximated velocity field
    #================================================================
    # Velocity field is stacked in the wall-normal direction.
    tmp = np.zeros((ffg.Nx, ffg.Nd*ffg.modes, ffg.Nz), dtype=np.complex128)
    # Note that the number of gridpoints in the wall-normal direction
    # are equal to the number of interior Chebyshev nodes. This is 
    # becasue when reshaped, the velocity field is missing its 
    # wall-boundaries, i.e. top and bottom xz-planes.

    #================================================================
    #### Loop through wavenumber triplets
    #================================================================
    for triplet in range(0, len(ffg.kx)):
        #============================================================
        #### Get fundamental wavenumbers from wavenumber triplets
        #============================================================
        fund_alpha = ffg.kx[triplet]
        fund_beta  = ffg.kz[triplet]

        #============================================================
        #### Calculate Fourier modes
        #============================================================
        kx_array = fund_alpha * ffg.Mx
        kz_array = fund_beta * ffg.Mz

        #============================================================
        #### Loop through the Fourier modes
        #============================================================
        for mx in range(0, len(kx_array)):
            for mz in range(0, len(kz_array)):
                    #------------------------------------------------
                    # Wavenumbers
                    #------------------------------------------------
                    kx = kx_array[mx] # Streamwise Fourier mode
                    kz = kz_array[mz] # Spanwise Fourier mode

                    if kx == 0 or kz == 0:
                        continue

                    #------------------------------------------------
                    #### Calculate the state vectors
                    #------------------------------------------------
                    omega = kx * ffg.c
                    state_vecs = ps.get_state_vectors(kx, kz, ffg.Re, ffg.Ny, omega, ffg.baseflow, [])
                    state_vecs['y_cheb_full'] = np.squeeze(np.asarray(state_vecs['y_cheb_full']))
                    y_cheb = state_vecs['y_cheb_full']

                    vel_modes, singular_values, forcing_modes_h = np.linalg.svd(state_vecs['H'])


                    Tests.SVDNorm(vel_modes, singular_values, forcing_modes_h, state_vecs['H'])
                    Tests.orthogonality(vel_modes)
                    Tests.orthogonality(forcing_modes_h)


                    #------------------------------------------------
                    # Non-weighted resolvent modes (physical velocity modes)
                    #------------------------------------------------
                    resolvent_modes = np.linalg.solve(state_vecs['w'], vel_modes)


                    #------------------------------------------------
                    # Non-weighted forcing modes (physical forcing modes)
                    #------------------------------------------------
                    unweighted_forcing_modes = np.linalg.solve(state_vecs['w'], forcing_modes_h.conjugate().T)


                    #------------------------------------------------
                    # Fix phase of first non-zero point based on critical layer
                    #------------------------------------------------
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
#                    unweighted_forcing_modes = np.linalg.solve(state_vecs['w'], forcing_modes_h.conjugate().T)
#                    unweighted_forcing_modes *= phase_shift
#                    unweighted_H = np.linalg.inv(state_vecs['w']) * state_vecs['H'] * state_vecs['w']
#                    unweighted_vel_modes = unweighted_H * unweighted_forcing_modes * np.diag(1.0/singular_values)
#                    
#                    delta_r = unweighted_vel_modes.real - resolvent_modes.real
#                    delta_i = unweighted_vel_modes.imag - resolvent_modes.imag


                    #------------------------------------------------
                    # Test the divergence
                    #------------------------------------------------
                    Tests.divergence(resolvent_modes, alpha, beta, ffg.modes, state_vecs['D1'])
                    resolvent_modes = np.asmatrix(resolvent_modes)


                    #------------------------------------------------
                    # u_tilde = chi_tilde * Psi
                    #------------------------------------------------
                    u_tilde = resolvent_modes[:, 0] * ffg.chi_tilde[i] # Rank 1
                    u_tilde = np.asmatrix(u_tilde)


                    #------------------------------------------------
                    # Inverse fourier transform
                    #------------------------------------------------
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




def deconstruct_field(original_ff_spectral,
                     kx_array,
                     kz_array,
                     Nm,
                     c,
                     Re,
                     baseflow,
                     r,
                     mean_profile,
                     sparse,
                     fixXi):

    ''' 
    Deconstruct given flow field and return the resolvent modes,
    singular values and coefficients needed to reconstruct it.

    ================================================================
    INPUTS:
    ================================================================
    original_ff_spectral:           3D array of flow field to 
                                    approximate in xz spectral form 
                                    with Chebyshev nodes in wall-normal 
                                    direction (Chebyshev spacing).
                                    Stacked in wall-normal direction:
                                    dimensions: (Nx, Nd*Ny, Nz).

    kx_array:                       1D array of streamwise Fourier modes.

    kz_array:                       1D array of spanwise Fourier modes.

    Nm:                             Number of interior Chebyshev nodes,
                                    i.e. Ny - 2 (endpoints removed).

    c:                              Wavespeed.

    Re:                             Reynolds number based on laminar
                                    baseflow centreline velocity, i.e.
                                    Re = 1 / nu.

    baseflow:                       Base flow type: [laminar, Couette],
                                    laminar means Plane Poiseuille.

    r:                              Rank to approximate to.

    mean_profile:                   1D array of streamwise velocity 
                                    profile of the turbulent mean in
                                    wall-normal direction, 
                                    (endpoints included).

    sparse:                         Boolean that determines whether to 
                                    use sparse or full SVD algorithm.

    ================================================================
    OUTPUTS:
    ================================================================
    A dictionary named 'deconstructed_dict' which contains:
    
    resolvent_modes:                4D array of truncated resolvent 
                                    modes at each Fourier mode 
                                    combination:
                                        resolvent_modes[mx, mz, :, :] gives 
                                        column vectors of given rank, i.e.
                                        rank = resolvent_modes.shape[3].

    singular_values:                3D array of truncated singular
                                    values at each Fourier mode 
                                    combination:
                                        singular_values[mx, mz, :] gives
                                        1D array of singular values, where
                                        rank = len(singular_values[mx, mz, :]).

    coefficients:                   3D array of truncated complex 
                                    amplitude coefficients at each 
                                    Fourier mode combination:
                                        coefficients[mx, mz, :] gives
                                        1D array of coefficients, where
                                        rank = len(coefficients[mx, mz, :]).
    '''


    #================================================================
    #### Store the resolvent modes and amplitude coefficients 
    #    at each  Fourier mode pair
    #================================================================
    resolvent_modes_array = np.zeros((len(kx_array), len(kz_array), 3*Nm, r), dtype=np.complex128)
    coefficients_array = np.zeros((len(kx_array), len(kz_array), r), dtype=np.complex128)
    sing_vals_array = np.zeros((len(kx_array), len(kz_array), r), dtype=np.float64)

    original_rank = r

    #================================================================
    #### Loop through wavenumbers 
    #================================================================
    chebyshev_differentiation, mean_flow_derivatives = ps.calculate_derivatives(Nm+2, mean_profile, baseflow)

    startTime = datetime.now()
    #================================================================
    #### Loop through wavenumbers 
    #================================================================
    for mx in range(0, len(kx_array)):
        kx = kx_array[mx]
        print('\n\nkx:'+ str(kx))

        for mz in range(0, len(kz_array)):
            kz  = kz_array[mz]
            sys.stdout.write(".")
            sys.stdout.flush()

            if kx == 0 or kz == 0: # Zero Fourier modes
                # Save zeroth modes to the resolvent modes...
                resolvent_modes_array[mx, mz, :, 0] = original_ff_spectral[mx,:,mz]
                continue # Start the loop again

            #------------------------------------------------
            #### Calculate the state vectors
            #------------------------------------------------
            omega = kx * c
            wegihted_transfer_function, w = ps.calculate_transfer_function(kx, kz, Re, Nm, omega, chebyshev_differentiation, mean_flow_derivatives)


            #------------------------------------------------
            #### Perform SVD
            #------------------------------------------------
            if sparse:
#                if r == min(wegihted_transfer_function.shape):
#                    vel_modes, singular_values, forcing_modes = svds(wegihted_transfer_function)
#                    
#                elif r < min(wegihted_transfer_function.shape):
#                    vel_modes, singular_values, forcing_modes = svds(wegihted_transfer_function, r)
#                r = len(singular_values)

                vel_modes, singular_values, forcing_modes = svd(wegihted_transfer_function, full_matrices=False)
            else:
                vel_modes, singular_values, forcing_modes = svd(wegihted_transfer_function)

            vel_modes = np.asmatrix(vel_modes)
            forcing_modes = np.asmatrix(forcing_modes)

            #------------------------------------------------
            #### Check SVD
            #------------------------------------------------
            Tests.SVD(vel_modes, singular_values, forcing_modes, wegihted_transfer_function, sparse)


            #------------------------------------------------
            #### Retrieve non-grid-weighted resolvent modes (physical modes)
            #------------------------------------------------
            resolvent_modes = np.linalg.solve(w, vel_modes)
            resolvent_modes = resolvent_modes[:,:r]

            #------------------------------------------------
            #### Check that the singular value matrix is invertible
            #------------------------------------------------
            S = np.diag(singular_values[:r])
            Tests.invertible(S)
            
            
            #------------------------------------------------
            #### Check that the continuity condition is satisfied
            #------------------------------------------------
            Tests.continuity(resolvent_modes, S, kx, kz, Nm, chebyshev_differentiation['D1'])
#            elif sparse:
#                Tests.continuity(resolvent_modes[: , :r], kx, kz, Nm, state_vecs['D1'])

            #------------------------------------------------
            #### Check that the weighted resovlent and forcing modes are orthogonal.
            #------------------------------------------------
            Tests.orthogonality(vel_modes)
            Tests.orthogonality(forcing_modes)

            Tests.orthogonality(vel_modes[:,:r])
            Tests.orthogonality(forcing_modes[:,:r])

            #------------------------------------------------
            #### Fix phase of resolvent modes based on critical layer or centreline
            #------------------------------------------------
            phase_shift = np.zeros((resolvent_modes.shape[1], resolvent_modes.shape[1]), dtype=np.complex128)
            ind0 = Nm/2 + 1 # Use centreline, unless
            if c < 1.0:
                inds = Tests.indices(mean_flow_derivatives['U'].diagonal(), lambda x: x > c)
                if len(inds) > 0:
                    ind0 = inds[0]

            np.fill_diagonal(phase_shift, np.exp(-1j * np.angle(resolvent_modes[ind0,:])))
            resolvent_modes *= phase_shift

#            # What do my modes look like?
#            resolvent_modes_1Ru = resolvent_modes[0:Nm,0].real
#            resolvent_modes_1Rv = resolvent_modes[Nm:Nm*2,0].real
#            resolvent_modes_1Rw = resolvent_modes[Nm*2:Nm*3,0].real


            #------------------------------------------------
            #### Project resolvent modes to get amplitude coefficients
            #------------------------------------------------
            # denoted xi, defined as
            # chi  = singular_values * xi

            # Initialize the scalars vector
            xi = np.zeros((r, 1), dtype=np.complex128)    

            # Projection
            xi[:r] = inv(S) * resolvent_modes.H * w.H * w * np.asmatrix(original_ff_spectral[mx, :, mz]).T
            test_Xi = resolvent_modes * S * xi[:r] # this variable should be equal to the original vector
            o = np.asmatrix(original_ff_spectral[mx, :, mz]).T          
            d = o - test_Xi
            d = np.linalg.norm(d)
#            phase_test = False
#            norm_test  = False
#            xi_norm = np.linalg.norm(xi)
#            if kz == 2.0 or kz == -2.0: 
#            #xi_norm >= 1e-10:
#                # Phase test
#                if phase_test:
#                    print("\nApplying phase shift at: " + str(kz) + "\n ")
#                    # shift the field by pi/3
#                    xi += 1.0j*(np.pi/3.0)
#
#                # Norm test
#                if norm_test:
#                    print("\nApplying norm doubling.")
#                    # Double the energy of the field
#                    xi *= np.sqrt(2.0) + 1.0j*np.sqrt(2.0)
#
#                # Fix xi to see if same coefficients are retrieved from projection
#                if fixXi:
#                    xi[0] = 1.0
#                    xi[1] = 0.0
#                    print("\nXI:")
#                    print("norm: " + str(xi_norm))
#                    print(xi)
#                    print("")
#                print(xi)

            #### Store the resolvent modes, amplitudesand coefficients
            # for reconstruction at a later date
            resolvent_modes_array[mx, mz, :, :r] = resolvent_modes[: , :r]
            coefficients_array[mx, mz, :r] = np.squeeze(np.asarray(xi[:r]))
            sing_vals_array[mx, mz, :r] = singular_values[:r]
            
            psi = np.asmatrix(resolvent_modes_array[mx,mz,:,:])
            sigma = np.asmatrix(np.diag(sing_vals_array[mx,mz,:]))
            xi = np.asmatrix(coefficients_array[mx,mz,:]).T

            test_full = psi * sigma * xi
            d2 = test_full - test_Xi
            d2 = np.linalg.norm(d2)
            a = d2
            
            
            
            

    calcTime = datetime.now() - startTime
    print("\n\n\n")
    print(calcTime)
    print("\n\n\n")
    deconstructed_dict = {}
    deconstructed_dict['resolvent_modes'] = resolvent_modes_array
    deconstructed_dict['singular_values'] = sing_vals_array
    deconstructed_dict['coefficients'] = coefficients_array

    return deconstructed_dict



 
def deconstruct_field_testing(original_ff_spectral,
                             kx_array,
                             kz_array,
                             Nm,
                             c,
                             Re,
                             baseflow,
                             mean_profile,
                             sparse,
                             fixXi):

    #================================================================
    #### Store the resolvent modes and amplitude coefficients 
    #    at each  Fourier mode pair
    #================================================================
    resolvent_modes_array = np.zeros((len(kx_array), len(kz_array), 3*Nm, 3*Nm), dtype=np.complex128)
    coefficients_array = np.zeros((len(kx_array), len(kz_array), 3*Nm), dtype=np.complex128)
    sing_vals_array = np.zeros((len(kx_array), len(kz_array), 3*Nm), dtype=np.float64)

    #================================================================
    #### Loop through wavenumbers 
    #================================================================
    chebyshev_differentiation, mean_flow_derivatives = ps.calculate_derivatives(Nm+2, mean_profile, baseflow)

    startTime = datetime.now()
    #================================================================
    #### Loop through wavenumbers 
    #================================================================
    for mx in range(0, len(kx_array)):
        kx = kx_array[mx]
        print('\n\nkx:'+ str(kx))

        for mz in range(0, len(kz_array)):
            kz  = kz_array[mz]
            sys.stdout.write(".")
            sys.stdout.flush()

            if kx == 0 or kz == 0: # Zero Fourier modes
                # Save zeroth modes to the resolvent modes...
                resolvent_modes_array[mx, mz, :, 0] = original_ff_spectral[mx,:,mz]
                continue # Start the loop again

            #------------------------------------------------
            #### Calculate the state vectors
            #------------------------------------------------
            omega = kx * c
            wegihted_transfer_function, w = ps.calculate_transfer_function(kx, kz, Re, Nm, omega, chebyshev_differentiation, mean_flow_derivatives)


            #------------------------------------------------
            #### Perform SVD
            #------------------------------------------------
            if sparse:
                vel_modes, singular_values, forcing_modes = svd(wegihted_transfer_function, full_matrices=False)
            else:
                vel_modes, singular_values, forcing_modes = svd(wegihted_transfer_function)

            vel_modes = np.asmatrix(vel_modes)
            forcing_modes = np.asmatrix(forcing_modes)

            #------------------------------------------------
            #### Check SVD
            #------------------------------------------------
            Tests.SVD(vel_modes, singular_values, forcing_modes, wegihted_transfer_function, sparse)


            #------------------------------------------------
            #### Retrieve non-grid-weighted resolvent modes (physical modes)
            #------------------------------------------------
            resolvent_modes = np.linalg.solve(w, vel_modes)

            #------------------------------------------------
            #### Check that the singular value matrix is invertible
            #------------------------------------------------
            S = np.diag(singular_values)
            Tests.invertible(S)
            
            
            #------------------------------------------------
            #### Check that the continuity condition is satisfied
            #------------------------------------------------
            Tests.continuity(resolvent_modes, S, kx, kz, Nm, chebyshev_differentiation['D1'])
#            elif sparse:
#                Tests.continuity(resolvent_modes[: , :r], kx, kz, Nm, state_vecs['D1'])

            #------------------------------------------------
            #### Check that the weighted resovlent and forcing modes are orthogonal.
            #------------------------------------------------
            Tests.orthogonality(vel_modes)
            Tests.orthogonality(forcing_modes)

            #------------------------------------------------
            #### Fix phase of resolvent modes based on critical layer or centreline
            #------------------------------------------------
            phase_shift = np.zeros((resolvent_modes.shape[1], resolvent_modes.shape[1]), dtype=np.complex128)
            ind0 = Nm/2 + 1 # Use centreline, unless
            if c < 1.0:
                inds = Tests.indices(mean_flow_derivatives['U'].diagonal(), lambda x: x > c)
                if len(inds) > 0:
                    ind0 = inds[0]

            np.fill_diagonal(phase_shift, np.exp(-1j * np.angle(resolvent_modes[ind0,:])))
            resolvent_modes *= phase_shift

            #------------------------------------------------
            #### Project resolvent modes to get amplitude coefficients
            #------------------------------------------------
            # denoted xi, defined as
            # chi  = singular_values * xi

            # Initialize the scalars vector
            xi = np.zeros((3*Nm, 1), dtype=np.complex128)    

            # Projection
            xi = inv(S) * resolvent_modes.H * w.H * w * np.asmatrix(original_ff_spectral[mx, :, mz]).T
            test_Xi = resolvent_modes * S * xi # this variable should be equal to the original vector
            o = np.asmatrix(original_ff_spectral[mx, :, mz]).T          
            d = o - test_Xi
            d = np.linalg.norm(d)
#            phase_test = False
#            norm_test  = False
#            xi_norm = np.linalg.norm(xi)
#            if kz == 2.0 or kz == -2.0: 
#            #xi_norm >= 1e-10:
#                # Phase test
#                if phase_test:
#                    print("\nApplying phase shift at: " + str(kz) + "\n ")
#                    # shift the field by pi/3
#                    xi += 1.0j*(np.pi/3.0)
#
#                # Norm test
#                if norm_test:
#                    print("\nApplying norm doubling.")
#                    # Double the energy of the field
#                    xi *= np.sqrt(2.0) + 1.0j*np.sqrt(2.0)
#
#                # Fix xi to see if same coefficients are retrieved from projection
#                if fixXi:
#                    xi[0] = 1.0
#                    xi[1] = 0.0
#                    print("\nXI:")
#                    print("norm: " + str(xi_norm))
#                    print(xi)
#                    print("")
#                print(xi)

            #### Store the resolvent modes, amplitudesand coefficients
            # for reconstruction at a later date
            resolvent_modes_array[mx, mz, :, :] = resolvent_modes
            coefficients_array[mx, mz, :] = np.squeeze(np.asarray(xi))
            sing_vals_array[mx, mz, :] = singular_values
            
            psi = np.asmatrix(resolvent_modes_array[mx,mz,:,:])
            sigma = np.asmatrix(np.diag(sing_vals_array[mx,mz,:]))
            xi = np.asmatrix(coefficients_array[mx,mz,:]).T

            test_full = psi * sigma * xi
            d2 = test_full - test_Xi
            d2 = np.linalg.norm(d2)
            a = d2
            
            
            
            

    calcTime = datetime.now() - startTime
    print("\n\n\n")
    print(calcTime)
    print("\n\n\n")
    deconstructed_dict = {}
    deconstructed_dict['resolvent_modes'] = resolvent_modes_array
    deconstructed_dict['singular_values'] = sing_vals_array
    deconstructed_dict['coefficients'] = coefficients_array

    return deconstructed_dict



def construct_field(resolvent_modes,
                    singular_values,
                    coefficients,
                    mean_ff_spectral,
                    kx_array,
                    kz_array,
                    Nm,
                    r):

    ''' 
    Construct a flow field using the resolvent modes, singular values and 
    amplitude coefficients.

    ================================================================
    INPUTS:
    ================================================================
    resolvent_modes:                4D array of truncated resolvent 
                                    modes at each Fourier mode 
                                    combination:
                                        resolvent_modes[mx, mz, :, :] gives 
                                        column vectors of given rank, i.e.
                                        rank = resolvent_modes.shape[3].

    singular_values:                3D array of truncated singular
                                    values at each Fourier mode 
                                    combination:
                                        singular_values[mx, mz, :] gives
                                        1D array of singular values, where
                                        rank = len(singular_values[mx, mz, :]).

    coefficients:                   3D array of truncated complex 
                                    amplitude coefficients at each 
                                    Fourier mode combination:
                                        coefficients[mx, mz, :] gives
                                        1D array of coefficients, where
                                        rank = len(coefficients[mx, mz, :]).

    mean_ff_spectral:               3D array of mean flow field  
                                    in xz spectral form with 
                                    Chebyshev nodes in wall-normal 
                                    direction (Chebyshev spacing).
                                    Stacked in wall-normal direction:
                                    dimensions: (Nx, Nd*Ny, Nz).

    kx_array:                       1D array of streamwise Fourier modes.

    kz_array:                       1D array of spanwise Fourier modes.

    Nm:                             Number of interior Chebyshev nodes,
                                    i.e. Ny - 2 (endpoints removed).

    r:                              Rank to approximate to.

    ================================================================
    OUTPUTS:
    ================================================================
    approximated_ff_spectral:       3D array of approximated flow field  
                                    in xz spectral form with 
                                    Chebyshev nodes in wall-normal 
                                    direction (Chebyshev spacing).
                                    Stacked in wall-normal direction:
                                    dimensions: (Nx, Nd*Ny, Nz).
    '''

    #================================================================
    #### Initialize empty 3D array to store approximated velocity field
    #================================================================
    # Approximated velocity field is stacked in the wall-normal direction
    approximated_ff_spectral  = np.zeros((len(kx_array), 3*Nm, len(kz_array)), dtype=np.complex128)

    for mx in range(0, len(kx_array)):
        kx = kx_array[mx]
        print('\n\nkx:'+ str(kx))

        for mz in range(0, len(kz_array)):
            kz  = kz_array[mz]
            sys.stdout.write(".")
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
            #### Construct approximated flow field
            #------------------------------------------------
            psi = np.asmatrix(resolvent_modes[mx,mz,:,:])
            sigma = np.asmatrix(np.diag(singular_values[mx,mz,:]))
            xi = np.asmatrix(coefficients[mx,mz,:]).T

            tmp = psi * sigma * xi
            approximated_ff_spectral[mx, :, mz] = np.squeeze(np.asarray(tmp))
    
    return approximated_ff_spectral




def construct_field_testing(resolvent_modes,
                            singular_values,
                            coefficients,
                            mean_ff_spectral,
                            kx_array,
                            kz_array,
                            Nm,
                            r,
                            original_field):

    #================================================================
    #### Initialize empty 3D array to store approximated velocity field
    #================================================================
    # Approximated velocity field is stacked in the wall-normal direction
    approximated_ff_spectral  = np.zeros((len(kx_array), 3*Nm, len(kz_array)), dtype=np.complex128)

    for mx in range(0, len(kx_array)):
        kx = kx_array[mx]
        print('\n\nkx:'+ str(kx))

        for mz in range(0, len(kz_array)):
            kz  = kz_array[mz]
            sys.stdout.write(".")
            sys.stdout.flush()

            if kx == 0 or kz == 0: # Zero Fourier modes
                #------------------------------------------------
                #### Set the zero Fourier modes to equal the mean flow
                #------------------------------------------------
                approximated_ff_spectral[mx, :, mz] = resolvent_modes[mx,mz,:,0]

                # [uvw] of approximation at the zero Fourier modes 
                # equals the [uvw] of mean at the zero Fourier modes.
                
                # Therefore the projection will be of the total flow,
                # i.e. original + mean flow

                continue # Start the loop again

            #------------------------------------------------
            #### Construct approximated flow field
            #------------------------------------------------
            psi = np.asmatrix(resolvent_modes[mx,mz,:,:r])
            sigma = np.asmatrix(np.diag(singular_values[mx,mz,:r]))
            xi = np.asmatrix(coefficients[mx,mz,:r]).T

            tmp = psi * sigma * xi
            
            a = original_field[mx,:,mz]
            d = a - tmp
            
            approximated_ff_spectral[mx, :, mz] = np.squeeze(np.asarray(tmp))
    
    return approximated_ff_spectral
