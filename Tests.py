import numpy as np
import Utils as ut

def SVD(U, s, V, A, sparse):
    svd_passed = False
    if sparse:
        S = np.diag(s)
        svd_passed = np.allclose(A, np.dot(U, np.dot(S, V)))
        if not svd_passed:
            norm = str(np.linalg.norm(np.dot( np.dot(U, S), V) - A))
            err = 'Something went wrong with the SVD, norm is %.2E' % norm 
            ut.error(err)

    elif not sparse:
        S = np.zeros((U.shape[1], V.shape[0]), dtype=complex)
        S[:V.shape[0], :V.shape[0]] = np.diag(s)
        svd_passed = np.allclose(A, np.dot(U, np.dot(S, V)))
        if not svd_passed:
            norm = str(np.linalg.norm(np.dot( np.dot(U, S), V) - A))
            err = 'Something went wrong with the SVD, norm is %.2E' % norm 
            ut.error(err)

#    if np.linalg.norm(np.dot( np.dot(U, np.diag(S)), V) - A) >= 1e-9:
#        nrm = str(np.linalg.norm(np.dot( np.dot(U, np.diag(S)), V) - A))
#        err = 'Something went wrong with the SVD, norm is ' + str(nrm)
#        ut.error(err)
        
    return


def continuity(resolvent_modes, S, kx, kz, Nm, D1):
    projected_field = resolvent_modes * S
#    for column in range(0, resolvent_modes.shape[1]):
    u = projected_field[   0:Nm,   :]
    v = projected_field[  Nm:2*Nm, :]
    w = projected_field[2*Nm:3*Nm, :]
    continuty = 1.0j*kx*u + np.dot(D1, v) + 1.0j*kz*w
    norm = np.linalg.norm(continuty)
    if norm >= 1e-7:
        err = 'Something went wrong with the continuity condition, norm is %.2E' % norm 
        ut.error(err)
    
    return


def orthogonality(A):
    B = A * A.H
    I = np.eye(B.shape[0])
    C = B.real - I
    C_norm = np.linalg.norm(C)
    if C_norm >= 1e-10:
        err = 'Modes are not orthogonal, norm is ' + str(C_norm)
        ut.error(err)
    
    return


def indices(a, func):
    return [i for (i, val) in enumerate(a) if func(val)]


def invertible(A):
    # A must be a square matrix
    if A.shape[0] == A.shape[1]:
        I = np.identity(A.shape[0])
        Z = I - (np.linalg.inv(A) * A)
        norm = np.linalg.norm(Z)
        if norm >= 1e-10:
            err = 'Matrix is not invertible, ||I - inv(A)A|| = %.2E' % norm 
            ut.error(err)
    else:
        err = "Non-square matrix provided."
        ut.error(err)

    return


def checkHermitianSymmetry(velocityField, Nx, Nz):
    # Check for conjugacy in spectral flow field
    # otherwise known as Hermitian symmetry
    Nx = int(Nx/2+1)
    Nz = int(Nz/2+1)
#    for kx in range(0, Nx):
#        for kz in range(0, Nz):
#            u_pos = velocityField[kx, :, kz]  # positive spanwise frequencies
#            u_neg = velocityField[kx, :, -kz] # negative spanwise frequencies
#            # Are they complex conjugates?
#            delta_real = u_pos.real - u_neg.real
#            delta_imag = u_pos.imag + u_neg.imag
#            delta_real_norm = np.linalg.norm(delta_real)
#            delta_imag_norm = np.linalg.norm(delta_imag)
#            if delta_real_norm >= 1e-8 and delta_imag_norm >= 1e-8:
#                # not conjugate
#                print("Spectral velocity field has no complex conjugacy")
#                print("between positive and negative spanwise frequencies at ")
#                print("kx:\t"+str(kx))
#                print("kz:\t"+str(kz))
#                print(str(delta_real_norm))
#                print(str(delta_imag_norm))
#                print("")
#                
#            u_pos = velocityField[kx, :, kz]  # positive stream-spanwise frequencies
#            u_neg = velocityField[-kx, :, -kz] # negative stream-spanwise frequencies
#            # Are they complex conjugates?
#            delta_real = u_pos.real - u_neg.real
#            delta_imag = u_pos.imag + u_neg.imag
#            delta_real_norm = np.linalg.norm(delta_real)
#            delta_imag_norm = np.linalg.norm(delta_imag)
#            if delta_real_norm >= 1e-8 and delta_imag_norm >= 1e-8:
#                # not conjugate
#                print("Spectral velocity field has no complex conjugacy")
#                print("between positive and negative spanwise frequencies at ")
#                print("kx:\t"+str(kx))
#                print("kz:\t"+str(kz))
#                print(str(delta_real_norm))
#                print(str(delta_imag_norm))
#                print("")

    for kz in range(0, Nz):
        u_pos = velocityField[1:, :, kz]  # positive spanwise frequencies
        u_neg = velocityField[1:, :, -kz] # negative spanwise frequencies
        # Are they complex conjugates?
        u_posR = u_pos.real
        u_negR = u_neg.real
        
        u_posI = u_pos.imag
        u_negI = u_pos.imag
        
        delta_real = u_posR + u_negR
        delta_imag = u_posI - u_negI
        
        delta_real_norm = np.linalg.norm(delta_real)
        delta_imag_norm = np.linalg.norm(delta_imag)
        if delta_real_norm >= 1e-8 and delta_imag_norm >= 1e-8:
            # not conjugate
            print("Spectral velocity field has no complex conjugacy")
            print("between positive and negative spanwise frequencies at ")
            print("kz:\t"+str(kz))
            print(str(delta_real_norm))
            print(str(delta_imag_norm))
            print("")

    return
    
    
def fft_ifft(A):
    original = A.velocityField

    # Remove wall boundaries
    A.remove_wall_boundaries()
    
    # FFT
    A.make_xz_spectral()
    
    # Stack
    A.stack_ff_in_y()
    
    # Unstack
    A.unstack_ff()
    
    # IFFT
    A.make_xz_physical()
    
    # Add wall boundaries
    A.add_wall_boundaries()
    
    # The difference between the flow field after having gone through
    # above operations and the original should be the same.
    diff = np.linalg.norm(original - A.velocityField.real)
    if diff >= 1e-12:
        message = "FFT <=> IFFT test failed."
        ut.error(message)
    else:
#        print("Passed FFT <=> IFFT test.\nThe norm of the difference is %.2E \n" % difference)
        print("")
    return


def difference(A, B, tolerance, message):
    diff = np.linalg.norm(A - B)
    if diff >= tolerance:
        print(message)
        err = "Matrices are not equal, ||delta|| = %.2E" % diff 
        ut.error(err)
    return
    
def projection(chi, sigma, psi, u_hat, tolerance):
    test = psi * sigma * chi    
    d = u_hat - test
    d = np.linalg.norm(d)
    if d >= tolerance:
        err = "Projection failed."
        ut.error(err)
    return
    