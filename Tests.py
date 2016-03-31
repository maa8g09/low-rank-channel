import numpy as np
import Utils as ut

def SVD(U, s, V, A, sparse):
    svd_passed = False
    if sparse:
        S = np.diag(s)
        svd_passed = np.allclose(A, np.dot(U, np.dot(S, V)))
        if not svd_passed:
            nrm = str(np.linalg.norm(np.dot( np.dot(U, S), V) - A))
            err = 'Something went wrong with the SVD, norm is ' + str(nrm)
            ut.error(err)

    elif not sparse:
        S = np.zeros((U.shape[1], V.shape[0]), dtype=complex)
        S[:V.shape[0], :V.shape[0]] = np.diag(s)
        svd_passed = np.allclose(A, np.dot(U, np.dot(S, V)))
        if not svd_passed:
            nrm = str(np.linalg.norm(np.dot( np.dot(U, S), V) - A))
            err = 'Something went wrong with the SVD, norm is ' + str(nrm)
            ut.error(err)

#    if np.linalg.norm(np.dot( np.dot(U, np.diag(S)), V) - A) >= 1e-9:
#        nrm = str(np.linalg.norm(np.dot( np.dot(U, np.diag(S)), V) - A))
#        err = 'Something went wrong with the SVD, norm is ' + str(nrm)
#        ut.error(err)
        
    return 0


def continuity(resolvent_modes, kx, kz, Nm, D1):
    for column in range(0, resolvent_modes.shape[1]):
        u = resolvent_modes[   0:Nm,   column]
        v = resolvent_modes[  Nm:2*Nm, column]
        w = resolvent_modes[2*Nm:3*Nm, column]
        continuty = 1.0j*kx*u + np.dot(D1, v) + 1.0j*kz*w
        norm = np.linalg.norm(continuty)
        if norm >= 1e-10:
            err = 'Something went wrong with the divergence criteria, norm is ' + str(norm)
            ut.error(err)
    
    return 0


def orthogonality(A):
    B = A * A.H
    I = np.eye(B.shape[0])
    C = B.real - I
    C_norm = np.linalg.norm(C)
    if C_norm >= 1e-10:
        err = 'Modes are not orthogonal, norm is ' + str(C_norm)
        ut.error(err)
    
    return 0


def indices(a, func):
    return [i for (i, val) in enumerate(a) if func(val)]


def invertible(A):
    # A must be a square matrix
    I = np.identity(A.shape[0])
    Z = I - (np.linalg.inv(A) * A)
    Znorm = np.linalg.norm(Z)
    if Znorm >= 1e-10:
        err = 'Matrix is not invertible, ||I - inv(A)A|| = ' + str(Znorm)
        ut.error(err)

    return 0


def checkHermitianSymmetry(velocityField, Nx, Nz):
    # Check for conjugacy in spectral flow field
    # otherwise known as Hermitian symmetry
    Nx = int(Nx/2+1)
    Nz = int(Nz/2+1)
    for kx in range(0, Nx):
        for kz in range(0, Nz):
            u_pos = velocityField[kx, :, kz]  # positive streamwise frequencies
            u_neg = velocityField[-kx, :, kz] # negative streamwise frequencies
            # Are they complex conjugates?
            delta_real = u_pos.real - u_neg.real
            delta_imag = u_pos.imag + u_neg.imag
            delta_real_norm = np.linalg.norm(delta_real)
            delta_imag_norm = np.linalg.norm(delta_imag)
            if delta_real_norm >= 1e-8 and delta_imag_norm >= 1e-8:
                # not conjugate
                print("Spectral velocity field has no complex conjugacy")
                print("between positive and negative streamwise frequencies at ")
                print("kx:\t"+str(kx))
                print("kz:\t"+str(kz))
                print(str(delta_real_norm))
                print(str(delta_imag_norm))
                print("")
    
            u_pos = velocityField[kx, :, kz]  # positive spanwise frequencies
            u_neg = velocityField[kx, :, -kz] # negative spanwise frequencies
            # Are they complex conjugates?
            delta_real = u_pos.real - u_neg.real
            delta_imag = u_pos.imag + u_neg.imag
            delta_real_norm = np.linalg.norm(delta_real)
            delta_imag_norm = np.linalg.norm(delta_imag)
            if delta_real_norm >= 1e-8 and delta_imag_norm >= 1e-8:
                # not conjugate
                print("Spectral velocity field has no complex conjugacy")
                print("between positive and negative spanwise frequencies at ")
                print("kx:\t"+str(kx))
                print("kz:\t"+str(kz))
                print(str(delta_real_norm))
                print(str(delta_imag_norm))
                print("")
                
            u_pos = velocityField[kx, :, kz]  # positive stream-spanwise frequencies
            u_neg = velocityField[-kx, :, -kz] # negative stream-spanwise frequencies
            # Are they complex conjugates?
            delta_real = u_pos.real - u_neg.real
            delta_imag = u_pos.imag + u_neg.imag
            delta_real_norm = np.linalg.norm(delta_real)
            delta_imag_norm = np.linalg.norm(delta_imag)
            if delta_real_norm >= 1e-8 and delta_imag_norm >= 1e-8:
                # not conjugate
                print("Spectral velocity field has no complex conjugacy")
                print("between positive and negative spanwise frequencies at ")
                print("kx:\t"+str(kx))
                print("kz:\t"+str(kz))
                print(str(delta_real_norm))
                print(str(delta_imag_norm))
                print("")
    
    return 0
    