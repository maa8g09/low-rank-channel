import numpy as np
import Utils as ut


def SVDNorm(U, S, V, A):
    if np.linalg.norm(np.dot( np.dot(U, np.diag(S)), V) - A) >= 1e-10:
        nrm = str(np.linalg.norm(np.dot( np.dot(U, np.diag(S)), V) - A))
        err = 'Something went wrong with the SVD, norm is ' + str(nrm)
        ut.error(err)
        
    return 0
    
    
def divergence(resolvent_modes, alpha, beta, m, D1):
    divgnce = 1.0j*alpha*resolvent_modes[0:m, 0] + np.dot(D1, resolvent_modes[m:2*m, 0]) + 1.0j*beta*resolvent_modes[2*m:3*m, 0]
    div_norm = np.linalg.norm(divgnce)
    
    if div_norm >= 1e-10:
        err = 'Something went wrong with the divergence criteria, norm is ' + str(div_norm)
        ut.error(err)
    
    return 0
    
def orthogonality(A):
    B = A * A.H
    I = np.eye(B.shape[0])
    C = B.real - I
    C_norm = np.linalg.norm(C)
    if C_norm >= 1e-10:
        err = 'Resolvent modes are not orthogonal, norm is ' + str(C_norm)
        ut.error(err)
    
    