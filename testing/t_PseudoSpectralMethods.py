import PseudoSpectralMethods as ps
from numpy.linalg import svd
from numpy import diag


alpha = 5.0                     # Streamwise Fourier mode
beta  = 1.0                     # Spanwise Fourier mode
c = 0.5                         # Wavespeed
omega = c * alpha       # Temporal frequency
Re = 1000.0                     # Reynolds number
N  = 31                         # Wall-normal resolution
Nm = N-2                        # Wall-normal resolution (endpoints removed)
bf = "lam"                      # Baseflow
vel_profile = []                # Mean velocity profile (if known)



chebyshev_differentiation, mean_flow_derivatives = ps.calculate_derivatives(N, vel_profile, bf)

H, w = ps.calculate_transfer_function(alpha, beta, Re, Nm, omega, chebyshev_differentiation, mean_flow_derivatives)

U, S, V = svd(H)
S = diag(S)
print("Fin")