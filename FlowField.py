import numpy as np
from Wavepackets import wavepackets_dictionary

class FlowField(object):
    def __init__(self, Nx, Ny, Nz, Lx, Lz, packet, c, theta, Re, baseflow,  ff, state):
        """
            A class representing a flow field generated using the 
            Resolvent Formulation. 
        """
        #============================================================
        #### Geometry (Physical)
        self.Nd = 3
        self.Nx = Nx    # Number of grid point in the streamwise direction
        self.Ny = Ny    # Number of grid point in the wall-normal direction
        self.Nz = Nz    # Number of grid point in the spanwise direction
        self.Lx = Lx    # Length of domain
        self.Ly = 2.0   # Height of domain
        self.Lz = Lz    # Width of domain
        self.x = np.linspace(0.0, self.Lx, self.Nx)
        self.y_uniform = np.linspace(1.0, -1.0, self.Ny)
        self.y = np.linspace(1.0, -1.0, self.Ny)
        for ny in range(0, Ny): # Chebyshev nodes
            self.y[ny] = np.cos(ny * np.pi/ (Ny-1))
        self.z = np.linspace(0.0, self.Lz, self.Nz)
        #============================================================
        #### Geometry (Fourier)
        # All Fourier modes
        Mx_tmp = np.arange(-np.ceil(Nx/2)+1, np.floor(Nx/2)+1)
        self.Mx_full = np.zeros(Nx)
        if Nx % 2 == 0:
            # even Nx
            self.Mx_full[:np.ceil(Nx/2)+1] = Mx_tmp[np.floor(Nx/2)-1:]
            self.Mx_full[-np.ceil(Nx/2)+1:] = Mx_tmp[:np.ceil(Nx/2)-1] 
        else:
            # odd Nx
            self.Mx_full[:np.ceil(Nx/2)] = Mx_tmp[np.floor(Nx/2):]
            self.Mx_full[-np.floor(Nx/2):] = Mx_tmp[:np.ceil(Nx/2)-1]
        Mz_tmp = np.arange(-np.ceil(Nz/2)+1, np.floor(Nz/2)+1)
        self.Mz_full = np.zeros(Nz)
        if Nz % 2 == 0:
            # even Nz
            self.Mz_full[:np.ceil(Nz/2)+1] = Mz_tmp[np.floor(Nz/2)-1:]
            self.Mz_full[-np.ceil(Nz/2)+1:] = Mz_tmp[:np.ceil(Nz/2)-1] 
        else:
            # odd Nz
            self.Mz_full[:np.ceil(Nz/2)] = Mz_tmp[np.floor(Nz/2):]
            self.Mz_full[-np.floor(Nz/2):] = Mz_tmp[:np.ceil(Nz/2)-1]
        # First Fourier modes
        self.Mx = np.arange(-1.0, 2.0)
        self.Mz = np.arange(-1.0, 2.0)
        # Chebyshev nodes (Resolvent/Forcing mode column length)
        self.Nm = Ny - 2
        #============================================================
        #### Wavepacket
        self.Wavepacket = packet
        self.kx  = wavepackets_dictionary[packet]['x']
        self.kz  = wavepackets_dictionary[packet]['z']
        self.chi = wavepackets_dictionary[packet]['a']
        self.theta = theta
        self.chi_tilde = 10.0**(self.theta) * self.chi
        #============================================================
        #### Flow details
        self.baseflow = baseflow
        self.Re = Re
        self.c = c
        self.t = 0
        self.velocityField = ff
        self.state = state

    def set_ff(self, ff, state):
        self.velocityField = ff
        self.state = state

    def make_xz_spectral(self):
        self.velocityField = np.fft.fft(self.velocityField, axis=3)
        self.velocityField = np.fft.fft(self.velocityField, axis=1)
        self.state = "sp"

    def make_xz_physical(self):
        self.velocityField = np.fft.ifft(self.velocityField, axis=1)
        self.velocityField = np.fft.ifft(self.velocityField, axis=3)
        self.state = "pp"


    def unstack_ff(self, tmp):
        self.is_stacked_in_y = False
        self.velocityField = self.velocityField.reshape((self.Nd,
                                                         self.Nx,
                                                         self.Nm,
                                                         self.Nz))
        self.velocityField[0, :, :, :] = tmp[:,         0:self.Nm  , :].real
        self.velocityField[1, :, :, :] = tmp[:,   self.Nm:self.Nm*2, :].real
        self.velocityField[2, :, :, :] = tmp[:, 2*self.Nm:self.Nm*3, :].real


    def add_wall_boundaries(self):
        if self.velocityField.ndim == 4:
            wall_boundary = np.zeros((self.Nd, self.Nx, 1, self.Nz)) # no-slip boundary condition
            self.velocityField = np.concatenate((wall_boundary[:,:,:,:],
                                                 self.velocityField[:,:,:,:],
                                                 wall_boundary[:,:,:,:]),
                                                 axis=2)











class FlowFieldChannelFlow(object):
    def __init__(self, Nd, Nx, Ny, Nz, Lx, Lz, alpha, beta, c, baseflow, Re, ff, state):
        """
        A class representing a flow field.
        """
        self.Nd = Nd
        self.Nx = Nx
        self.Ny = Ny
        self.Nm = Ny - 2
        self.Nz = Nz
        self.Lx = Lx
        self.Lz = Lz

        Mx_tmp = np.arange(-np.ceil(Nx/2)+1, np.floor(Nx/2)+1)
        self.Mx = np.zeros(Nx)
        if Nx % 2 == 0:
            # even Nx
            self.Mx[:np.ceil(Nx/2)+1] = Mx_tmp[np.floor(Nx/2)-1:]
            self.Mx[-np.ceil(Nx/2)+1:] = Mx_tmp[:np.ceil(Nx/2)-1] 
        else:
            # odd Nx
            self.Mx[:np.ceil(Nx/2)] = Mx_tmp[np.floor(Nx/2):]
            self.Mx[-np.floor(Nx/2):] = Mx_tmp[:np.ceil(Nx/2)-1]

        Mz_tmp = np.arange(-np.ceil(Nz/2)+1, np.floor(Nz/2)+1)
        self.Mz = np.zeros(Nz)
        if Nz % 2 == 0:
            # even Nz
            self.Mz[:np.ceil(Nz/2)+1] = Mz_tmp[np.floor(Nz/2)-1:]
            self.Mz[-np.ceil(Nz/2)+1:] = Mz_tmp[:np.ceil(Nz/2)-1] 
        else:
            # odd Nz
            self.Mz[:np.ceil(Nz/2)] = Mz_tmp[np.floor(Nz/2):]
            self.Mz[-np.floor(Nz/2):] = Mz_tmp[:np.ceil(Nz/2)-1]


        self.alpha = alpha
        self.beta = beta
        self.c = c
        self.baseflow = baseflow
        self.Re = Re
        self.velocityField = ff

        if self.velocityField.ndim == 4:
            self.is_stacked_in_y = False

        elif self.velocityField.ndim == 3:
            self.is_stacked_in_y = True

        self.state = state
        
        self.x = np.linspace(0.0, self.Lx, Nx)
        self.z = np.linspace(0.0, self.Lz, Nz)
        self.y = np.linspace(1.0, -1.0, Ny) # uniform grid spacing
        for ny in range(0, Ny): # Chebyshev nodes
            self.y[ny] = np.cos(ny * np.pi/ (Ny-1))




    def set_ff(self, ff, state):
        self.velocityField = ff
        self.state = state




    def make_xz_spectral(self):
        if self.velocityField.ndim == 4:
            self.velocityField = np.fft.fft(self.velocityField, axis=3)
            self.velocityField = np.fft.fft(self.velocityField, axis=1)
            self.state = "sp"
        elif self.velocityField.ndim == 3:
            self.velocityField = np.fft.fft(self.velocityField, axis=2)
            self.velocityField = np.fft.fft(self.velocityField, axis=0)
            self.state = "sp"




    def make_xz_physical(self):
        if self.velocityField.ndim == 4:
            self.velocityField = np.fft.ifft(self.velocityField, axis=1)
            self.velocityField = np.fft.ifft(self.velocityField, axis=3)
            self.velocityField = self.velocityField.real
            self.state = "pp"
        elif self.velocityField.ndim == 3:
            self.velocityField = np.fft.ifft(self.velocityField, axis=0)
            self.velocityField = np.fft.ifft(self.velocityField, axis=2)
            self.state = "sp"




    def stack_ff_in_y(self):
        self.is_stacked_in_y = True
#        self.velocityField = np.concatenate((self.velocityField[0, :, :, :],
#                                             self.velocityField[1, :, :, :],
#                                             self.velocityField[2, :, :, :]),
#                                             axis=1)
        self.velocityField = self.velocityField.reshape((self.Nx,
                                                         self.Ny*self.Nd,
                                                         self.Nz))



    def unstack_ff(self):
        self.is_stacked_in_y = False
        self.velocityField = self.velocityField.reshape((self.Nd,
                                                         self.Nx,
                                                         self.Ny,
                                                         self.Nz))




    def remove_wall_boundaries(self):
        if self.velocityField.ndim == 4:
            self.Ny -= 2 # end points removed
            self.numModes = self.Ny # since end points removed, Nm = Ny
            self.velocityField = self.velocityField[:, :, 1:-1, :]



    def add_wall_boundaries(self):
        if self.velocityField.ndim == 4:
            self.Ny += 2
            self.numModes = self.Ny-2
            if self.baseflow == "lam": # Plane Poiseuille base flow
                wall_boundary = np.zeros((self.Nd, self.Nx, 1, self.Nz)) # no-slip boundary condition 
                self.velocityField = np.concatenate((wall_boundary[:,:,:,:],
                                                     self.velocityField[:,:,:,:],
                                                     wall_boundary[:,:,:,:]),
                                                     axis=2)

            elif self.baseflow == "cou": # Couette base flow
                top_boundary = np.zeros((self.Nd, self.Nx, 1, self.Nz))
                bot_boundary =top_boundary # no-slip boundary condition
                self.velocityField = np.concatenate((top_boundary[:,:,:,:],
                                                     self.velocityField[:,:,:,:],
                                                     bot_boundary[:,:,:,:]),
                                                     axis=2)


    def make_real(self):
        for nx in range(self.Nx):
            for ny in range(0, self.Ny):
                for nz in range(0, self.Nz):
                    for i in range(0, self.Nd):
                        tmp = self.velocityField[i, nx, ny, nz]
                        tmp = np.linalg(tmp)
                        self.velocityField[i, nx, ny, nz] = tmp
