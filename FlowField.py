####################################################################################################
## This is an editable file where the settings for the simulations can be set.                    ##
####################################################################################################
import numpy as np
import Wavepackets as wp

class FlowFieldGeometry(object):
    def __init__(self, baseflow, packet, Nd, Nx, Ny, Nz, Re, c, theta):
        """
        A class representing a flow field's geometric attributes (The flow field is not generated).
        """
        self.baseflow = baseflow
        # Wave number triplet_______________________________________________________________________
        # The mode combinations are taken from Sharma & McKeon (2013) paper
        self.wavepacket = packet
        xlabel = packet+'_x'
        zlabel = packet+'_z'
        kx = wp.wavepackets[xlabel]
        self.kx = kx
        kz = wp.wavepackets[zlabel]
        self.kz = kz

        # Amplitudes________________________________________________________________________________
        alabel = packet+'_a'
        chi = wp.wavepackets[alabel]
        self.chi = chi
        self.chi_tilde = 10.0**(theta) * chi
        self.theta = theta

        # Resolution in wall-normal direction_______________________________________________________
        self.Nd = Nd
        self.Nx = Nx
        self.Ny = Ny
        self.Nz = Nz
        lamda_x = 2.0*np.pi / kx[0]
        lamda_z = 2.0*np.pi / kz[0]
        required_Lx = 4.0 * np.pi
        required_Lz = 2.0 * np.pi / 3
#        required_Lx = lamda_x
#        required_Lz = lamda_z
        cx = required_Lx / lamda_x
        cz = required_Lz / lamda_z
        self.Lx = cx * lamda_x # Length is 4pi (twice the streamwise wavelength)
        self.Lz = cz * lamda_z # Length is 2pi/3 (twice the spanwise wavelength)

        # Stationary nodes along each axis__________________________________________________________
        # X axis
        # Full range of Fourier modes
        self.Mx_full = np.arange((-Nx/2.0), (Nx/2.0)+1)
        # 1 harmonic 
        self.Mx = np.arange(-1.0, 2.0)

        # Z axis____________________________________________________________________________________
        # Full range of Fourier modes
        self.Mz_full = np.arange(0, 0.5*(Nz) + 1.0)
        # 1 harmonic
        self.Mz = np.arange(2.0)

        # X & Z axes________________________________________________________________________________
        self.x = np.linspace(0.0, self.Lx, Nx)
        self.y_uniform = np.linspace(1.0, -1.0, Ny)
        self.y = np.linspace(1.0, -1.0, Ny)
        self.z = np.linspace(-self.Lz/2.0, self.Lz/2.0, Nz)

        # Total number of modes_____________________________________________________________________
        self.modes = self.Ny - 2 # take away the first and last grid points
        self.t = 0

        # Reynolds number___________________________________________________________________________
        self.Re = Re

        # Wave speed________________________________________________________________________________
        self.c = c


    def set_y(self, y_cheb_full):
        self.y = y_cheb_full

    def set_Ny(self, m):
        self.Ny = m

class FlowField(FlowFieldGeometry):
    def __init__(self, ffgeom, ff, state):
        FlowFieldGeometry.__init__(self,
                                   ffgeom.baseflow,
                                   ffgeom.wavepacket,
                                   ffgeom.Nd,
                                   ffgeom.Nx,
                                   ffgeom.Ny,
                                   ffgeom.Nz,
                                   ffgeom.Re,
                                   ffgeom.c,
                                   ffgeom.theta)
        self.velocityField = ff
        self.state = state


class FlowFieldChannelFlow(object):
    def __init__(self, Nd, Nx, Ny, Nz, Lx, Lz, alpha, beta, c, baseflow, Re, ff, state):
        """
        A class representing a flow field.
        """
        self.Nd = Nd
        self.Nx = Nx
        self.Ny = Ny
        self.modes = Ny - 2
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
        
        self.Mz = np.arange(0, np.floor(Nz/2) + 1)
        
        self.alpha = alpha
        self.beta = beta
        self.c = c
        self.baseflow = baseflow
        self.Re = Re
        self.velocityField = ff
        self.state = state

        # X & Z axes________________________________________________________________________________
        self.x = np.linspace(0.0, self.Lx, Nx)
        self.y = np.linspace(1.0, -1.0, Ny)
        for ny in range(0, Ny):
            self.y[ny] = np.cos(ny * np.pi/ (Ny-1) )
        self.z = np.linspace(-self.Lz/2.0, self.Lz/2.0, Nz)
        
    def set_modes(self, m):
        self.modes = m

    def set_ff(self, ff, state):
        self.velocityField = ff
        self.state = state

    def set_Ny(self, Ny):
        self.Ny = Ny
        self.modes = Ny - 2

    def set_rank(self, rank):
        self.rank = rank

class FlowFieldChannelFlow2(object):
    def __init__(self, Nd, Nx, Ny, Nz, Lx, Lz, alpha, beta, c, baseflow, Re, ff, state):
        """
        A class representing a flow field.
        """
        self.Nd = Nd
        self.Nx = Nx
        self.Ny = Ny
        self.modes = Ny - 2
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
        self.state = state

        self.x = np.linspace(0.0, self.Lx, Nx)
        self.y = np.linspace(1.0, -1.0, Ny)
        for ny in range(0, Ny):
            self.y[ny] = np.cos(ny * np.pi/ (Ny-1) )
        self.z = np.linspace(-self.Lz/2.0, self.Lz/2.0, Nz)
        
    def set_modes(self, m):
        self.modes = m

    def set_ff(self, ff, state):
        self.velocityField = ff
        self.state = state

    def set_Ny(self, Ny):
        self.Ny = Ny
        self.modes = Ny - 2

    def set_rank(self, rank):
        self.rank = rank

    def make_xz_spectral(self):
        self.velocityField = np.fft.fft(self.velocityField, axis=3)
        self.velocityField = np.fft.fft(self.velocityField, axis=1)
        self.state = "sp"
        # shift the flow field 

    def make_xz_physical(self):
        self.velocityField = np.fft.ifft(self.velocityField, axis=1)
        self.velocityField = np.fft.ifft(self.velocityField, axis=3)
        self.state = "pp"