import numpy as np
import time
import FlowField as ffClass
import ChannelResolvent as cr
import Utils as ut

def main(wp, Nx, Ny, Nz, Lx, Lz, Re, c, th, bf):
    #====================================================================
    #### Construct an instance of the FlowFieldGeometry class
    #====================================================================
    tmp = np.zeros((3, Nx, Ny, Nz),dtype=float)
    ff = ffClass.FlowField(Nx,
                            Ny,
                            Nz,
                            Lx,
                            Lz,
                            wp,
                            c,
                            th,
                            Re,
                            bf,
                            tmp,
                            "pp")
    #================================================================
    #### Generate flow field using resolvent formulation
    #================================================================
    velocity_field = cr.resolvent_formulation(ff)
    ff.set_ff(velocity_field, "pp")
    #================================================================
    #### ---- Unstack velocity field in the wall-normal direction
    #================================================================
    ff.unstack_ff(ff.velocityField)
    #================================================================
    #### ---- Add wall boundaries
    #================================================================
    ff.add_wall_boundaries()

    u = ff.velocityField[0,:,:,:].real
    v = ff.velocityField[1,:,:,:].real
    w = ff.velocityField[2,:,:,:].real
    
    #================================================================
    # End of file
    #================================================================
    ut.print_EndMessage()
    return 0






wavepacket="KB"
nx=36
ny=65
nz=36
lx=2.0*np.pi
lz=2.0*np.pi/3.0
re=400
c=0.5
th=0.0
bf="lam"
main("KB",nx,ny,nz,lx,lz,re,c,th,bf)
# NX MUST EQUAL NZ