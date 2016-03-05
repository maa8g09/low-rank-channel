import sys
import os
import numpy as np
#import pylab 
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
#from matplotlib import animation
from matplotlib import cm as cm


def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


def make_FlowField_output_directory(output_directory, flowFieldGeometry, date):
    # Make sure that the trailing character is a forward slash
    if output_directory[-1] != '/':
        output_directory += '/'

    # Create the output directory name
    # It is based on the date
    main_direc = 'Re' + str(flowFieldGeometry.Re) + '/' + flowFieldGeometry.wavepacket + '/' + date + '/theta_' + format(flowFieldGeometry.theta, '.4f') + '/'
    output_directory += main_direc

    # if directory doesn't exist, create it
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
        print('\nMade the following directory:')
        print(output_directory,)

    return output_directory


def make_FlowField_output_directory_wIteration(output_directory, flowFieldGeometry, date, iteration):
    # Make sure that the trailing character is a forward slash
    if output_directory[-1] != '/':
        output_directory += '/'

    # Create the output directory name
    # It is based on the date
    main_direc = 'Re' + str(flowFieldGeometry.Re) + '/' + flowFieldGeometry.wavepacket + '/' + date + '/' + str(iteration).zfill(3) +'_theta_' + format(flowFieldGeometry.theta, '.4f') + '/'
    output_directory += main_direc

    # if directory doesn't exist, create it
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
#        print('\nMade the following directory:')
#        print(output_directory)

    return output_directory


def read_ASC_PP(directory, fileName):

    var = read_GEOM(directory, fileName)

    file = directory + fileName + "_pp.asc"
    file = open(file, 'r')

    # Create an empty 4D array to store the velocity data
    U = np.zeros((var['Nd'],
                  var['Nx'],
                  var['Ny'],
                  var['Nz']),
                  dtype=np.float128)

    for i, line in enumerate(file):
        values = line.split()
        nd = int(values[0])
        nx = int(values[1])
        ny = int(values[2])
        nz = int(values[3])
        coeff = float(values[4])

        U[nd, nx, ny, nz] = coeff

    file.close()

#    U = np.concatenate((U[0,:,:,:],
#                        U[1,:,:,:],
#                        U[2,:,:,:]), 
#                        axis=1)
    var['ff'] = U

    return var


def read_ASC_SP(directory, fileName):

    var = read_GEOM(directory, fileName)

    file = directory + fileName + "_sp.asc"
    file = open(file, 'r')

    # Create an empty 4D array to store the velocity data
    U_hat = np.zeros((var['Nd'],
                      var['Nx'],
                      var['Ny'],
                      np.floor(var['Nz'] / 2.0) + 1.0),
                      dtype=np.complex128)

    for i, line in enumerate(file):
        values = line.split()
        nd = int(values[0])
        mx = int(values[1])
        ny = int(values[2])
        mz = int(values[3])
        coeff = complex(values[4])

        U_hat[nd, mx, ny, mz] = coeff

    file.close()

    U_hat = np.concatenate((U_hat[0,:,:,:],
                            U_hat[1,:,:,:],
                            U_hat[2,:,:,:]),
                            axis=1)
    var['ff'] = U_hat

    return var


def read_Details(directory, fileName):
    file = directory + fileName + "_Details.txt"
    file = open(file, 'r')

    var = {}

    # i, kx, y, kz
    for i, line in enumerate(file):
        values = line.split()
        if len(values) != 0:
            if values[0] == 'Nx:':
                var['Nx'] = int(values[1])
    
            elif values[0] == 'Ny:':
                var['Ny'] = int(values[1])
    
            elif values[0] == 'Nz:':
                var['Nz'] = int(values[1])
    
            elif values[0] == 'Nd:':
                var['Nd'] = int(values[1])
    
            elif values[0] == 'Lx:':
                var['Lx'] = float(values[1])
    
            elif values[0] == 'Lz:':
                var['Lz'] = float(values[1])
    
            elif values[0] == 'Re:':
                var['Re'] = float(values[1])
    
            elif values[0] == 'c:':
                var['c'] = float(values[1])
    
            elif values[0] == 'wp:':
                var['wp'] = str(values[1])
    
            elif values[0] == 'bf:':
                var['bf'] = str(values[1])
    
            elif values[0] == 'theta:':
                var['theta'] = float(values[1])

    file.close()

    return var


def read_FF(directory, fileName):
    command = "field2ascii -p " + fileName+".ff " + fileName
    print("\n",command,"\n")
    os.system(command)

    var = read_ASC_PP(directory, fileName)

    command = "rm " + fileName + "_pp.asc"
    command  += " " + fileName + "_pp_cpp.asc"
    command  += " " + fileName + "_sp.asc"
    command  += " " + fileName + "_sp_cpp.asc"
    command  += " " + fileName + ".geom"
    os.system(command)
    
    return var

def read_GEOM(directory, fileName):
    
    file = directory + fileName + ".geom"
    file = open(file, 'r')

    var = {}

    # i, kx, y, kz
    for i, line in enumerate(file):
        values = line.split()

        if values[1] == '%Nx':
            var['Nx'] = int(values[0])

        elif values[1] == '%Ny':
            var['Ny'] = int(values[0])

        elif values[1] == '%Nz':
            var['Nz'] = int(values[0])

        elif values[1] == '%Nd':
            var['Nd'] = int(values[0])

        elif values[1] == '%Lx':
            var['Lx'] = float(values[0])

        elif values[1] == '%Lz':
            var['Lz'] = float(values[0])

        elif values[1] == '%alpha=2pi/Lx':
            var['alpha'] = float(values[0])

        elif values[1] == '%gamma=2pi/Lz':
            var['beta'] = float(values[0])

    file.close()

    return var

def read_Output_DNS(fileName, T0, T1):    
    '''
    Sample output, you can decide which variables you want ot plot from here.
    
               t == 42.8
       L2Norm(u) == 0.704647
    chebyNorm(u) == 1.0705
     dissip(u+U) == 2.59421
      input(u+U) == 1.00164
      divNorm(u) == 0.0407906
           Ubulk == 1.33333
           ubulk == 0.666667
            dPdx == -0.0040663
      Ubulk*h/nu == 1600
     Uparab*h/nu == 2400
             CFL == 0.0220104
    '''


    file = open(fileName, 'r')

    data = {}
    data['t'] = []
    data['L2Norm(u)'] = []
#    data['dissip(u+U)'] = []
#    data['CFL'] = []

    for k, line in enumerate(file):
        values = line.split()

        if len(values) == 0:
            print("Reading...")

        else:
            if values[0] == 't':
                if values[2] != '-nan' or values[2] != 'nan' or values[2] != 'done!':
                    data['t'].append(float(values[2]))

            elif values[0] == 'L2Norm(u)':
                if values[2] != '-nan' or values[2] != 'nan' or values[2] != 'done!':
                    data['L2Norm(u)'].append(float(values[2]))
                    continue

#            elif values[0] == 'dissip(u+U)':
#                if values[2] != '-nan' or values[2] != 'nan' or values[2] != 'done!':
#                    data['dissip(u+U)'].append(float(values[2]))
#
#            elif values[0] == 'CFL':
#                if values[2] != '-nan' or values[2] != 'nan' or values[2] != 'done!':
#                    data['CFL'].append(float(values[2]))

    file.close()

    data['t']           = np.asarray(data['t'])
    data['L2Norm(u)']   = np.asarray(data['L2Norm(u)'])
#    data['dissip(u+U)'] = np.asarray(data['dissip(u+U)'])
#    data['CFL']         = np.asarray(data['CFL'])

    T0 = T0*2 +1
    T1 = T1*2 +1

    data['t']         = data['t'][T0:T1]
    data['L2Norm(u)'] = data['L2Norm(u)'][T0:T1]

    return data


def read_Grid(directory, fileName):
    file = open(fileName, 'r')
    aray = []

    for j, line in enumerate(file):
        if j!=0: # Don't read the first line.
            values = line.split()
            aray.append(float(values[0]))

    aray = np.asanyarray(aray)
    
    return aray


def read_Vel_Profile(fileName):
    file = open(fileName, 'r')
    aray = []

    for j, line in enumerate(file):
        values = line.split()
        aray.append(float(values[0]))

    aray = np.asanyarray(aray)

    return aray


def write_approximated_ASC(flowField, output_directory, rank):

    fileName = "u_rank.asc"
    file = open(output_directory + fileName, "w")

    # Write it in 2 chunks
    # first write from 0 -> Nx/2 + 1 
    print("Mx and Mz lengths")
    print(str(len(flowField.Mx)))
    print(str(len(flowField.Mz)))

    for mx in range(0, len(flowField.Mx)):
        for ny in range(0, flowField.Ny):
            for mz in range(0, len(flowField.Mz)):
                for nd in range(0, flowField.Nd):

#                    ikx = flowField.Mx[mx]
#                    ikz = flowField.Mz[mz]
#                    print(str(mx) + "\t\t" + str(mz))

                    tmpReal = flowField.velocityField[nd, mx, ny, mz].real
                    tmpImag = flowField.velocityField[nd, mx, ny, mz].imag

                    output = '(' + format( tmpReal, '.16f') + ', ' + format( tmpImag, '.16f') + ')'

                    file.write(output + "\n")

    file.close()

    print('\nSaved approximated ASCII file')

#
#    fileName = "python_"+str(rank)+".txt"
#    file = open(output_directory + fileName, "w")
#
#
#    for mx in range(0, len(flowField.Mx)):
#        for ny in range(0, flowField.Ny):
#            for mz in range(0, len(flowField.Mz)):
#                for nd in range(0, flowField.Nd):
#
#                    ikx = flowField.Mx[mx]
#                    ikz = flowField.Mz[mz]
#
#                    tmpReal = flowField.ff[nd, mx, ny, mz].real
#                    tmpImag = flowField.ff[nd, mx, ny, mz].imag
#
#                    output = str(nd) + '\t'
#                    output+= str(mx) + '\t'
#                    output+= str(ny) + '\t'
#                    output+= str(mx) + '\t'
#                    output+= '(' + format( tmpReal, '.16f') + ', ' + format( tmpImag, '.16f') + ')'
#                    output+= '\t' + str(ikx) + '\t'
#                    output+= str(ikz)
#
#                    file.write(output + "\n")
#
#    file.close()

    return fileName


def write_ASC(flowField, output_directory, fileName):
    fileName += '.asc'
    file = open(output_directory + fileName, "w")

    for nx in range(0, flowField.Nx):
        for ny in range(0, flowField.Ny):
            for nz in range(0, flowField.Nz):
                for nd in range(0, flowField.Nd):
                    tmp = flowField.velocityField[nd, nx, ny, nz]
                    tmp = format(tmp, '.16f')
                    file.write(tmp + "\n")

    file.close()

    print('\nSaved ASCII file')

    return 0


def write_ASC_Py(flowField, output_directory, fileName):
    fileName += '_pp.asc'
    file = open(output_directory + fileName, "w")

    for nx in range(0, flowField.Nx):
        for ny in range(0, flowField.Ny):
            for nz in range(0, flowField.Nz):
                for nd in range(0, flowField.Nd):
                    tmp = flowField.velocityField[nd, nx, ny, nz]
                    line = str(nd) + '\t'
                    line+= str(nx) + '\t'
                    line+= str(ny) + '\t'
                    line+= str(nz) + '\t'
                    line+= format(tmp, '.16f')
                    file.write(line + "\n")

    file.close()

    print('\nSaved ASCII Py file')

    return 0


def write_Details(flowField, output_directory, fileName):

    fileName += "_Details.txt"
    file = open(output_directory + fileName, "w")

    output = "\nDetails of the velocity field generated:\n"
    output+= "\nNd:\t\t" + str(flowField.Nd)
    output+= "\nNx:\t\t" + str(flowField.Nx)
    output+= "\nNy:\t\t" + str(flowField.Ny)
    output+= "\nNz:\t\t" + str(flowField.Nz)
    output+= "\nLx:\t\t" + str(flowField.Lx)
    output+= "\nLz:\t\t" + str(flowField.Lz) + "\n"
    output+= "\nRe:\t\t" + str(flowField.Re)
    output+= "\nc:\t\t"  + str(flowField.c) + "\n"
    output+= "\nwp:\t\t" + str(flowField.wavepacket)
    output+= "\nkx:\t\t" + str(flowField.kx)
    output+= "\nkz:\t\t" + str(flowField.kz)
    output+= "\nchi:\t\t"  + str(flowField.chi)
    output+= "\ntheta:\t\t" + str(flowField.theta)
    output+= "\nchi_tilde:\t" + str(flowField.chi_tilde) + "\n"
    #output+= "\nMx:" + str(flowField.Mx_full)
    #output+= "\nMz:" + str(flowField.Mz_full) + "\n"
    output+= "\nbf:\t\t" + str(flowField.baseflow) + "\n"
    file.write(output)
    file.close()
    print("\nSaved details file")

    return 0

def write_DAT(ff, output_directory, fileName):
    fileName += '.dat'
    file = open(output_directory + fileName, "w")

    title   = 'TITLE= "Initial flow field at Re = ' + str(ff.Re) + '"\n'
    columns = 'VARIABLES = "X", "Y", "Z", "U", "V", "W"\n'
    zones   = 'ZONE I=' + str(int(ff.Nx)) + ', J=' + str(int(ff.Ny)) + ', K=' + str(int(ff.Nz)) + ', DATAPACKING=POINT\n'
#    zones   = 'ZONE I=' + str(int(ff.Nx)) + ', J=' + str(int(ff.Ny)) + ', DATAPACKING=POINT\n'

    file.write(title)
    file.write(columns)
    file.write(zones)

    for nx in range(0, ff.Nx):
        for ny in range(0, ff.Ny):
            for nz in range(0, ff.Nz):
                string  = format(ff.x[nx], '.2f') + ' '
                string += format(ff.y[ny], '.2f') + ' '
                string += format(ff.z[nz], '.2f') + ' '
                string += format(ff.velocityField[0, nx, ny, nz], '.16f') + ' '
                string += format(ff.velocityField[1, nx, ny, nz], '.16f') + ' '
                string += format(ff.velocityField[2, nx, ny, nz], '.16f')

                file.write(string+'\n')

    file.close()

    print('\nSaved DAT file')

    return 0


def write_FF(output_directory, fileName):
    
    fileName = output_directory + fileName
    
    command = "ascii2field -p false -ge " + fileName + ".geom " + fileName + ".asc " + fileName
    
    os.system(command)

    return 0


def write_GEOM(flowField, output_directory, fileName):
    fileName += '.geom'
    file = open(output_directory + fileName, "w")

    file.write(str( int(flowField.Nx) ) + '\t\t\t\t\t\t%Nx' + "\n")
    file.write(str( int(flowField.Ny) ) + '\t\t\t\t\t\t%Ny' + "\n")
    file.write(str( int(flowField.Nz) ) + '\t\t\t\t\t\t%Nz' + "\n")
    file.write(str( int(flowField.Nd) ) + '\t\t\t\t\t\t%Nd' + "\n")


    Lx = format( flowField.Lx, '.16f')
    Lz = format( flowField.Lz, '.16f')
    file.write(Lx + '\t\t%Lx' + "\n")
    file.write(Lz + '\t\t%Lz' + "\n")

    lx = flowField.Lx / (2. * np.pi)
    lz = flowField.Lz / (2. * np.pi)
    file.write(str( lx ) + '\t\t\t\t\t\t%lx=Lx/(2pi)' + "\n")
    file.write(str( lz ) + '\t\t\t\t\t\t%lz=Lz/(2pi)' + "\n")

    alpha = float((2.* np.pi) / flowField.Lx)
    file.write(str( alpha ) + '\t\t\t\t\t\t%alpha=2pi/Lx' + "\n")

    gamma = float((2.* np.pi) / flowField.Lz)
    file.write(str( gamma ) + '\t\t\t\t\t\t%gamma=2pi/Lz' + "\n")

    file.close()

    print('\nSaved GEOM file')

    return 0


def write_Symms_File(directory, fileName, N, symStrAry):
    """
    The FieldSymmetry uses ASCII input-output. The storage format is
    c sx sy sz ax az

    where
    
    (c sx sy sz ax az) [u, v, w](x, y, z)  →  c [sx u, sy v, sz w](sx x + ax Lx, sy y, sz z + az Lz)
    """

    file = open(directory + fileName, "w")

    header = '% ' + str(N) + '\n'
    file.write(header)
    for i in range(0, N):
        string = symStrAry[i] + '\n'
        file.write(string)

#    if N == 1:
#        # simply write the sigma file
#        file.write(symStrAry[0])
#        
#    elif N > 1:
#        # Change format of file
#        header = '% ' + str(N) + '\n'
#        file.write(header)
#        for i in range(0, N):
#            string = symStrAry[i] + '\n'
#            file.write(string)

    file.close()

    return 0


def write_Vel_Profile(vel_profile, output_directory, fileName):
    
    fileName += ".txt"
    
    file = open(output_directory + fileName, "w")
    
    for i in range(0, len(vel_profile)):
        vel = format(vel_profile[i], '.16f')
        vel += "\n"
        file.write(vel)

    file.close()

    return 0


def plot_Contour(output_directory, fileName, 
                 xAxis, yAxis, data, 
                 xCoordStr, yCoordStr, zCoordStr, 
                 Re, 
                 nxCoordStr, nyCoordStr, nzCoordStr,
                 xAxisVel, yAxisVel,
                 xAxisName, yAxisName, 
                 velName, VL_max, VL_min,
                 quiv):

    xticks = np.linspace(xAxis[0], xAxis[-1], 4)
    yticks = np.linspace(yAxis[0], yAxis[-1], 5)

    x, y = np.meshgrid(xAxis, yAxis)
    v_min = np.amin(data)
    v_max = np.amax(data)
    v = np.linspace(VL_min, VL_max, 21, endpoint=True)
    ticks_at = [VL_min, v_min, 0, v_max, VL_max]
#    print(ticks_at)
    origin = 'lower'

    my_dpi = 150
    figx = abs(xAxis[0] - xAxis[-1]) * 300.0
    figy = abs(yAxis[0] - yAxis[-1]) * 150.0
    fig = plt.figure(figsize=(figx/my_dpi, figy/my_dpi), dpi=my_dpi)
#    fig = plt.figure(dpi=my_dpi)
    CS = plt.contourf(x, 
                      y, 
                      data, 
                      v,
                      cmap=cm.seismic,
                      origin=origin)

    if quiv:
        plt.quiver(x, y,
           xAxisVel,
           yAxisVel,
           color='k'
           )
                   

#    if xAxisName == "x" and yAxisName == "y":
#        plt.quiver(x, y,
#                   xAxisVel,
#                   yAxisVel,
#                   color='k',
#                   scale=1.5
#                   )
#    else:
#        plt.quiver(x, y,
#                   xAxisVel,
#                   yAxisVel,
#                   color='k'
#                   )


    plt.xlabel(xAxisName)
    plt.ylabel(yAxisName)

    title = "Re" + str(Re) + "\nt = " + fileName[1:] + "\n[ " + velName + " , " + xCoordStr + " , " + yCoordStr + " , " + zCoordStr + " ]"
    plt.title(title)

#    cbar = fig.colorbar(CS)
    cbar = fig.colorbar(CS,ticks=ticks_at,format='%1.2g')#,fraction=0.046, pad=0.04)
    cbar.ax.set_ylabel(velName)

    plt.axes().set_aspect('equal')

    plt.xticks(xticks)#, fontsize = 15)
    plt.yticks(yticks)#, fontsize = 15)

    fileName = output_directory + fileName + "_" + velName + "_x" + nxCoordStr + "_y" + nyCoordStr + "_z" + nzCoordStr + ".png"
    plt.savefig(fileName, bbox_inches='tight', dpi=96)
    plt.close(fig)

    return 0


def plot_Convergence_DNS(data, T0, T1): # include T0 and T1 in the name of the file

    fileName = "convergence_DNS_"+str(T0)+"-"+str(T1)+".png"

    x = data['t']
    y = data['L2Norm(u)']

    ymax = max(y) * 1.01 # 1% above max
    ymin = min(y) * 0.99 # 1% below min
    plt.figure()
    plt.plot(x, y, 'b-')
    plt.xlabel('t')
    plt.ylabel('$||u||_{2}$')
    plt.ylim([ymin, ymax])
    plt.grid(True)
    plt.title("DNS Convergence")
    plt.savefig(fileName)

    return 0


def plot_Vel_Profile(output_directory, fileName, vel_profile, y,
                     xAxisName, yAxisName):

    my_dpi = 96
    plt.figure(figsize=(450/my_dpi, 800/my_dpi), dpi=my_dpi)
    plt.plot(vel_profile, y, 'k-')
    plt.xlabel(xAxisName)
    plt.ylabel(yAxisName)
    plt.grid(True)
    plt.savefig(output_directory + fileName + ".png", dpi=my_dpi)

    return 0


def calculate_Difference(ff1, ff2): # 4D velocity fields

    delta = ff1 - ff2

    return delta

def calculate_Mean(dns_data_directory, tmp_directory, T0, T1):

    #### Change into the tmp directory
    os.chdir(tmp_directory)

    print("Currently in: " + tmp_directory)
    
    # range = T0:T1
    TSteps = T1 - T0 + 1
    TRange = np.linspace(T0, T1, TSteps)

    # Read geometry file to make empty array to store cumulative data
    var = read_GEOM(dns_data_directory[:dns_data_directory.find("data")], "u0")
    mean_ff = np.zeros((var['Nd'],
                        var['Nx'],
                        var['Ny'],
                        var['Nz']),
                        dtype=np.float128)

    #### Loop DNS data directory
    # Loop through all files in the DNS data folder (avoiding couette.args)
    files = [fi for fi in os.listdir(dns_data_directory) if os.path.isfile(os.path.join(dns_data_directory,fi))]

    files = sorted(files)

    for k in files:
        k = str(k)
        
        if k[0] == 'u' and k[-3:] == '.ff':
            timeUnit = float(k[1:-3])

            bool_start_pt = isclose(T0, timeUnit) # for the end points
            bool_end_pt   = isclose(T1, timeUnit) # for the end points
            bool_in_range = False

            if timeUnit > T0 and timeUnit < T1:
                bool_in_range = True

            if bool_in_range or bool_start_pt or bool_end_pt:

                # if the file is within the range specified
                if timeUnit in TRange:

                    # convert ff into ascii in the temporary folder
                    # (remember we are in the temporary folder within the data directory)
                    command = "field2ascii -p ../" + k + " " + k[:-3]
                    os.system(command)

                    # read the ascii file and add it to the flowField vector
                    var = read_ASC_PP(tmp_directory, k[:-3])
                    mean_ff += var['ff']

                    # remove contents of temporary folder
                    command = "rm -rf *.asc *.geom"
                    os.system(command)



    # divide the cumulative by the number of steps to get mean
    mean_ff = (1.0 / TSteps) * mean_ff

    var['ff'] = mean_ff

    return var


def calculate_Vel_Profiles(ff):

    cumulative = np.zeros((ff.Nd, ff.Nx, ff.Ny))
    vel_profile = np.zeros((ff.Nd, ff.Ny))
    z_avgd_mean = cumulative
    
    # Average in spanwise direction
    for i in range(0, ff.Nd):
        for nz in range(0, ff.Nz):
            cumulative[i, :, :] += ff.velocityField[i, :, :, nz]

    z_avgd_mean = (1.0/ff.Nz) * cumulative
    cumulative = vel_profile
    
    # Average in streamwise direction
    for i in range(0, ff.Nd):
        for nx in range(0, ff.Nx):
            cumulative[i, :] += z_avgd_mean[i, nx, :]

    vel_profile = (1.0/ff.Nx) * cumulative

    var = {}
    var['u'] = vel_profile[0, :]
    var['v'] = vel_profile[1, :]
    var['w'] = vel_profile[2, :]
    var['U'] = np.sqrt(var['u']**2 + var['v']**2 + var['w']**2)
    var['y'] = np.linspace(-1.0, 1.0, ff.Ny)

    return var


def print_Start_Bar():
    print('___________________________________________________________________________________\n\n')
    return 0


def print_Convergence_DNS():
    print('    ____  _   _______')
    print('   / __ \/ | / / ___/')
    print('  / / / /  |/ /\__ \ ')
    print(' / /_/ / /|  /___/ / ')
    print('/_____/_/ |_//____/  ')
    print("     _______  ___ _  _____ _______ ____ ___  _______ ")
    print("    / __/ _ \/ _ \ |/ / -_) __/ _ `/ -_) _ \/ __/ -_)")
    print("    \__/\___/_//_/___/\__/_/  \_, /\__/_//_/\__/\__/ ")
    print("                             /___/                   ")
    return 0

def print_DNSHeader():
    print('    ____  _   _______')
    print('   / __ \/ | / / ___/')
    print('  / / / /  |/ /\__ \ ')
    print(' / /_/ / /|  /___/ / ')
    print('/_____/_/ |_//____/  ')
    print('                     ')
    return 0


def print_DNSSubHeader():
    print('\nRunning DNS')
    print('___________________________________________________________________________________')
    return 0


def print_EndMessage():
    print('\nAll finished!')
    print('___________________________________________________________________________________\n\n')
    return 0


def print_ResolventHeader():
    print('___________________________________________________________________________________')
    print('                                                                                   ')
    print('   ________                           __                                           ')
    print('  / ____/ /_  ____ _____  ____  ___  / /                                           ')
    print(' / /   / __ \/ __ `/ __ \/ __ \/ _ \/ /                                            ')
    print('/ /___/ / / / /_/ / / / / / / /  __/ /                                             ')
    print('\____/_/ /_/\__,_/_/ /_/_/ /_/\_____/                __                 __         ')
    print('                              / __ \___  _________  / /   _____  ____  / /_        ')
    print('                             / /_/ / _ \/ ___/ __ \/ / | / / _ \/ __ \/ __/        ')
    print('                            / _, _/  __(__  ) /_/ / /| |/ /  __/ / / / /_          ')
    print('                           /_/ |_|\___/____/\____/_/ |___/\___/_/ /_/\__/          ')
    print('                                                                                   ')
    print(' Muhammad Arslan Ahmed                                                             ')
    print(' University of Southampton\n                                                       ') 
    print('___________________________________________________________________________________')
    return 0

def print_ResolventSubHeader():
    print('\nMaking a velocity field using the resolvent formulation')
    print('___________________________________________________________________________________')
    return 0

def print_Time(timeToPrint):
    print("\n\n")
    print("Time taken: " + str(timeToPrint))
    print("\n'n")

    return 0

def error(str):
    # Print the error and then exit from the program entirely
    print('\n\n')
    print('!!!!====!!!!====!!!!====!!!!====!!!!====')
    print('ERROR: ', str)
    print('!!!!====!!!!====!!!!====!!!!====!!!!====')
    print('\n')
    print('  ______    ______   __       __  ________       ')
    print(' /      \  /      \ |  \     /  \|        \      ')
    print('|  $$$$$$\|  $$$$$$\| $$\   /  $$| $$$$$$$$      ')
    print('| $$ __\$$| $$__| $$| $$$\ /  $$$| $$__          ')
    print('| $$|    \| $$    $$| $$$$\  $$$$| $$  \         ')
    print('| $$ \$$$$| $$$$$$$$| $$\$$ $$ $$| $$$$$         ')
    print('| $$__| $$| $$  | $$| $$ \$$$| $$| $$_____       ')
    print(' \$$    $$| $$  | $$| $$  \$ | $$| $$     \      ')
    print('  \$$$$$$  \$$   \$$ \$$      \$$ \$$$$$$$$      ')
    print('                                                 ')
    print('  ______   __     __  ________  _______          ')
    print(' /      \ |  \   |  \|        \|       \         ')
    print('|  $$$$$$\| $$   | $$| $$$$$$$$| $$$$$$$\        ')
    print('| $$  | $$| $$   | $$| $$__    | $$__| $$        ')
    print('| $$  | $$ \$$\ /  $$| $$  \   | $$    $$        ')
    print('| $$  | $$  \$$\  $$ | $$$$$   | $$$$$$$\        ')
    print('| $$__/ $$   \$$ $$  | $$_____ | $$  | $$        ')
    print(' \$$    $$    \$$$   | $$     \| $$  | $$        ')
    print('  \$$$$$$      \$     \$$$$$$$$ \$$   \$$        ')
    print('\n')
    print('!!!!====!!!!====!!!!====!!!!====!!!!====')
    print('ERROR: ', str)
    print('!!!!====!!!!====!!!!====!!!!====!!!!====')
    print('\n\n')
    sys.exit('')
    return
