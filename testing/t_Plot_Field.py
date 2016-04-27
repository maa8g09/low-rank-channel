#!/usr/bin/env python3
import os
import numpy as np
import sys
sys.path.append("/".join(sys.path[0].split("/")[:-1]))
import Utils as ut
import FlowField as ffClass
import Tests

def main(File, Details, VelComponent, SpatialComponent, Coordinate,
         Full, SpatiallyAveraged, Quiver, Levels):
    #===================================================================#
    #### Format current directory path                               ####
    #===================================================================#
    pwd = ut.format_Directory_Path(os.getcwd())
    #===================================================================#
    #### Read file and details file                                  ####
    #===================================================================#
    file_info, original_attrs = ut.read_H5(File)
    details = ut.read_Details(Details)
    #===================================================================#
    #### Make output directory for images and change into it         ####
    #===================================================================#
    images_directory = ut.make_Folder(pwd, "images_" + str(File)[:-3], False)
    os.chdir(images_directory)
    #===================================================================#
    #### Initialise flow field class with given file                 ####
    #===================================================================#
    ff = ffClass.FlowFieldChannelFlow(  file_info['Nd'],
                                        file_info['Nx'],
                                        file_info['Ny'],
                                        file_info['Nz'],
                                        file_info['Lx'],
                                        file_info['Lz'],
                                        file_info['alpha'],
                                        file_info['beta'],
                                        details['c'],
                                        details['bf'],
                                        details['Re'],
                                        file_info['ff'],
                                        "pp")
    #===================================================================#
    #### Set velocity and spatial components                         ####
    #===================================================================#
    i = VelComponent
    n = SpatialComponent
    m=n+1; p=n-1
    if m==3:
        m=0
    if p==-1:
        p=2
    velName = ""
    if i==0:
        velName = "u"
    elif i==1:
        velName = "v"
    elif i==2:
        velName = "w"
    #===================================================================#
    #### Start plotting                                              ####
    #===================================================================#
    print("Plotting...")
    vl_max = np.amax(ff.velocityField[i, :, :, :])
    vl_min = np.amin(ff.velocityField[i, :, :, :])
    
    printFullTitle = True
    if Levels:
        contour_levels = Levels
    else:
        contour_levels = 20
    #===================================================================#
    #### If Full selected                                            ####
    #===================================================================#
    if Full:
        if n == 0:
            x_directory = ut.make_Folder(images_directory, "x", False)
            for j in range(0, ff.Nx):
                ut.plot_Contour(x_directory, str(File[:-3]), 
                                ff.z, ff.y, ff.velocityField[i, j, :, :], 
                                format(ff.x[j], '.2f'), ":", ":", 
                                ff.Re, 
                                str(j),  "-", "-",
                                ff.velocityField[m, j, :, :], 
                                ff.velocityField[p, j, :, :], 
                                "z", "y", velName,
                                vl_max, vl_min, Quiver, contour_levels,
                                printFullTitle)
        elif n == 1:
            y_directory = ut.make_Folder(images_directory, "y", False)
            for j in range(0, ff.Ny):
                ut.plot_Contour(y_directory, str(File[:-3]), 
                                ff.z, ff.x, ff.velocityField[i, :, j, :], 
                                ":", format(ff.y[j], '.2f'), ":", 
                                ff.Re, 
                                "-", str(j), "-", 
                                ff.velocityField[m, :, j, :], 
                                ff.velocityField[p, :, j, :], 
                                "z", "x", velName,
                                vl_max, vl_min, Quiver, contour_levels,
                                printFullTitle)
        elif n == 2:
            z_directory = ut.make_Folder(images_directory, "z", False)
            for j in range(0, ff.Nz):
                ut.plot_Contour(z_directory, str(File[:-3]), 
                                ff.x, ff.y, ff.velocityField[i, :, :, j].T, 
                                ":", ":", format(ff.z[j], '.2f'), 
                                ff.Re, 
                                "-", "-", str(j), 
                                ff.velocityField[m, :, :, j].T, 
                                ff.velocityField[p, :, :, j].T, 
                                "x", "y", velName,
                                vl_max, vl_min, Quiver, contour_levels,
                                printFullTitle)
    #===================================================================#
    #### If Co-Ordinate given                                        ####
    #===================================================================#
    printFullTitle = True
    if Coordinate:
        if n == 0:
            coords = Tests.indices(ff.x, lambda m: m > float(Coordinate))
            x_coord = coords[0]
            x_directory = ut.make_Folder(images_directory, "x", False)
            ut.plot_Contour(x_directory, str(File[:-3]), 
                            ff.z, ff.y, ff.velocityField[i, x_coord, :, :], 
                            format(ff.x[x_coord], '.2f'), ":", ":", 
                            ff.Re, 
                            str(x_coord),  "-", "-",
                            ff.velocityField[m, x_coord, :, :], 
                            ff.velocityField[p, x_coord, :, :], 
                            "z", "y", velName,
                            vl_max, vl_min, Quiver, contour_levels,
                            printFullTitle)
        elif n == 1:
            coords = Tests.indices(ff.y, lambda m: m > float(Coordinate))
            y_coord = coords[-1]
            y_directory = ut.make_Folder(images_directory, "y", False)
            ut.plot_Contour(y_directory, str(File[:-3]), 
                            ff.z, ff.x, ff.velocityField[i, :, y_coord, :], 
                            ":", format(ff.y[y_coord], '.2f'), ":", 
                            ff.Re, 
                            "-", str(y_coord), "-", 
                            ff.velocityField[m, :, y_coord, :], 
                            ff.velocityField[p, :, y_coord, :], 
                            "z", "x", velName,
                            vl_max, vl_min, Quiver, contour_levels,
                            printFullTitle)
        elif n == 2:
            coords = Tests.indices(ff.z, lambda m: m > float(Coordinate))
            z_coord = coords[0]
            z_directory = ut.make_Folder(images_directory, "z", False)
            ut.plot_Contour(z_directory, str(File[:-3]), 
                            ff.x, ff.y, ff.velocityField[i, :, :, z_coord].T, 
                            ":", ":", format(ff.z[z_coord], '.2f'), 
                            ff.Re, 
                            "-", "-", str(z_coord), 
                            ff.velocityField[m, :, :, z_coord].T, 
                            ff.velocityField[p, :, :, z_coord].T, 
                            "x", "y", velName,
                            vl_max, vl_min, Quiver, contour_levels,
                            printFullTitle)
    #===================================================================#
    #### If Spatially averaging selected                             ####
    #===================================================================#
    printFullTitle = False
    if SpatiallyAveraged:
        if n == 0:
            # Average the flow field in teh streamwise direction.
            x_averaged_ff = np.zeros((ff.Nd, ff.Ny, ff.Nz))
            for nd in range(0, ff.Nd):
                for nx in range(0, ff.Nx):
                    x_averaged_ff[nd, :, :] += ff.velocityField[nd, nx, :, :]
            x_averaged_ff *= (1.0 / ff.Nx)
            vl_max = np.amax(x_averaged_ff[i, :, :])
            vl_min = np.amin(x_averaged_ff[i, :, :])
            fileName = str(File)[:-3] + "_x_avgd_"
            ut.plot_Contour(images_directory, fileName, 
                            ff.z, ff.y, x_averaged_ff[i, :, :], 
                            "avg", ":", ":", 
                            ff.Re, 
                            "x-avg",  "-", "-",
                            x_averaged_ff[m, :, :], 
                            x_averaged_ff[p, :, :], 
                            "z", "y", velName,
                            vl_max, vl_min, Quiver, contour_levels,
                            printFullTitle)
        elif n == 1:
            # Average the flow field in teh streamwise direction.
            y_averaged_ff = np.zeros((ff.Nd, ff.Nx, ff.Nz))
            for nd in range(0, ff.Nd):
                for ny in range(0, ff.Ny):
                    y_averaged_ff[nd, :, :] += ff.velocityField[nd, :, ny, :]
            y_averaged_ff *= (1.0 / ff.Ny)
            vl_max = np.amax(y_averaged_ff[i, :, :])
            vl_min = np.amin(y_averaged_ff[i, :, :])
            fileName = str(File)[:-3] + "_y_avgd_"
            ut.plot_Contour(images_directory, fileName, 
                            ff.z, ff.x, y_averaged_ff[i, :, :], 
                            ":", "avg", ":", 
                            ff.Re, 
                            "-", "y-avg", "-", 
                            y_averaged_ff[m, :, :], 
                            y_averaged_ff[p, :, :], 
                            "z", "x", velName,
                            vl_max, vl_min, Quiver, contour_levels,
                            printFullTitle)
        elif n == 2:
            # Average the flow field in teh streamwise direction.
            z_averaged_ff = np.zeros((ff.Nd, ff.Nx, ff.Ny))
            for nd in range(0, ff.Nd):
                for nz in range(0, ff.Nz):
                    z_averaged_ff[nd, :, :] += ff.velocityField[nd, :, :, nz]
            z_averaged_ff *= (1.0 / ff.Nz)
            vl_max = np.amax(z_averaged_ff[i, :, :])
            vl_min = np.amin(z_averaged_ff[i, :, :])
            fileName = str(File)[:-3] + "_z_avgd_"
            ut.plot_Contour(images_directory, fileName,
                            ff.x, ff.y, z_averaged_ff[i, :, :].T, 
                            ":", ":", "avg", 
                            ff.Re, 
                            "-", "-", "z-avg", 
                            z_averaged_ff[m, :, :].T, 
                            z_averaged_ff[p, :, :].T, 
                            "x", "y", velName,
                            vl_max, vl_min, Quiver, contour_levels,
                            printFullTitle)
    ut.print_EndMessage()


directory = "/home/arslan/Documents/work/cfd-channelflow_solutions/w03_TWs/tw2"
os.chdir(directory)

file             = "tw2.h5"
details          = "tw2_Details.txt"
velComponent     = 0
spatialComponent = 0
coordinate       = 2.5
full             = False
spatiallyAveraged= False
quiver           = False
levels           = 50

main(file,
     details,
     velComponent,
     spatialComponent,
     coordinate,
     full,
     spatiallyAveraged,
     quiver,
     levels)
