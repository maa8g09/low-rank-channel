import os
import sys
sys.path.append("/".join(sys.path[0].split("/")[:-1]))
import Utils as ut

def main(File, Directory):
    # If a convergence file given:
    # plot it
    
    # If a directory given:
    # plot each convergence file in each subdirectory
    # 
    # check to see if a convergence file exists in the directory given, if there is... use it as a constant line.
    
    if Directory: # if a directory is given...
        directory = ut.format_Directory_Path(Directory)
        os.chdir(directory)
        # the last directory name is the case name.
        fileName = directory.split("/")
        fileName = fileName[-2]
        # Read convergence files in this directory and all sub directories...
        all_convergence_data = {}
        for root, sub_dirs, files in os.walk(directory):
            for file in files:
                if file.endswith("convergence.asc"):
                     print(os.path.join(root, file))
                     # read in convergence file
                     # initialise convergence dictionary
                     all_convergence_data[fileName] = {}
                     # populate the dictionary
                     all_convergence_data[fileName] = ut.read_Convergence_NKH(file)
            print("")        
            for sub_dir in sub_dirs:
                sub_fileName = sub_dir
                sub_dir = os.path.join(root, sub_dir)
                os.chdir(sub_dir)
                for file in os.listdir(sub_dir):
                    if file.endswith("convergence.asc"):
                        print(os.path.join(sub_dir, file))
                        # read in convergence file
                        # initialise convergence dictionary
                        all_convergence_data[sub_fileName] = {}
                        # populate the dictionary
                        all_convergence_data[sub_fileName] = ut.read_Convergence_NKH(file)
            break
        os.chdir(directory)
        ut.plot_Convergence_NKH_multi(all_convergence_data, "Newton_Steps", "L2Norm(u)")
        ut.plot_Convergence_NKH_multi(all_convergence_data, "ftotal", "L2Norm(u)")
        ut.plot_Convergence_NKH_multi(all_convergence_data, "fnewt", "L2Norm(u)")
        ut.plot_Convergence_NKH_multi(all_convergence_data, "fhook", "L2Norm(u)")
        
        ut.plot_Convergence_NKH_multi(all_convergence_data, "Newton_Steps", "L2Norm(G)")
        ut.plot_Convergence_NKH_multi(all_convergence_data, "ftotal", "L2Norm(G)")
        ut.plot_Convergence_NKH_multi(all_convergence_data, "fnewt", "L2Norm(G)")
        ut.plot_Convergence_NKH_multi(all_convergence_data, "fhook", "L2Norm(G)")
        
        ut.plot_Convergence_NKH_multi(all_convergence_data, "Newton_Steps", "L2Norm(du)")
        ut.plot_Convergence_NKH_multi(all_convergence_data, "ftotal", "L2Norm(du)")
        ut.plot_Convergence_NKH_multi(all_convergence_data, "fnewt", "L2Norm(du)")
        ut.plot_Convergence_NKH_multi(all_convergence_data, "fhook", "L2Norm(du)")
        
        ut.plot_Convergence_NKH_multi(all_convergence_data, "Newton_Steps", "GMRESerr")
        
    
di = "/Users/arslan/temp/eq1"
f = ""
main(f, di)
