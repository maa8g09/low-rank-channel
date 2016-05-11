import numpy as np
import os
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
plt.title(r'ABC123 vs $\mathrm{ABC123}^{123}$')
import matplotlib.ticker as mtick

#### Keep all files in one directory just to make this plot...
directory = "/home/arslan/Documents/work/cfd-channelflow_solutions/w03_EQ/ALL"
os.chdir(directory)
ranks = [2,4,10,20,40,60, 'full']
equilibria = ["EQ1","EQ2","EQ3","EQ4","EQ5","EQ6","EQ7","EQ8","EQ9","EQ10","EQ11"]
norms_rank_02 = np.zeros((len(equilibria)))
norms_rank_04 = np.zeros((len(equilibria)))
norms_rank_10 = np.zeros((len(equilibria)))
norms_rank_20 = np.zeros((len(equilibria)))
norms_rank_40 = np.zeros((len(equilibria)))
norms_rank_60 = np.zeros((len(equilibria)))
norms_rank_full = np.zeros((len(equilibria)))
#### Loop through all files in the directory:
for root, sub_dirs, files in os.walk(directory):
    # for each file:
    files = sorted(files)
    for file in files:
        # Get norm
        command = "fieldprops -n " + str(file)
        output = os.popen(command).readlines()
        norm = float(output[2].split(" ")[-1][:-1])
        file_info = file.split("_")
        print(file_info)
        if len(file_info) > 1:
            index = equilibria.index(file_info[0].upper())
            rank = file_info[2].split(".")[0]
            if rank == "02":
                norms_rank_02[index] = norm
            elif rank == "04":
                norms_rank_04[index] = norm
            elif rank == "10":
                norms_rank_10[index] = norm
            elif rank == "20":
                norms_rank_20[index] = norm
            elif rank == "40":
                norms_rank_40[index] = norm
            elif rank == "60":
                norms_rank_60[index] = norm
            elif rank == "full":
                norms_rank_full[index] = norm
percentage_norm_rank_02   = norms_rank_02 / norms_rank_full * 100.0
percentage_norm_rank_04   = norms_rank_04 / norms_rank_full * 100.0
percentage_norm_rank_10   = norms_rank_10 / norms_rank_full * 100.0
percentage_norm_rank_20   = norms_rank_20 / norms_rank_full * 100.0
percentage_norm_rank_40   = norms_rank_40 / norms_rank_full * 100.0
percentage_norm_rank_60   = norms_rank_60 / norms_rank_full * 100.0
percentage_norm_rank_full = norms_rank_full / norms_rank_full * 100.0
stacked_norm_rank_02    = percentage_norm_rank_02
stacked_norm_rank_04    = percentage_norm_rank_04 - percentage_norm_rank_02
stacked_norm_rank_10    = percentage_norm_rank_10 - percentage_norm_rank_04
stacked_norm_rank_20    = percentage_norm_rank_20 - percentage_norm_rank_10
stacked_norm_rank_40    = percentage_norm_rank_40 - percentage_norm_rank_20
stacked_norm_rank_60    = percentage_norm_rank_60 - percentage_norm_rank_40
stacked_norm_rank_full  = percentage_norm_rank_full - percentage_norm_rank_60
#### Now we start plotting
fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(111)
colors ='cmybrgw'
equilibria = ("EQ1","EQ2","EQ3","EQ4","EQ5","EQ6","EQ7","EQ8","EQ9","EQ10","EQ11")
y_pos = np.arange(len(equilibria))
bar_width = 0.5
plt.axis([0,100,0,len(equilibria)])
plt.title("Norms")
index = np.arange(len(equilibria))
plt.barh(index,stacked_norm_rank_02,color=colors[0])
plt.barh(index,stacked_norm_rank_04,color=colors[1],left=stacked_norm_rank_02)
plt.barh(index,stacked_norm_rank_10,color=colors[2],left=stacked_norm_rank_02+stacked_norm_rank_04)
plt.barh(index,stacked_norm_rank_20,color=colors[3],left=stacked_norm_rank_02+stacked_norm_rank_04+stacked_norm_rank_10)
plt.barh(index,stacked_norm_rank_40,color=colors[4],left=stacked_norm_rank_02+stacked_norm_rank_04+stacked_norm_rank_10+stacked_norm_rank_20)
plt.barh(index,stacked_norm_rank_60,color=colors[5],left=stacked_norm_rank_02+stacked_norm_rank_04+stacked_norm_rank_10+stacked_norm_rank_20+stacked_norm_rank_40)
plt.barh(index,stacked_norm_rank_full,color=colors[6],left=stacked_norm_rank_02+stacked_norm_rank_04+stacked_norm_rank_10+stacked_norm_rank_20+stacked_norm_rank_40+stacked_norm_rank_60)
plt.yticks(index+0.5,equilibria)
ax.invert_yaxis()
fmt = "%.0f%%"
xticks = mtick.FFormatStrFormatter(fmt)
ax.xaxis.set_major_formatter(xticks)
ax.invert_yaxis()
plt.savefig(directory + "/plot.png", bbox_inches='tight')
plt.close(fig)
print("Finished")