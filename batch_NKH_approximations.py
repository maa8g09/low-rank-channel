#!/usr/bin/env python3
import os
#### Eq1
name="eq1"
directory="/home/arslan/Documents/work/cfd-channelflow_solutions/w03_EQ/"+name+"/"
os.chdir(directory)
# Deconstruct
deconstruct = "eDeconstructField.py -d eq1_Details.txt -f "+name+".h5"
os.system(deconstruct)
# Construct
# Rank: Full
rankfull = "eConstructField.py -f "+name+"_deconstructed.h5 -r 1000"
os.system(rankfull)
fileName = name+"_rank_full"
rankfull_d = directory + fileName
os.chdir(rankfull_d)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)
# Rank: 2
r=2
os.chdir(directory)
rank = "eConstructField.py -f "+name+"_deconstructed.h5 -r " + str(r)
os.system(rank)
fileName = name+"_rank_" + str(r)
rank_dir = directory + fileName
os.chdir(rank_dir)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)
# Rank: 4
r=4
os.chdir(directory)
rank = "eConstructField.py -f "+name+"_deconstructed.h5 -r " + str(r)
os.system(rank)
fileName = name+"_rank_" + str(r)
rank_dir = directory + fileName
os.chdir(rank_dir)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)
# Rank: 10
r=10
os.chdir(directory)
rank = "eConstructField.py -f "+name+"_deconstructed.h5 -r " + str(r)
os.system(rank)
fileName = name+"_rank_" + str(r)
rank_dir = directory + fileName
os.chdir(rank_dir)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)
# Rank: 20
r=20
os.chdir(directory)
rank = "eConstructField.py -f "+name+"_deconstructed.h5 -r " + str(r)
os.system(rank)
fileName = name+"_rank_" + str(r)
rank_dir = directory + fileName
os.chdir(rank_dir)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)

#### Eq2
name="eq2"
directory="/home/arslan/Documents/work/cfd-channelflow_solutions/w03_EQ/"+name+"/"
os.chdir(directory)
# Deconstruct
deconstruct = "eDeconstructField.py -d eq1_Details.txt -f "+name+".h5"
os.system(deconstruct)
# Construct
# Rank: Full
rankfull = "eConstructField.py -f "+name+"_deconstructed.h5 -r 1000"
os.system(rankfull)
fileName = name+"_rank_full"
rankfull_d = directory + fileName
os.chdir(rankfull_d)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)
# Rank: 2
r=2
os.chdir(directory)
rank = "eConstructField.py -f "+name+"_deconstructed.h5 -r " + str(r)
os.system(rank)
fileName = name+"_rank_" + str(r)
rank_dir = directory + fileName
os.chdir(rank_dir)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)
# Rank: 4
r=4
os.chdir(directory)
rank = "eConstructField.py -f "+name+"_deconstructed.h5 -r " + str(r)
os.system(rank)
fileName = name+"_rank_" + str(r)
rank_dir = directory + fileName
os.chdir(rank_dir)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)
# Rank: 10
r=10
os.chdir(directory)
rank = "eConstructField.py -f "+name+"_deconstructed.h5 -r " + str(r)
os.system(rank)
fileName = name+"_rank_" + str(r)
rank_dir = directory + fileName
os.chdir(rank_dir)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)
# Rank: 20
r=20
os.chdir(directory)
rank = "eConstructField.py -f "+name+"_deconstructed.h5 -r " + str(r)
os.system(rank)
fileName = name+"_rank_" + str(r)
rank_dir = directory + fileName
os.chdir(rank_dir)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)

#### Eq3
name="eq3"
directory="/home/arslan/Documents/work/cfd-channelflow_solutions/w03_EQ/"+name+"/"
os.chdir(directory)
# Deconstruct
deconstruct = "eDeconstructField.py -d eq1_Details.txt -f "+name+".h5"
os.system(deconstruct)
# Construct
# Rank: Full
rankfull = "eConstructField.py -f "+name+"_deconstructed.h5 -r 1000"
os.system(rankfull)
fileName = name+"_rank_full"
rankfull_d = directory + fileName
os.chdir(rankfull_d)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)
# Rank: 2
r=2
os.chdir(directory)
rank = "eConstructField.py -f "+name+"_deconstructed.h5 -r " + str(r)
os.system(rank)
fileName = name+"_rank_" + str(r)
rank_dir = directory + fileName
os.chdir(rank_dir)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)
# Rank: 4
r=4
os.chdir(directory)
rank = "eConstructField.py -f "+name+"_deconstructed.h5 -r " + str(r)
os.system(rank)
fileName = name+"_rank_" + str(r)
rank_dir = directory + fileName
os.chdir(rank_dir)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)
# Rank: 10
r=10
os.chdir(directory)
rank = "eConstructField.py -f "+name+"_deconstructed.h5 -r " + str(r)
os.system(rank)
fileName = name+"_rank_" + str(r)
rank_dir = directory + fileName
os.chdir(rank_dir)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)
# Rank: 20
r=20
os.chdir(directory)
rank = "eConstructField.py -f "+name+"_deconstructed.h5 -r " + str(r)
os.system(rank)
fileName = name+"_rank_" + str(r)
rank_dir = directory + fileName
os.chdir(rank_dir)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)

#### Eq4
name="eq4"
directory="/home/arslan/Documents/work/cfd-channelflow_solutions/w03_EQ/"+name+"/"
os.chdir(directory)
# Deconstruct
deconstruct = "eDeconstructField.py -d eq1_Details.txt -f "+name+".h5"
os.system(deconstruct)
# Construct
# Rank: Full
rankfull = "eConstructField.py -f "+name+"_deconstructed.h5 -r 1000"
os.system(rankfull)
fileName = name+"_rank_full"
rankfull_d = directory + fileName
os.chdir(rankfull_d)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)
# Rank: 2
r=2
os.chdir(directory)
rank = "eConstructField.py -f "+name+"_deconstructed.h5 -r " + str(r)
os.system(rank)
fileName = name+"_rank_" + str(r)
rank_dir = directory + fileName
os.chdir(rank_dir)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)
# Rank: 4
r=4
os.chdir(directory)
rank = "eConstructField.py -f "+name+"_deconstructed.h5 -r " + str(r)
os.system(rank)
fileName = name+"_rank_" + str(r)
rank_dir = directory + fileName
os.chdir(rank_dir)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)
# Rank: 10
r=10
os.chdir(directory)
rank = "eConstructField.py -f "+name+"_deconstructed.h5 -r " + str(r)
os.system(rank)
fileName = name+"_rank_" + str(r)
rank_dir = directory + fileName
os.chdir(rank_dir)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)
# Rank: 20
r=20
os.chdir(directory)
rank = "eConstructField.py -f "+name+"_deconstructed.h5 -r " + str(r)
os.system(rank)
fileName = name+"_rank_" + str(r)
rank_dir = directory + fileName
os.chdir(rank_dir)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)


#### Eq5
name="eq5"
directory="/home/arslan/Documents/work/cfd-channelflow_solutions/w03_EQ/"+name+"/"
os.chdir(directory)
# Deconstruct
deconstruct = "eDeconstructField.py -d eq1_Details.txt -f "+name+".h5"
os.system(deconstruct)
# Construct
# Rank: Full
rankfull = "eConstructField.py -f "+name+"_deconstructed.h5 -r 1000"
os.system(rankfull)
fileName = name+"_rank_full"
rankfull_d = directory + fileName
os.chdir(rankfull_d)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)
# Rank: 2
r=2
os.chdir(directory)
rank = "eConstructField.py -f "+name+"_deconstructed.h5 -r " + str(r)
os.system(rank)
fileName = name+"_rank_" + str(r)
rank_dir = directory + fileName
os.chdir(rank_dir)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)
# Rank: 4
r=4
os.chdir(directory)
rank = "eConstructField.py -f "+name+"_deconstructed.h5 -r " + str(r)
os.system(rank)
fileName = name+"_rank_" + str(r)
rank_dir = directory + fileName
os.chdir(rank_dir)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)
# Rank: 10
r=10
os.chdir(directory)
rank = "eConstructField.py -f "+name+"_deconstructed.h5 -r " + str(r)
os.system(rank)
fileName = name+"_rank_" + str(r)
rank_dir = directory + fileName
os.chdir(rank_dir)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)
# Rank: 20
r=20
os.chdir(directory)
rank = "eConstructField.py -f "+name+"_deconstructed.h5 -r " + str(r)
os.system(rank)
fileName = name+"_rank_" + str(r)
rank_dir = directory + fileName
os.chdir(rank_dir)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)

#### Eq6
name="eq6"
directory="/home/arslan/Documents/work/cfd-channelflow_solutions/w03_EQ/"+name+"/"
os.chdir(directory)
# Deconstruct
deconstruct = "eDeconstructField.py -d eq1_Details.txt -f "+name+".h5"
os.system(deconstruct)
# Construct
# Rank: Full
rankfull = "eConstructField.py -f "+name+"_deconstructed.h5 -r 1000"
os.system(rankfull)
fileName = name+"_rank_full"
rankfull_d = directory + fileName
os.chdir(rankfull_d)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)
# Rank: 2
r=2
os.chdir(directory)
rank = "eConstructField.py -f "+name+"_deconstructed.h5 -r " + str(r)
os.system(rank)
fileName = name+"_rank_" + str(r)
rank_dir = directory + fileName
os.chdir(rank_dir)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)
# Rank: 4
r=4
os.chdir(directory)
rank = "eConstructField.py -f "+name+"_deconstructed.h5 -r " + str(r)
os.system(rank)
fileName = name+"_rank_" + str(r)
rank_dir = directory + fileName
os.chdir(rank_dir)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)
# Rank: 10
r=10
os.chdir(directory)
rank = "eConstructField.py -f "+name+"_deconstructed.h5 -r " + str(r)
os.system(rank)
fileName = name+"_rank_" + str(r)
rank_dir = directory + fileName
os.chdir(rank_dir)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)
# Rank: 20
r=20
os.chdir(directory)
rank = "eConstructField.py -f "+name+"_deconstructed.h5 -r " + str(r)
os.system(rank)
fileName = name+"_rank_" + str(r)
rank_dir = directory + fileName
os.chdir(rank_dir)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)

#### Eq7
name="eq7"
directory="/home/arslan/Documents/work/cfd-channelflow_solutions/w03_EQ/"+name+"/"
os.chdir(directory)
# Deconstruct
deconstruct = "eDeconstructField.py -d eq1_Details.txt -f "+name+".h5"
os.system(deconstruct)
# Construct
# Rank: Full
rankfull = "eConstructField.py -f "+name+"_deconstructed.h5 -r 1000"
os.system(rankfull)
fileName = name+"_rank_full"
rankfull_d = directory + fileName
os.chdir(rankfull_d)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)
# Rank: 2
r=2
os.chdir(directory)
rank = "eConstructField.py -f "+name+"_deconstructed.h5 -r " + str(r)
os.system(rank)
fileName = name+"_rank_" + str(r)
rank_dir = directory + fileName
os.chdir(rank_dir)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)
# Rank: 4
r=4
os.chdir(directory)
rank = "eConstructField.py -f "+name+"_deconstructed.h5 -r " + str(r)
os.system(rank)
fileName = name+"_rank_" + str(r)
rank_dir = directory + fileName
os.chdir(rank_dir)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)
# Rank: 10
r=10
os.chdir(directory)
rank = "eConstructField.py -f "+name+"_deconstructed.h5 -r " + str(r)
os.system(rank)
fileName = name+"_rank_" + str(r)
rank_dir = directory + fileName
os.chdir(rank_dir)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)
# Rank: 20
r=20
os.chdir(directory)
rank = "eConstructField.py -f "+name+"_deconstructed.h5 -r " + str(r)
os.system(rank)
fileName = name+"_rank_" + str(r)
rank_dir = directory + fileName
os.chdir(rank_dir)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)


#### Eq8
name="eq8"
directory="/home/arslan/Documents/work/cfd-channelflow_solutions/w03_EQ/"+name+"/"
os.chdir(directory)
# Deconstruct
deconstruct = "eDeconstructField.py -d eq1_Details.txt -f "+name+".h5"
os.system(deconstruct)
# Construct
# Rank: Full
rankfull = "eConstructField.py -f "+name+"_deconstructed.h5 -r 1000"
os.system(rankfull)
fileName = name+"_rank_full"
rankfull_d = directory + fileName
os.chdir(rankfull_d)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)
# Rank: 2
r=2
os.chdir(directory)
rank = "eConstructField.py -f "+name+"_deconstructed.h5 -r " + str(r)
os.system(rank)
fileName = name+"_rank_" + str(r)
rank_dir = directory + fileName
os.chdir(rank_dir)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)
# Rank: 4
r=4
os.chdir(directory)
rank = "eConstructField.py -f "+name+"_deconstructed.h5 -r " + str(r)
os.system(rank)
fileName = name+"_rank_" + str(r)
rank_dir = directory + fileName
os.chdir(rank_dir)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)
# Rank: 10
r=10
os.chdir(directory)
rank = "eConstructField.py -f "+name+"_deconstructed.h5 -r " + str(r)
os.system(rank)
fileName = name+"_rank_" + str(r)
rank_dir = directory + fileName
os.chdir(rank_dir)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)
# Rank: 20
r=20
os.chdir(directory)
rank = "eConstructField.py -f "+name+"_deconstructed.h5 -r " + str(r)
os.system(rank)
fileName = name+"_rank_" + str(r)
rank_dir = directory + fileName
os.chdir(rank_dir)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)


#### Eq9
name="eq9"
directory="/home/arslan/Documents/work/cfd-channelflow_solutions/w03_EQ/"+name+"/"
os.chdir(directory)
# Deconstruct
deconstruct = "eDeconstructField.py -d eq1_Details.txt -f "+name+".h5"
os.system(deconstruct)
# Construct
# Rank: Full
rankfull = "eConstructField.py -f "+name+"_deconstructed.h5 -r 1000"
os.system(rankfull)
fileName = name+"_rank_full"
rankfull_d = directory + fileName
os.chdir(rankfull_d)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)
# Rank: 2
r=2
os.chdir(directory)
rank = "eConstructField.py -f "+name+"_deconstructed.h5 -r " + str(r)
os.system(rank)
fileName = name+"_rank_" + str(r)
rank_dir = directory + fileName
os.chdir(rank_dir)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)
# Rank: 4
r=4
os.chdir(directory)
rank = "eConstructField.py -f "+name+"_deconstructed.h5 -r " + str(r)
os.system(rank)
fileName = name+"_rank_" + str(r)
rank_dir = directory + fileName
os.chdir(rank_dir)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)
# Rank: 10
r=10
os.chdir(directory)
rank = "eConstructField.py -f "+name+"_deconstructed.h5 -r " + str(r)
os.system(rank)
fileName = name+"_rank_" + str(r)
rank_dir = directory + fileName
os.chdir(rank_dir)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)
# Rank: 20
r=20
os.chdir(directory)
rank = "eConstructField.py -f "+name+"_deconstructed.h5 -r " + str(r)
os.system(rank)
fileName = name+"_rank_" + str(r)
rank_dir = directory + fileName
os.chdir(rank_dir)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)



#### Eq10
name="eq10"
directory="/home/arslan/Documents/work/cfd-channelflow_solutions/w03_EQ/"+name+"/"
os.chdir(directory)
# Deconstruct
deconstruct = "eDeconstructField.py -d eq1_Details.txt -f "+name+".h5"
os.system(deconstruct)
# Construct
# Rank: Full
rankfull = "eConstructField.py -f "+name+"_deconstructed.h5 -r 1000"
os.system(rankfull)
fileName = name+"_rank_full"
rankfull_d = directory + fileName
os.chdir(rankfull_d)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)
# Rank: 2
r=2
os.chdir(directory)
rank = "eConstructField.py -f "+name+"_deconstructed.h5 -r " + str(r)
os.system(rank)
fileName = name+"_rank_" + str(r)
rank_dir = directory + fileName
os.chdir(rank_dir)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)
# Rank: 4
r=4
os.chdir(directory)
rank = "eConstructField.py -f "+name+"_deconstructed.h5 -r " + str(r)
os.system(rank)
fileName = name+"_rank_" + str(r)
rank_dir = directory + fileName
os.chdir(rank_dir)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)
# Rank: 10
r=10
os.chdir(directory)
rank = "eConstructField.py -f "+name+"_deconstructed.h5 -r " + str(r)
os.system(rank)
fileName = name+"_rank_" + str(r)
rank_dir = directory + fileName
os.chdir(rank_dir)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)
# Rank: 20
r=20
os.chdir(directory)
rank = "eConstructField.py -f "+name+"_deconstructed.h5 -r " + str(r)
os.system(rank)
fileName = name+"_rank_" + str(r)
rank_dir = directory + fileName
os.chdir(rank_dir)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)



#### Eq11
name="eq11"
directory="/home/arslan/Documents/work/cfd-channelflow_solutions/w03_EQ/"+name+"/"
os.chdir(directory)
# Deconstruct
deconstruct = "eDeconstructField.py -d eq1_Details.txt -f "+name+".h5"
os.system(deconstruct)
# Construct
# Rank: Full
rankfull = "eConstructField.py -f "+name+"_deconstructed.h5 -r 1000"
os.system(rankfull)
fileName = name+"_rank_full"
rankfull_d = directory + fileName
os.chdir(rankfull_d)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)
# Rank: 2
r=2
os.chdir(directory)
rank = "eConstructField.py -f "+name+"_deconstructed.h5 -r " + str(r)
os.system(rank)
fileName = name+"_rank_" + str(r)
rank_dir = directory + fileName
os.chdir(rank_dir)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)
# Rank: 4
r=4
os.chdir(directory)
rank = "eConstructField.py -f "+name+"_deconstructed.h5 -r " + str(r)
os.system(rank)
fileName = name+"_rank_" + str(r)
rank_dir = directory + fileName
os.chdir(rank_dir)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)
# Rank: 10
r=10
os.chdir(directory)
rank = "eConstructField.py -f "+name+"_deconstructed.h5 -r " + str(r)
os.system(rank)
fileName = name+"_rank_" + str(r)
rank_dir = directory + fileName
os.chdir(rank_dir)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)
# Rank: 20
r=20
os.chdir(directory)
rank = "eConstructField.py -f "+name+"_deconstructed.h5 -r " + str(r)
os.system(rank)
fileName = name+"_rank_" + str(r)
rank_dir = directory + fileName
os.chdir(rank_dir)
nkh = "findsoln -eqb -log nkh-"+fileName+".log "+fileName+".h5"
print(nkh)
os.system(nkh)

