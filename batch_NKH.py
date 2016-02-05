#!/usr/bin/env python

import os
print("")
print("Running Newton searches at Re 1800")
# RANK 10
print("\nRank 10")

d_1 = "/home/arslan/Documents/work/channelflow-related/set01/Re1800/KB/ampls-DNS-2015_11_20/wavepacket_001/data-skew/u600.000_rank_10"
os.chdir(d_1)
command_1 = "findsoln -orb -R 1800 -T 24.0 -sigma ../../sigma.asc -symms ../../symms.asc -nl skew -mc bulkv -Ubulk 1.33333333 -Uwall 0 -vdt -CFLmin 0.01 -CFLmax 1.0 -sn -log nkh-log.txt u600.000_rank_10.ff"
os.system(command_1)


# RANK 50
print("\nRank 50")
d_1 = "/home/arslan/Documents/work/channelflow-related/set01/Re1800/KB/ampls-DNS-2015_11_20/wavepacket_001/data-skew/u600.000_rank_50"
os.chdir(d_1)
command_1 = "findsoln -orb -R 1800 -T 24.0 -sigma ../../sigma.asc -symms ../../symms.asc -nl skew -mc bulkv -Ubulk 1.33333333 -Uwall 0 -vdt -CFLmin 0.01 -CFLmax 1.0 -sn -log nkh-log.txt u600.000_rank_50.ff"
os.system(command_1)


# RANK 100
print("\nRank 100")
d_1 = "/home/arslan/Documents/work/channelflow-related/set01/Re1800/KB/ampls-DNS-2015_11_20/wavepacket_001/data-skew/u600.000_rank_100"
os.chdir(d_1)
command_1 = "findsoln -orb -R 1800 -T 24.0 -sigma ../../sigma.asc -symms ../../symms.asc -nl skew -mc bulkv -Ubulk 1.33333333 -Uwall 0 -vdt -CFLmin 0.01 -CFLmax 1.0 -sn -log nkh-log.txt u600.000_rank_100.ff"
os.system(command_1)


# RANK 104
print("\nRank 104")
d_1 = "/home/arslan/Documents/work/channelflow-related/set01/Re1800/KB/ampls-DNS-2015_11_20/wavepacket_001/data-skew/u600.000_rank_104"
os.chdir(d_1)
command_1 = "findsoln -orb -R 1800 -T 24.0 -sigma ../../sigma.asc -symms ../../symms.asc -nl skew -mc bulkv -Ubulk 1.33333333 -Uwall 0 -vdt -CFLmin 0.01 -CFLmax 1.0 -sn -log nkh-log.txt u600.000_rank_104.ff"
os.system(command_1)