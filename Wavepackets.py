import numpy as np

#=================== K1
kx1a = 6.0
kx1b = 6.0

kz1a = 6.0
kz1b = -6.0

a1a = 1.0j
a1b = 1.0j



#=================== K2
kx2a = 1.0
kx2b = 1.0

kz2a = 6.0
kz2b = -6.0

a2a = -4.5
a2b = -4.5



#=================== K3
kx3a = 7.0
kx3b = 7.0

kz3a = 12.0
kz3b = -12.0

a3a = 0.83j
a3b = 0.83j




#=================== K1-STAR
kx1a_star = 6.0
kx1b_star = 6.0

kz1a_star = 8.0
kz1b_star = -8.0

a1a_star = 1.0j
a1b_star = 1.0j



#=================== K3-STAR
kx3a_star = 7.0
kx3b_star = 7.0

kz3a_star = 14.0
kz3b_star = -14.0

a3a_star = 0.83j
a3b_star = 0.83j











#=================== K4
kx4a = 0.3
kx4b = 0.3

kz4a = 3.0
kz4b = -3.0

a4a = 0.3
a4b = 0.3

#=================== K5
kx5a = 1.5
kx5b = 1.5

kz5a = 4.0
kz5b = -4.0

a5a = 1.0
a5b = 1.0

#=================== K6
kx6a = 2.1
kx6b = 2.1

kz6a = 5.0
kz6b = -5.0

a6a = 3.0
a6b = 3.0

#=================== K7
kx7a = 1.0
kx7b = 1.0

kz7a = 6.0
kz7b = -6.0

a7a = 2.0
a7b = 2.0





wavepackets={}
wavepackets['KA_x'] = np.array([kx2a, kx2b])
wavepackets['KA_z'] = np.array([kz2a, kz2b])
wavepackets['KA_a'] = np.array([a2a, a2b])



wavepackets['KB_x'] = np.array([kx2a, kx2b, kx1a, kx1b])
wavepackets['KB_z'] = np.array([kz2a, kz2b, kz1a, kz1b])
wavepackets['KB_a'] = np.array([a2a, a2b, a1a, a1b])



wavepackets['KC_x'] = np.array([kx2a, kx2b, kx1a, kx1b, kx3a, kx3b])
wavepackets['KC_z'] = np.array([kz2a, kz2b, kz1a, kz1b, kz3a, kz3b])
wavepackets['KC_a'] = np.array([a2a, a2b, a1a, a1b, a3a, a3b])



wavepackets['KD_x'] = np.array([kx2a, kx2b, kx1a_star, kx1b_star, kx3a_star, kx3b_star])
wavepackets['KD_z'] = np.array([kz2a, kz2b, kz1a_star, kz1b_star, kz3a_star, kz3b_star])
wavepackets['KD_a'] = np.array([a2a, a2b, a1a_star, a1b_star, a3a_star, a3b_star])



wavepackets['KE_x'] = np.array([kx4a, kx4b, kx5a, kx5b, kx6a, kx6b, kx7a, kx7b])
wavepackets['KE_z'] = np.array([kz4a, kz4b, kz5a, kz5b, kz6a, kz6b, kz7a, kz7b])
wavepackets['KE_a'] = np.array([a4a, a4b, a5a, a5b, a6a, a6b, a7a, a7b])


## KA
#kx = np.array([kx1a, kx1b])
#kz = np.array([kz1a, kz1b])
#amplitudes = np.array([a1a, a1b])
#
#kx = np.array([kx2a, kx2b])
#kz = np.array([kz2a, kz2b])
#amplitudes = np.array([a2a, a2b])
#d = '/home/arslan/Documents/work/channelflow-related/set01/Re00400/kA/ampls01'






## KB
#kx = np.array([kx2a, kx2b, kx1a, kx1b])
#kz = np.array([kz2a, kz2b, kz1a, kz1b])
#amplitudes = np.array([a2a, a2b, a1a, a1b])
#d = '/home/arslan/Documents/work/channelflow-related/set01/Re01800/kB/ampls01'
##d = '/home/arslan/Desktop/asciibinarytests'





## KC
#kx = np.array([kx2a, kx2b, kx1a, kx1b, kx3a, kx3b])
#kz = np.array([kz2a, kz2b, kz1a, kz1b, kz3a, kz3b])
#amplitudes = np.array([a2a, a2b, a1a, a1b, a3a, a3b])
#d = '/home/arslan/Documents/work/channelflow-related/set01/Re01800/kC/ampls01'





## KD
#kx = np.array([kx2a, kx2b, kx1a_star, kx1b_star, kx3a_star, kx3b_star])
#kz = np.array([kz2a, kz2b, kz1a_star, kz1b_star, kz3a_star, kz3b_star])
#amplitudes = np.array([a2a, a2b, a1a_star, a1b_star, a3a_star, a3b_star])
#d = '/home/arslan/Documents/work/channelflow-related/set01/Re01800/kD/ampls01'





## KE
#kx = np.array([kx4a, kx4b, kx5a, kx5b, kx6a, kx6b, kx7a, kx7b])
#kz = np.array([kz4a, kz4b, kz5a, kz5b, kz6a, kz6b, kz7a, kz7b])
#amplitudes = np.array([a4a, a4b, a5a, a5b, a6a, a6b, a7a, a7b])
#d = '/home/arslan/Documents/work/channelflow-related/set01/Re01800/kE/ampls01'
