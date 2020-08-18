"""Plot time-space diagram of Otago faults active periods
"""
#import matplotlib
from matplotlib import pyplot as plt
import matplotlib.patches as patches
import numpy as np

Akatore = [[-1,13]]
Ak_n = ['     3']
#Titri = [[17.5, 18.5], [71, 128]] # real data
Titri = [[17.5, 18.5], [71, 80]] # Dummy data for nice plotting
T_n = ['?', '> 1']
Hyde = [[11., 25.]]
H_n = ['     2']
Dunstan = [[15,24], [54,56]]
D_n = [' 4 - 6', '?']
Pisa_L = [[25,25]] # Dummy data for plotting text                                                                               
P_L_n = ['No late Quaternary deformation']    
Pisa_GV = [[35,70],[200,250]]
P_GV_n = ['> 20?', '?']
Cardrona_BK = [[0,18]]
C_BK_n = ['    3 - 4']
Cardrona_Kawarau = [[0,7],[10,11],[20,23]]                                                                                        
C_K_n = ['1', '1', '1']
Nevis_BN = [[18,70]]
N_BN_n = ['1 - 2']
Nevis_CC = [[0,11],[13, 18]]
N_CC_n = ['1', '1']
Little_Nevis = [[0, 18]]
LN_n = ['      4']
#faults = [Akatore, Titri, Hyde, Dunstan, Pisa_L, Pisa_GV, Cardrona_BK, Cardrona_Kawarau, Nevis_BN, Nevis_CC,
#          Little_Nevis]
#event_numbers = [Ak_n, T_n, H_n, D_n, P_L_n, P_GV_n, C_BK_n, C_K_n, N_BN_n, N_CC_n, LN_n]
# Exclude Hyde Fault for now
faults = [Akatore, Titri, Dunstan, Pisa_L, Pisa_GV, Cardrona_BK, Cardrona_Kawarau, Nevis_BN, Nevis_CC,
          Little_Nevis]
event_numbers = [Ak_n, T_n, D_n, P_L_n, P_GV_n, C_BK_n, C_K_n, N_BN_n, N_CC_n, LN_n]
#plt.xlim([-5, 50])
for i, fault in enumerate(faults):
    for j, date_range in enumerate(fault):
        plt.plot(date_range, [i,i], c='k', linewidth=1)
        plt.text(date_range[0], i+0.1, event_numbers[i][j])
# plot Akatore open interval
plt.plot([125,150],[0,0], c='k', linewidth=1, linestyle=':')
plt.text(124, 0.1, '? > 125 ka')
plt.ylim([-0.2, len(faults)-0.7])
#plt.xlim([-5, 50]) 
plt.yticks(np.arange(0, len(faults)+1, 1),
           ['Akatore', 'Titri', #'Hyde',
            'Dunstan', 'Pisa-\nLindis', 'Pisa-\nGrandview', 'NW Cardrona\n(Branch Ck)',
            'NW Cardrona\n(Kawarau Gorge)', 'Nevis\n(Ben Nevis)', 'Nevis\n(Coal Creek)',
            'Little\nNevis'])
plt.xlabel('Thousands of years before present')
#plt.arrow(250,0,0,4)
plt.annotate('', xy=(1.05, 0), xycoords='axes fraction', xytext=(1.05, 1), 
            arrowprops=dict(arrowstyle="<-", color='k'))
plt.text(0.95, 0.7, 'Proximity to Alpine Fault', transform=plt.gcf().transFigure,
         rotation=-90)
#plt.plot([80,80],[-1,len(faults)+1], linestyle=':', linewidth=0.2)   
plt.tight_layout()
plt.savefig('Otago_time_space.png')
        
plt.xlim([-5, 100])
plt.xticks(np.arange(0, 90, 10))
# Create a Rectangle patch
rect = patches.Rectangle((-5,5.5),85,5.5,linewidth=0.5,edgecolor='k',facecolor='none',
                         linestyle='-.')
# Add the patch to the Axes
ax = plt.gca()
ax.add_patch(rect)
#plt.text(30, 3.1, 'No late Quaternary deformation')
plt.text(81, 0.1, '? > 125 ka')
plt.text(81, 1.1, '? < 128 ka')
#plt.text(55, 2.1, '? > 55 ka')
plt.text(81, 5.1, '? ~ 250 ka')
plt.text(40,8.6, 'Nevis-Cardrona\nFault System')
plt.plot([80,80],[-1,len(faults)+1], linestyle=':', linewidth=1, c='0.5')
plt.savefig('Otago_time_space_80ka.png') 
         
