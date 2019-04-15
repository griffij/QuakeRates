"""Plot long-term vs quiescent rates of earthquake occurrence
"""

import numpy as np
#import matplotlib
from matplotlib import pyplot as plt

akatore = [115000, 600, 4, 125000]
cadell = [1500000, 7500, 8, 4500000]
dunstan = [30000, 2000, 7, 55000]
lake_george = [750000, 3000, 30, 3000000]

wharekuri = [65000, 7500, 5, 150000]


galeen=[30000, 5000, 5, 99000]
roer = [95000, 600, 5, 101000]
alpine = [330, 330, 22, 8000]

data = np.array([akatore, cadell, dunstan, galeen, roer, wharekuri, alpine])

x = data[:,0]
y = data[:,1]
print(x)
print(y)

plt.scatter(x,y)
ax = plt.gca()
ax.set_xscale('log')
ax.set_xlabel('Quiescent recurrence (years)')
ax.set_ylabel('Active recurrence (years)')
              
plt.savefig('longterm_vs_quiescence.png')

x = data[:,3]                                                                                                                   
y = data[:,2]                                                                                                                    
print(x)                                                                                                                         
print(y)
plt.clf()
plt.scatter(x,y)                                                                                                                
ax = plt.gca()
ax.set_xscale('log')  
ax.set_xlabel('Time (years)')                                                                                    
ax.set_ylabel('Number of events')
plt.savefig('number_event_vs_time.png')  

x = data[:,0]                                                                                                                   
y = data[:,2]                                                                                                                   
print(x)                                                                                                                        
print(y)                                                                                                                        
plt.clf()                                                                                                                       
plt.scatter(x,y)                                                                                                                
ax = plt.gca()                                                                                                                  
ax.set_xlabel('Quiescent Recurrence (years)')                                                                                   
ax.set_ylabel('Number of events')
ax.set_xscale('log')  
plt.savefig('number_event_vs_quiescent_period.png')


x = data[:,1]                                                                                                                   
y = data[:,2]                                                                                                                   
print(x)                                                                                                                        
print(y)                                                                                                                        
plt.clf()                                                                                                                       
plt.scatter(x,y)                                                                                                                
ax = plt.gca()
ax.set_xscale('log')
ax.set_xlabel('Active Recurrence (years)')                                                                                   
ax.set_ylabel('Number of events')                                                                                               
plt.savefig('number_event_vs_active_period.png')   

x = data[:,0]
y = data[:,2]/data[:,3]
print(x)
print(y)
plt.clf()
plt.scatter(x,y)
ax = plt.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim([0.9*min(y), 1.1*max(y)])
ax.set_xlabel('Quiescent Recurrence (years)')
ax.set_ylabel('Events per year')
plt.tight_layout()
plt.savefig('event_per_year_vs_quiescent_period.png')

x = data[:,1]
y = data[:,2]/data[:,3]
print(x)
print(y)
plt.clf()
plt.scatter(x,y)
ax = plt.gca()
ax.set_xscale('log')
ax.set_yscale('log') 
ax.set_ylim([0.9*min(y), 1.1*max(y)])
ax.set_xlabel('Active Recurrence (years)')
ax.set_ylabel('Events per year')
plt.tight_layout()
plt.savefig('event_per_year_vs_active_period.png')
