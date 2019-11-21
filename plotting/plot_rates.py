"""Plot long-term vs quiescent rates of earthquake occurrence
"""

import numpy as np
from adjustText import adjust_text
#import matplotlib
from matplotlib import pyplot as plt

akatore = [115000, 600, 4, 125000, 'Akatore']
cadell = [1500000, 7500, 8, 4500000,'Cadell']
dunstan = [30000, 2000, 7, 55000, 'Dunstan']
hyde = [25000, 10000, 2, 52000, 'Hyde'] # Just estimated for now
lake_george = [750000, 3000, 30, 3000000, 'Lake George']
lake_edgar = [1000000, 16500, 5, 2000000, 'Lake Edgar']
hyden = [140000, 30000, 3, 250000, 'Hyden']

wharekuri = [65000, 7500, 5, 150000, 'Wharekuri']
wairau_cloudy_bay = [2200, 1000, 5, 6000, 'Wairau CB'] # Nicol and Van Dissen 2018
wairau_full_rupture = [1077, 1077, 5, 6000, 'Wairau Full'] 
wellington = [1575, 715, 4, 4900, 'Wellington'] # Langridge et al 2011

galeen=[30000, 5000, 5, 99000, 'Galeen']
roer = [95000, 600, 5, 101000, 'Roer']

alpine = [330, 330, 27, 8000, 'Alpine']
san_jacinto_anza = [254, 254, 15, 4000, 'San Jacinto'] # Rockwell et al 2015

meers = [120000, 1700, 2, 120000, 'Meers']

data = np.array([akatore, cadell, dunstan, hyde, lake_edgar, hyden,
                 galeen, roer, wharekuri, lake_george,
                 wairau_cloudy_bay, wairau_full_rupture,
                 alpine, san_jacinto_anza,
                 wellington, meers])
labels = data[:,4]
data = data[:,0:-1]
data = data.astype(np.float)
#print(data)
#print(labels)
x = data[:,0]
y = data[:,1]
#print(x)
#print(y)
plt.scatter(x,y)
ax = plt.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Quiescent recurrence (years)', fontsize=18)
ax.set_ylabel('Active recurrence (years)', fontsize=18)
texts = []                                                                                                                      
for i, label in enumerate(labels):                                                                                              
    texts.append(plt.text(x[i], y[i], label))                                                                                   
adjust_text(texts)                                                                                                              
plt.tight_layout()
plt.savefig('active_vs_quiescence.png')

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
texts = []                                                                                                                      
for i, label in enumerate(labels):                                                                                              
    texts.append(plt.text(x[i], y[i], label))                                                                                   
adjust_text(texts)                                                                                                               
plt.tight_layout()
plt.savefig('number_event_vs_time.png')  

x = data[:,0]                                                                                                                   
y = data[:,2]                                                                                                                   
print(x)                                                                                                                        
print(y)                                                                                                                        
plt.clf()                                                                                                                       
plt.scatter(x,y)
ax = plt.gca()                                                                                                                  
ax.set_xlabel('Quiescent Recurrence (years)', fontsize=18)                                                                                   
ax.set_ylabel('Number of events', fontsize=18)
ax.set_xscale('log')  
texts = []                                                                                                                      
for i, label in enumerate(labels):                                                                                               
    texts.append(plt.text(x[i], y[i], label))                                                                                   
adjust_text(texts)
plt.tight_layout() 
plt.savefig('number_event_vs_quiescent_period.png')


x = data[:,1]                                                                                                                   
y = data[:,2]                                                                                                                   
print(x)                                                                                                                        
print(y)                                                                                                                        
plt.clf()                                                                                                                       
plt.scatter(x,y)                                                                                                                
ax = plt.gca()
ax.set_xscale('log')
ax.set_xlabel('Active Recurrence (years)', fontsize=18)                                                                                   
ax.set_ylabel('Number of events', fontsize=18)
texts = []                                                                                                                       
for i, label in enumerate(labels):                                                                                              
    texts.append(plt.text(x[i], y[i], label))                                                                                    
adjust_text(texts)                                                                                                               
plt.tight_layout()
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
ax.set_xlabel('Quiescent Recurrence (years)', fontsize=18)
ax.set_ylabel('Events per year', fontsize=18)
texts = []                                                                                                                     
for i, label in enumerate(labels):                                                                                              
    texts.append(plt.text(x[i], y[i], label))
adjust_text(texts)
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
ax.set_xlabel('Active Recurrence (years)', fontsize=18)
ax.set_ylabel('Events per year', fontsize=18)
texts = []
for i, label in enumerate(labels):
    texts.append(plt.text(x[i], y[i], label))
#    ax.annotate(label, (x[i], y[i]))
adjust_text(texts)#, arrowprops=dict(arrowstyle="->", color='r', lw=0.5))
plt.tight_layout()
plt.savefig('event_per_year_vs_active_period.png')

x = data[:,1]/data[:,0]
y = data[:,2]/data[:,3]                                                                                                                                   
print(x)
print(y)
plt.clf()
plt.scatter(x,y)
ax = plt.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim([0.9*min(y), 1.1*max(y)])                                                                                                    
ax.set_xlim([0.9*min(x), 1.1*max(x)]) 
ax.set_xlabel('Ratio Active:Quiescent', fontsize=18)                                                                                                  
ax.set_ylabel('Events per year', fontsize=18)                                                                                                             
texts = []
for i, label in enumerate(labels):
    if label == 'Akatore' or label == 'Dunstan' or label == 'Hyde':         
        texts.append(plt.text(x[i], y[i], label, color='r'))
    else:
        texts.append(plt.text(x[i], y[i], label, color='k'))
#    ax.annotate(label, (x[i], y[i]))                                                                                                                 
adjust_text(texts)#, arrowprops=dict(arrowstyle="->", color='r', lw=0.5))
plt.tight_layout()                                                                                                                                     
plt.savefig('event_per_year_vs_ratio_active_quiescent.png') 
