# Estimate Be10/Be9 ratio assuming a carrier dose of 300 um Be9
import numpy as np

# Estimated age of sample (years)
age = 2000000
# Estimated production rate (atoms g^-1 yr^-1)
rate = 4 #3.84 # Sea level rate from Putnam et al. 2010
# In situ cosmogenic 10Be production-rate calibration from
# the Southern Alps, New Zealand. Quaternary Geochronology 5(4).
fsp = 0.974 #Ignore for now as near 1 # Fraction of production rate due to spallation from
#Stone, J. O. (2000). Air pressure and cosmogenic isotope production,
#J. Geophys. Res., 105, no. B10, 23753–23759.
# Remaining fraction (1-fsp) is due to muon production

# Estimated quartz content (g)
quartz_weight = 20
h = 1360 # Elevation above sea level (m)
#Er = 10 # Erosion rate mm/kyr (Schist)
Er = 0.4 # Erosion rate mm/kyr (Silcrete) Values from Bennett et al 2005 
Ps = 1013.5 # Mean sea level pressure hPa
Ts = 288.15 # Mean sea level temperature K
eta = 0.0065 # Adiabatic lapse rate K m^-1
gM_r = 0.03417 # Combined molar weight, gravity and gas constants
Pz = Ps*np.exp((-1/eta)*gM_r*(np.log(Ts) - np.log(Ts-eta*h)))
print('Pressure at elevation %.1f masl is %.1f hPa' % (h, Pz))

# Other constants
decay_constant = np.log(2)/1360000 #1/yr
attenuation_length = 160 #g/cm^2
density = 2.8 #g/cm^3
burial_tau = 0.0 # Assumed for now

# Scaling constants for 40 degrees latitude
a40 = 56.7733
b40 = 649.1314
c40 = -0.160859
d40 = 1.5463e-04
e40 = -5.0330e-08

# Scaling constants for 50 degrees latitude
a50 = 69.0720
b50 = 832.4566
c50 = -0.199252
d50 = 1.9391e-4
e50 = -6.3653e-8

S_lambda40 = a40+b40*np.exp(-1*Pz/150) + c40*Pz + d40*Pz**2 + e40*Pz**3
S_lambda50 = a50+b50*np.exp(-1*Pz/150) + c50*Pz + d50*Pz**2 + e50*Pz**3
S_lambda45 = (S_lambda40 + S_lambda50)/2 # Linear interpolation

print('Spallation scaling factor %.2f' % S_lambda45)
#rate_s = rate*S_lambda45

# Calulate scaled Be10 concentration per gram
be10_s = (S_lambda45*rate/(decay_constant + (density*Er/10/1000)/attenuation_length)) * \
    (1 - np.exp(-1*(decay_constant + (density*Er/10/1000)/attenuation_length)*age)) * \
    np.exp(-1*(decay_constant*burial_tau))

#print('Total scaling ', (rate_s/rate))
print('Scaled concentration', be10_s)
# Avagadro's number
ava = 6.022e23
# Weight of single Be9 atom based on molar mass of 81 g
be9w=81.109638/ava
# Number of atoms in 3 um of carrier Be9
be9n = 0.0003/be9w
# Number of atoms in 250 ug of carrier Be9
be9n = 0.00025/be9w
# Try Klaus's calculation
be9n = 0.00025 / 9.012 * ava
be10n = rate*age*quartz_weight
ratio = be10n/be9n
print('Ratio of Be10 to Be9 atoms ', ratio)
be10n = be10_s*quartz_weight
print('Expected B10 atoms, scaled for altitude', be10n)
ratio = be10n/be9n
print('Ratio of Be10 to Be9 atoms, scaled for altitude, latitude, decay and erosion ', ratio) 
