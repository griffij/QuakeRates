# Estimate Be10/Be9 ratio assuming a carrier dose of 300 um Be9
import numpy as np

# Estimated age of sample (years)
age = 20000
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
h = 350 # Elevation above sea level (m)

Ps = 1013.5 # Mean sea level pressure hPa
Ts = 288.15 # Mean sea level temperature K
eta = 0.0065 # Adiabatic lapse rate K m^-1
gM_r = 0.03417 # Combined molar weight, gravity and gas constants
Pz = Ps*np.exp((-1/eta)*gM_r*(np.log(Ts) - np.log(Ts-eta*h)))
print('Pressure at elevation %.1f masl is %.1f hPa' % (h, Pz))

# Scaling constants for 40 degrees latitud
a40 = 56.7733
b40 = 649.1314
c40 = -0.160859
d40 = 1.5463e-04
e40 = -5.0330e-08

S_lambda = a40+b40*np.exp(-1*Pz/150) + c40*Pz + d40*Pz**2 + e40*Pz**3
print('Spallation scaling factor %.2f' % S_lambda)
rate_s = rate*S_lambda


# Avagadro's number
ava = 6.022e23
# Weight of single Be9 atom based on molar mass of 81 g
be9w=81.109638/ava
# Number of atoms in 3 um of carrier Be9
be9n = 0.0003/be9w
be10n = rate*age*quartz_weight
ratio = be10n/be9n
print('Ratio of Be10 to Be9 atoms ', ratio)
be10n = rate_s*age*quartz_weight
print('Expected B10 atoms, scaled for altitude', be10n)
ratio = be10n/be9n
print('Ratio of Be10 to Be9 atoms, scaled for altitude ', ratio) 
