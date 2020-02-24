# Estimate Be10/Be9 ratio assuming a carrier dose of 300 um Be9

# Estimated age of sample (years)
age = 10000
# Estimated production rate (atoms g^-1 yr^-1)
rate = 3.84 # Sea level rate from Putnam et al. 2010
# In situ cosmogenic 10Be production-rate calibration from
# the Southern Alps, New Zealand. Quaternary Geochronology 5(4).

# Estimated quartz content (g)
quartz_weight = 20
# Avagadro's number
ava = 6.022e23
# Weight of single Be9 atom based on molar mass of 81 g
be9w=81.109638/ava
# Number of atoms in 3 um of carrier Be9
be9n = 0.0003/be9w
be10n = rate*age*quartz_weight
ratio = be10n/be9n
print('Ratio of Be10 to Be9 atoms ', ratio)
