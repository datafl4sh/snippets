## A little program to compute the budget of a microwave link
## between two sites, with references for who wants to know
## from where all those formulas arise.
##
## Matteo Cicuttin 2017 - IV3IWE - This code is in public domain.

from math import log10, sqrt

frequency       = 2450  # Frequency in MHz
distance        = 10    # Distance in kilometers
foliage_depth   = 0     # ITU-R foliage model, 0-400 meters

#### Site A parameters ####
pwr_tx_a        = 27    # Transmitter power in dBm
gain_a          = 10    # Antenna gain in dB
cable_loss_a    = 0     # Cable losses in dB
other_losses_a  = 0     # Other losses in dB
sensitivity_a   = -88   # Receiver sensitivity in dBm
reqd_margin_a   = 12    # Required margin in dB

#### Site B parameters ####
pwr_tx_b        = 27    # Transmitter power in dBm
gain_b          = 10    # Antenna gain in dB
cable_loss_b    = 0     # Cable losses in dB
other_losses_b  = 0     # Other losses in dB
sensitivity_b   = -88   # Receiver sensitivity in dBm
reqd_margin_b   = 12    # Required margin in dB

# https://en.wikipedia.org/wiki/Free-space_path_loss
fspl = 20*log10(distance) + 20*log10(frequency) + 32.45

# Y. S. Meng, Y. H. Lee, Investigations of foliage effect
# on modern wireless communication systems: a review
foliage_loss = 0.2*(frequency**0.3)*(foliage_depth**0.6)

# https://en.wikipedia.org/wiki/Link_budget
gains       = gain_a + gain_b
losses      = fspl + foliage_loss + cable_loss_a + other_losses_a + \
              cable_loss_b + other_losses_b

pwr_at_a    = pwr_tx_b + gains - losses
pwr_at_b    = pwr_tx_a + gains - losses

margin_at_a = pwr_at_a - sensitivity_a
margin_at_b = pwr_at_b - sensitivity_b

# https://en.wikipedia.org/wiki/Earth_bulge
earth_bulge = (distance/4.12)**2

# https://en.wikipedia.org/wiki/Fresnel_zone
f1_radius = 8.656*sqrt(distance/(frequency/1000))

min_antenna_height = earth_bulge + f1_radius

print("LINK BUDGET COMPUTATIONS")
print("  Total gains:       %8.3f dB" % gains)
print("  Total losses:      %8.3f dB" % losses)
print("    Path loss:       %8.3f dB" % fspl)
print("    Foliage loss:    %8.3f dB" % foliage_loss)
print("\n                    SITE A             SITE B")
print("Avail. power:   %8.3f dBm" % pwr_at_a,  "      %8.3f dBm" % pwr_at_b)
print("Link margin:    %8.3f dBm" % margin_at_a,  "      %8.3f dBm" % margin_at_b)

if (margin_at_a < reqd_margin_a) and (margin_at_a > 0):
    print("WARNING: margin at site A is less than required.")
if margin_at_a <= 0:
    print("WARNING: insufficient signal at site A. Link is likely to NOT work.")

if (margin_at_b < reqd_margin_b) and (margin_at_b > 0):
    print("WARNING: margin at site B is less than required.")
if margin_at_b <= 0:
    print("WARNING: insufficient signal at site B. Link is likely to NOT work.")

print("\nANTENNA HEIGHT COMPUTATIONS")
print("  Earth bulge:               %.3f meters" % earth_bulge)
print("  F1 Fresnel zone radius:    %.3f meters" % f1_radius)
print("  Minimum antenna height:    %.3f meters" % min_antenna_height)



