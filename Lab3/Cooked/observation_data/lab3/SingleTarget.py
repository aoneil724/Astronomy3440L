##############################################################################
# Import packages
#

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from photutils.aperture import CircularAnnulus, CircularAperture
from photutils.aperture import ApertureStats
from photutils.aperture import aperture_photometry


##############################################################################
# Get the image files
#

#M39 Image
image_file = "cookedfilB_B_x2_7-5_005.FIT"
image_data, image_hdr  = fits.getdata(image_file, ext=0, header=True)
# this is required for images obtained by the ISU CCDs
image_data  = np.float64(image_data)

##############################################################################
# Enter Aperture and Annulus Information
#
# Remember to subtract 1 from the aperture x and y pixel values if you 
#    are using ds9 to estimate aperture positions

#M39 Positions
positions = [(203,20), (257, 273)]
aperture = CircularAperture(positions, r=8)
annulus_aperture = CircularAnnulus(positions, r_in=12, r_out=16)

##############################################################################
# Do aperture photometry and output results
#

# Do photometry in circular aperture
#
phot_table = aperture_photometry(image_data,aperture)

#Get mean background level in annulus
aperstats = ApertureStats(image_data, annulus_aperture)
bkg_mean = aperstats.mean

# Can also use median
#bkg_median = aperstats.median
#print("Annulus Median", bkg_median)

#Get aperture area
# 
aperture_area = aperture.area_overlap(image_data)

# Calculate Total Background
#
total_bkg = bkg_mean*aperture_area

# Do Background Subtracted Photometry
#
phot_bksub = phot_table['aperture_sum'] - total_bkg

# Add items to the photometry table for output
#
phot_table['Background'] = total_bkg
phot_table['BG Sub Phot'] = phot_bksub
phot_table['Inst. Mag'] = -2.5*np.log10(phot_bksub)

# Calculate the zeropoint then add V mags to table and print
#

zp = 6.786 - phot_table['Inst. Mag'][1]
print("Zero Point:", zp)


phot_table['V Mag'] = phot_table['Inst. Mag'] + zp

print("")
print(phot_table)


##############################################################################
# Display image with apertures and annuli added
#
# Can used ds9 to estimate vmin and vmax stretch values
#
# Adjust xlim and ylim to match your image size
#

plt.imshow(image_data, cmap='gray', vmin=0, vmax=500)
plt.xlim(0, 764)
plt.ylim(0, 509)

aperture.plot(color='white', lw=2)
annulus_aperture.plot(color='red', lw=2)

plt.show()