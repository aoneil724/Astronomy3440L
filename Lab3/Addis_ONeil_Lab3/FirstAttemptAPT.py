import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from photutils.aperture import CircularAnnulus, CircularAperture
from photutils.aperture import ApertureStats
from photutils.aperture import aperture_photometry
from photutils.datasets import make_4gaussians_image

from photutils.centroids import centroid_2dg, centroid_sources, centroid_com

# Use a dark frame median subtracted fit file

image_file = "filV_A_x6_7-5_005.FIT"
image_data, image_hdr  = fits.getdata(image_file, ext=0, header=True)
image_data  = np.float64(image_data)


# Remember to subtract 1 from the aperture x and y pixel values if you are using ds9 to estimate aperture positions

centroid_data = image_data

# Specify a region where no stars exist and its just background, y1:y2, x1:x2
blank_region =centroid_data[350:450, 0:200]

centroid_data -= np.median(blank_region)

# Rough position of stars you want to measure (use ds9)
x_init = (400, 56)

y_init = (400, 205)

x, y = centroid_sources(centroid_data, x_init, y_init, box_size=25,

                        centroid_func=centroid_2dg)
x-=1
y-=1

#turn x and y into integers
rounded_x = np.round(x)
rounded_y = np.round(y)

# find the radius of the aperature
starting_val = centroid_data[int(rounded_x[0]), int(rounded_y[0])]
rad = 0
for z in range(1000):
    val = centroid_data[int(rounded_x[0])+z, int(rounded_y[0])]
    if val < (starting_val / 4):
        rad += z-1
        break

# Position of Stars
position = []
for c in range(len(x)):
    position.append((x[c],y[c]))


aperture = CircularAperture(position, r=rad)
annulus_aperture = CircularAnnulus(position, r_in=(rad + rad/2), r_out=(2*rad))

phot_table = aperture_photometry(image_data,aperture)

#Get mean background level in annulus
aperstats = ApertureStats(image_data, annulus_aperture)
bkg_mean = aperstats.mean

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

zp = zeropoint_value - phot_table['Inst. Mag'][1]
print("Zero Point:", zp)


phot_table['V Mag'] = phot_table['Inst. Mag'] + zp

print("")
print(phot_table)

#plot
plt.imshow(image_data, cmap='gray', vmin=1640, vmax=2230)

plt.xlim(0, 764)
plt.ylim(0, 509)

aperture.plot(color='white', lw=2)
annulus_aperture.plot(color='red', lw=2)

plt.show()