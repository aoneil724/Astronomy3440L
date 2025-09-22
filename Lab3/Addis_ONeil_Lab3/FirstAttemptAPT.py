import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from photutils.aperture import CircularAnnulus, CircularAperture
from photutils.aperture import ApertureStats
from photutils.aperture import aperture_photometry
from photutils.datasets import make_4gaussians_image

from photutils.centroids import centroid_2dg, centroid_sources, centroid_com

image_file = "filV_A_x6_7-5_005.FIT"
image_data, image_hdr  = fits.getdata(image_file, ext=0, header=True)
# this is required for images obtained by the ISU CCDs
image_data  = np.float64(image_data)

# Enter Aperture and Annulus Information
#
# Remember to subtract 1 from the aperture x and y pixel values if you are using ds9 to estimate aperture positions

centroid_data = image_data

# Specify a region where no stars exist and its just background, y1:y2, x1:x2
blank_region =centroid_data[350:450, 0:200]

centroid_data -= np.median(blank_region)

# Rough position of stars you want to measure
x_init = (400, 56)

y_init = (400, 205)

x, y = centroid_sources(centroid_data, x_init, y_init, box_size=25,

                        centroid_func=centroid_2dg)

x-=1
y-=1

rounded_x = np.round(x)
rounded_y = np.round(y)
starting_val = centroid_data[int(rounded_x[0]), int(rounded_y[0])]
rad = 0
ind = 0
for z in range(1000):
    val = centroid_data[int(rounded_x[0])+z, int(rounded_y[0])]
    print(starting_val, val)
    if val < (starting_val / 4):
        rad += z-1
        break
    print(z)
position = []
for c in range(len(x)):
    position.append((x[c],y[c]))
print(rad)
aperture = CircularAperture(position, r=rad)
annulus_aperture = CircularAnnulus(position, r_in=(rad + rad/2), r_out=(2*rad))
plt.imshow(image_data, cmap='gray', vmin=1640, vmax=2230)

plt.xlim(0, 764)
plt.ylim(0, 509)

aperture.plot(color='white', lw=2)
annulus_aperture.plot(color='red', lw=2)

plt.show()