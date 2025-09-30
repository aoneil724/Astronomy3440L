import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from photutils.centroids import centroid_2dg, centroid_sources, centroid_com
from photutils.aperture import CircularAnnulus, CircularAperture, ApertureStats, aperture_photometry

def find_star_pos(image_file, estimate):
    # Use a dark frame median subtracted fit file
    image_data, image_hdr  = fits.getdata(image_file, ext=0, header=True)
    image_data  = np.float64(image_data)

    # Remember to subtract 1 from the aperture x and y pixel values if you are using ds9 to estimate aperture positions
    centroid_data = image_data

    # Specify a region where no stars exist and its just background, y1:y2, x1:x2

    blank_region =centroid_data[400:500, 600:700]

    centroid_data -= np.median(blank_region)

    # Rough position of stars you want to measure (use ds9)
    x_init, y_init = map(list, zip(*estimate))

    #high box size?
    x, y = centroid_sources(centroid_data, x_init, y_init, box_size=35, centroid_func=centroid_2dg)
    x-=1
    y-=1

    #turn x and y into integers
    rounded_x = np.round(x)
    rounded_y = np.round(y)

    # find the radius of the aperature (this is just looking at one star and generalizing)
    starting_val = centroid_data[int(rounded_y[0]), int(rounded_x[0])]
    rad = 0
    for z in range(1000):
        val = centroid_data[int(rounded_y[0]), int(rounded_x[0]) + z]
        if val < (starting_val / 10):
            rad += z-1
            break

    # Position of Stars
    position = []
    for c in range(len(x)):
        position.append((x[c],y[c]))

    print(rad)
    if rad == 8:
        print(image_file)
        aperture = CircularAperture(position, r=rad + 2)
        annulus_aperture = CircularAnnulus(position, r_in=rad + 5, r_out=rad + 10)

        plt.figure(figsize=(15,20))
        plt.imshow(image_data, cmap='gray', vmin=10, vmax=1500)
        plt.xlim(0, 764)
        plt.ylim(0, 509)

        aperture.plot(color='white', lw=2)
        annulus_aperture.plot(color='red', lw=2)

        plt.show()
    return position, rad
