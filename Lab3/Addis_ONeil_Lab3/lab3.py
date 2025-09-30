import numpy as np
import matplotlib.pyplot as plt
import os
from astropy.io import fits
from photutils.aperture import CircularAnnulus, CircularAperture, ApertureStats, aperture_photometry
from centroids import find_star_pos

directory_path = 'astro344/observation_data/lab3'

def std(x):
    a = np.sum([(i - np.mean(x))**2 for i in x])
    return np.sqrt((1 / (len(x) - 1)) * a)

# 1 B filter, A, X5, X6     see 3 /home/colton/school/classes/astro344/observation_data/lab3/cookedfilB_A_x6_7-5_00x.FIT
# 2 B filter, B, X2, X3, X8 see 5 /home/colton/school/classes/astro344/observation_data/lab3/cookedfilB_B_x2_7-5_00x.FIT
# 3 V filter, A, X5, X6     see 1 /home/colton/school/classes/astro344/observation_data/lab3/cookedfilV_A_x6_7-5_00x.FIT
# 4 bad
# 5 V filter, B, X2, X3, X8 see 2 /home/colton/school/classes/astro344/observation_data/lab3/cookedfilV_B_xx_7-5_00x.FIT

reference_mag = {'A' : {'B' : 6.558, 'V' : 6.569}, 'B' : {'B' : 6.786, 'V' : 6.821}}

#                               A          x5        x6                   B           x2          x3         x8
reference_coords = {'A' : [[404, 404], [114, 216], [714, 278]], 'B' : [[244, 294], [472, 419], [644, 151], [180, 29]]}

results = {'B' : {'A' : {}, 'B' : {}}, 'V' : {'A' : {}, 'B' : {}}}

def find_magnitudes(file_name, guess, ref_mag, ref_filt):
    #get positions and radius
    positions, radius = find_star_pos(file_name, guess)
    aperture = CircularAperture(positions, r=radius)
    annulus_aperture = CircularAnnulus(positions, r_in=radius + 5, r_out=radius + 10)

    #get image data
    image_data, image_hdr  = fits.getdata(file_name, ext=0, header=True)
    image_data  = np.float64(image_data)

    #star
    phot_table = aperture_photometry(image_data,aperture)
    #background
    aperstats = ApertureStats(image_data, annulus_aperture)
    bkg_mean = aperstats.mean

    aperture_area = aperture.area_overlap(image_data)
    total_bkg = bkg_mean*aperture_area

    # Do Background Subtracted Photometry
    phot_bksub = phot_table['aperture_sum'] - total_bkg

    # Add items to the photometry table 
    phot_table['Background'] = total_bkg
    phot_table['BG Sub Phot'] = phot_bksub
    phot_table['Inst. Mag'] = -2.5*np.log10(phot_bksub)

    # Calculate the zeropoint then add mags to table and print
    zp = ref_mag - phot_table['Inst. Mag'][0]
    phot_table[f'{ref_filt} Mag'] = phot_table['Inst. Mag'] + zp
    #print(phot_table[f'{ref_filt} Mag'])
    return phot_table[f'{ref_filt} Mag']

#go through files
for file in os.listdir(os.fsencode(directory_path)):
    filename = os.fsdecode(file)
    filename_parts = filename.strip().split('_')
    star_loc = filename_parts[1]
    light_filter = filename_parts[0][-1]
    if filename_parts[0].find('cooked') > -1 and filename_parts[3] == '7-5':
        mag_table = find_magnitudes(directory_path + '/' + filename, reference_coords[star_loc], reference_mag[star_loc][light_filter], light_filter)
        for i in range(len(mag_table)):
            coords = reference_coords[star_loc][i][0]
            if coords in results[light_filter][star_loc]:
                results[light_filter][star_loc][coords].append(mag_table[i])
            else:
                results[light_filter][star_loc][coords] = [mag_table[i]]

for key, value in results.items():
    print(key)
    for key1, value1 in value.items():
        print(key1)
        for key2, value2 in value1.items():
            print(str(key2) + ': ' + str(np.array(value2)))
            print("   " + str(np.mean(value2)) + ":::" + str(std(value2) / np.sqrt(len(value2))))

