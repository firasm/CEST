# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 23:48:38 2013

@author: firas
"""

##### ROI testing

import sarpy
import sarpy.fmoosvi.analysis
import sarpy.fmoosvi.getters
import sarpy.fmoosvi.wrappers
import sarpy.ImageProcessing.resample_onto

import pylab
import numpy
import json

with open('/Volumes/Data/Dropboxes/PhD./Dropbox/studies/analysis/NecS3/NecS3.json','r') as infile:
    master_sheet = json.load(infile)

data_list = {}
T1_vals_day1 = {}
T1_vals_day2= {}

day1_B = []
day2_B = []


a = 2294.623591
2432.317817
2773.823654
2467.396795
2336.813094
2276.250102
2351.946978
2404.542488
2224.638317
3066.115236

    
    try:
        data_list = []
        key_list = []
        key_list.append(k) # 0
        key_list.append('T1map-') #1
        key_list.append('0h-LL') #2
        key_list.append('24h-LL') #3
        key_list.append('0h-IR_A') #4
        key_list.append('24h-IR_A') #5
        key_list.append('T1map_LL') #6
        key_list.append('IR_tumour_rois') #7
        
        ## I have  feeling that the animals that got Omniscan on Day 1 had some left over 24h later
#        if (k == 'NecS3Hs02' or
#            k == 'NecS3Hs04' or
#            k == 'NecS3Hs05' or
#            k == 'NecS3Hs06'):
#            
#            #Belay that, doesn't actually seem to matter too much
#            continue

        # Get the calculated T1 maps from the adata
        LL_1 = sarpy.Scan(master_sheet[k][key_list[2]][0]).adata[key_list[6]]
        LL_2 = sarpy.Scan(master_sheet[k][key_list[3]][0]).adata[key_list[6]]    
      
        # Get the stored ROIs
        roi1 = sarpy.Scan(master_sheet[k][key_list[4]][0]).adata[key_list[7]]
        roi2 = sarpy.Scan(master_sheet[k][key_list[5]][0]).adata[key_list[7]]

        # Turn the roi images into masks        
        roi_mask1 = sarpy.fmoosvi.analysis.h_image_to_mask(roi1)
        roi_mask2 = sarpy.fmoosvi.analysis.h_image_to_mask(roi2)

        # Resample the ROI image to match the T1 parameter map      
        roi_resample1 = sarpy.ImageProcessing.resample_onto.resample_onto_pdata(roi_mask1,LL_1)      
        roi_resample2 = sarpy.ImageProcessing.resample_onto.resample_onto_pdata(roi_mask2,LL_2)

        # Getting the T1s in an ROI and add them to the aggregate list
        curr_T1s1 = LL_1.data * roi_resample1 
        curr_T1s2 = LL_2.data * roi_resample2
        
        T1_vals_day1[k]=list(curr_T1s1[numpy.isfinite(curr_T1s1)])
        T1_vals_day2[k]=list(curr_T1s2[numpy.isfinite(curr_T1s2)])

        day1_B.append(numpy.mean(curr_T1s1[numpy.isfinite(curr_T1s1)]))
        day2_B.append(numpy.mean(curr_T1s2[numpy.isfinite(curr_T1s2)]))


        # Drawing and saving plots
#        pylab.close('all')
#        name = 'T1roi-' + k
#        bins = numpy.linspace(200,3500,50)
#        sarpy.fmoosvi.wrappers.roi_distribution(LL_1, roi_resample1, bins)
#        sarpy.fmoosvi.wrappers.roi_distribution(LL_2, roi_resample2, bins, save_histogram = True, save_name = name)
#        pylab.title(name)
 
        # Print the number of pixels used in each study
        print k
        print('\t {0} had {1} pixels on day1'.format(master_sheet[k][key_list[4]][0], curr_T1s1[numpy.isfinite(curr_T1s1)].shape))      
        print('\t {0} had {1} pixels on day2'.format(master_sheet[k][key_list[5]][0], curr_T1s2[numpy.isfinite(curr_T1s2)].shape))
        
    except KeyError:
        print('Key error, please ignore {0}'.format(k))
        
#Now squash the lists together, and turn into an array use an incomprehensible list comprehension, as named by someone on SE
#T1_vals_day1 = [item for sublist in T1_vals_day1 for item in sublist]
#T1_vals_day1 = numpy.array(T1_vals_day1)
#
#T1_vals_day2 = [item for sublist in T1_vals_day2 for item in sublist]
#T1_vals_day2 = numpy.array(T1_vals_day2)
#
#pylab.figure()
#bins = numpy.linspace(200,3500,50)
#pylab.hist(T1_vals_day1, bins, alpha=0.6)
#pylab.hist(T1_vals_day2, bins, alpha=0.6)
#pylab.title('RF1106 animals from NecS3')
#
### Export the data to a JSON file:
#
with open('/Volumes/Data/Dropboxes/PhD./Dropbox/studies/analysis/NecS3/NecS3_T1_vals_day1.json', 'w') as f:
    json.dump(T1_vals_day1, f)
with open('/Volumes/Data/Dropboxes/PhD./Dropbox/studies/analysis/NecS3/NecS3_T1_vals_day2.json', 'w') as f:
    json.dump(T1_vals_day2, f)

with open('/Volumes/Data/Dropboxes/PhD./Dropbox/studies/analysis/NecS3/day1_B.json', 'w') as f:
    json.dump(day1_B, f)
with open('/Volumes/Data/Dropboxes/PhD./Dropbox/studies/analysis/NecS3/day2_B.json', 'w') as f:
    json.dump(day2_B, f)





