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

data_list = []
T1_vals_day1 = []
T1_vals_day2= []

for k,v in master_sheet.iteritems():
    
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
    
#        LL_1 = sarpy.Scan(master_sheet[k][key_list[2]][0]).adata[key_list[6]]
#        LL_2 = sarpy.Scan(master_sheet[k][key_list[3]][0]).adata[key_list[6]]

        LL_1 = sarpy.Scan(master_sheet[k][key_list[2]][0]).adata[key_list[6]]
        LL_2 = sarpy.Scan(master_sheet[k][key_list[3]][0]).adata[key_list[6]]    

        roi1 = sarpy.Scan(master_sheet[k][key_list[4]][0]).adata[key_list[7]]
        roi2 = sarpy.Scan(master_sheet[k][key_list[5]][0]).adata[key_list[7]]

        roi_resample1 = sarpy.ImageProcessing.resample_onto.resample_onto_pdata(roi1,LL_1)      
        roi_resample2 = sarpy.ImageProcessing.resample_onto.resample_onto_pdata(roi2,LL_2)

        curr_T1s1 = LL_1 * sarpy.fmoosvi.analysis.h_image_to_mask(roi_resample1)
        curr_T1s2 = LL_2 * sarpy.fmoosvi.analysis.h_image_to_mask(roi_resample1)

        T1_vals_day1.append(numpy.isfinite(curr_T1s1.flatten()))
        T1_vals_day2.append(numpy.isfinite(curr_T1s2.flatten()))
        
    except:
        pylab.close('all')
        print('Unknown error {0}'.format(k))
        raise
        
# Now squash the lists together, and turn into an array

T1_vals_day1 = nump.array(T1_vals_day1).flatten()
T1_vals_day2 = nump.array(T1_vals_day2).flatten()        
        
        
        
        
        
        
        
        
        
        
        
