# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 23:48:38 2013

@author: firas
"""

##### ROI testing

import sarpy
import sarpy.fmoosvi.getters
import sarpy.fmoosvi.wrappers
import sarpy.ImageProcessing.resample_onto

import pylab
import numpy
import json
import nibabel

with open('/Volumes/Data/Dropboxes/PhD./Dropbox/studies/analysis/NecS1/NecS1.json','r') as infile:
    master_sheet = json.load(infile)

data_list = []
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

        LL_1 = sarpy.Scan(master_sheet[k][key_list[2]][0]).adata[key_list[6]]
        LL_2 = sarpy.Scan(master_sheet[k][key_list[3]][0]).adata[key_list[6]]    

        roi1 = sarpy.Scan(master_sheet[k][key_list[4]][0]).adata[key_list[7]]
        roi2 = sarpy.Scan(master_sheet[k][key_list[5]][0]).adata[key_list[7]]
        
        roi_mask1 = sarpy.fmoosvi.analysis.h_image_to_mask(roi1)
        roi_mask2 = sarpy.fmoosvi.analysis.h_image_to_mask(roi2)

        roi_resample1 = sarpy.ImageProcessing.resample_onto.resample_onto_pdata(roi_mask1,LL_1)      
        roi_resample2 = sarpy.ImageProcessing.resample_onto.resample_onto_pdata(roi_mask2,LL_2)

        curr_T1s1 = LL_1.data * roi_resample1 
        curr_T1s2 = LL_2.data * roi_resample2
              
        bins = numpy.linspace(500,3500,20)
        name1 = 'T1distrib1-' + k
        name2 = 'T1distrib2-' + k

        pylab.close('all')
        sarpy.fmoosvi.wrappers.roi_distribution(LL_1, roi_resample1, bins)
        sarpy.fmoosvi.wrappers.roi_distribution(LL_2, roi_resample2, bins, save_histogram = True, save_name = name1)
                      
        print('Success')
    except KeyError:
        print('Key error {0}'.format(k))


    except IOError:        
        print('Unknown error {0}'.format(k))
        raise
