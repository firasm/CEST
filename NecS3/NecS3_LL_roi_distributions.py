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

with open('/Volumes/Data/Dropboxes/PhD./Dropbox/studies/analysis/NecS3/NecS3.json','r') as infile:
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
    
        sarpy.Scan(master_sheet[k][key_list[2]][0]).adata[key_list[6]].export2nii('LL1.nii')
        sarpy.Scan(master_sheet[k][key_list[3]][0]).adata[key_list[6]].export2nii('LL2.nii')
        
        roi1 = sarpy.Scan(master_sheet[k][key_list[4]][0]).adata[key_list[7]]
        roi1.export2nii('roi1.nii')
        
        roi2 = sarpy.Scan(master_sheet[k][key_list[5]][0]).adata[key_list[7]]
        roi2.export2nii('roi2.nii')
      
        LLresample1 = sarpy.ImageProcessing.resample_onto.resample_onto('LL1.nii','roi1.nii')      
        LLresample2 = sarpy.ImageProcessing.resample_onto.resample_onto('LL2.nii','roi2.nii')
        
#        LL1n = nibabel.Nifti1Image(LLresample1,numpy.eye(4))
#        LL1n.to_filename('LL1n.nii')
        
        bins = numpy.linspace(500,3500,20)
        name1 = 'T1distrib-' + k
        #name2 = 'T1distrib2-' + k

        pylab.close('all')
        sarpy.fmoosvi.wrappers.roi_distribution(LLresample1, roi1, bins)
        sarpy.fmoosvi.wrappers.roi_distribution(LLresample2, roi2, bins, save_histogram = True, save_name = name1)
        

                
        print('Success')

    except:
        pylab.close('all')

        print('Unknown error {0}'.format(k))
        
