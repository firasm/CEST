# -*- coding: utf-8 -*-
"""
Created on Sat Apr 13 17:01:59 2013

@author: fmoosvi
"""

import sarpy 
import sarpy.fmoosvi.analysis
import sarpy.ImageProcessing.resample_onto
import pylab
import json

with open('/Volumes/Data/Dropboxes/PhD./Dropbox/studies/analysis/NecS3/NecS3.json','r') as infile:
    master_sheet = json.load(infile)

#k = 'NecS3Hs11'

for k,v in master_sheet.iteritems():

    try:
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
        roi1 = sarpy.Scan(master_sheet[k][key_list[4]][0]).adata[key_list[7]]
        roi1_mask = sarpy.fmoosvi.analysis.h_image_to_mask(roi1)
        roi_resample1 = sarpy.ImageProcessing.resample_onto.resample_onto_pdata(roi1_mask,LL_1)
        curr_T1s1 = LL_1.data * roi_resample1
        
        fig = pylab.figure(figsize = (14,5))
        
        fig.add_subplot(241)
        pylab.imshow(LL_1.parent.data[:,:,3,0])
        pylab.title('LL data (original)')
        
        fig.add_subplot(242)
        pylab.imshow(LL_1.data[:,:,3])
        pylab.title('T1 map')
        
        fig.add_subplot(243)
        pylab.imshow(roi1.parent.data[:,:,3])
        pylab.title('Original_RARE_image')
        
        fig.add_subplot(244)
        pylab.imshow(roi1.data[:,:,3])
        pylab.title('ROI_from RARE')
        
        fig.add_subplot(245)
        pylab.imshow(roi_resample1[:,:,3])
        pylab.title('Resampled ROI')
        
        fig.add_subplot(246)
        pylab.imshow(curr_T1s1[:,:,3])
        pylab.title('Multiplied LL and resampled ROI')
    
        pylab.savefig(k, bbox_inches=0, dpi=300)
        pylab.close('all')

    except KeyError:
        
        print('This did not work: %s' %k)
        continue


pylab.close('all')


