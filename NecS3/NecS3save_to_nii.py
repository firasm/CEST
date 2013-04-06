# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 05:23:58 2013

@author: firas
"""
import sarpy
import numpy
import sarpy.fmoosvi.wrappers
import json
import nibabel

## Script to save IR-RARE files as niftis


with open('/Volumes/Data/Dropboxes/PhD./Dropbox/Studies/analysis/NecS3/NecS3.json','r') as infile:
    master_sheet = json.load(infile)

NecS3 = sarpy.Experiment('NecS3').find_scan_by_protocol('04')

################################
# This is to export images for ROIs
################################

#for k,v in master_sheet.iteritems():
#    
#    try:
#        name = k+'-1'
#        
#        sarpy.Scan(master_sheet[k]['24h-IR_A'][0]).pdata[0].export2nii(name)
#        
#    except:
#        print k


################################################################
# This is to load them in after ROIs have been drawn using ImageJ
################################################################

for k,v in master_sheet.iteritems():
   
    try: 
        name = k+'-2'+'.nii'
       
        scan = sarpy.Scan(master_sheet[k]['24h-IR_A'][0])
        roi = nibabel.load(name)
        
        scan.store_adata(key='IR_tumour_rois', data = roi)

    except:
        print k