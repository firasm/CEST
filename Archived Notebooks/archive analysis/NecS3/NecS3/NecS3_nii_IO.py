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

#with open('/Users/firas/Desktop/Dropbox/Studies/analysis/NecS3/NecS3.json','r') as infile:
#    master_sheet = json.load(infile)

NecS3 = sarpy.Experiment('NecS3').find_scan_by_protocol('04')

################################
# This is to export images for ROIs
################################

export_path = '/Volumes/Data/Dropboxes/PhD./Dropbox/Studies/analysis/NecS1/roi/input1/'


for k,v in master_sheet.iteritems():
    
    try:
        name = export_path + k+'-2'
        
        sarpy.Scan(master_sheet[k]['24h-IR_A'][0]).pdata[0].export2nii(name)
        
    except:
        print 'Problem with' + k

 #Later realized I also need to get the day 1 images for LL comparisons

for k,v in master_sheet.iteritems():
    
    try:
        name = export_path + k+'-1'
        
        sarpy.Scan(master_sheet[k]['0h-IR_A'][0]).pdata[0].export2nii(name)
        
    except:
        print 'Problem with' + k


################################################################
# This is to load them in after ROIs have been drawn using ImageJ
################################################################
#
#import_path = '/Volumes/Data/Dropboxes/PhD./Dropbox/Studies/analysis/NecS3/roi/output1/'
#
#for k,v in master_sheet.iteritems():
#   
#    try: 
#        name1 = import_path + k+'-1'+'.nii'
#        name2 = import_path + k+'-2'+'.nii'
#
#        scan1 = sarpy.Scan(master_sheet[k]['0h-IR_A'][0])
#        scan2 = sarpy.Scan(master_sheet[k]['24h-IR_A'][0])              
#        
#        roi1 = nibabel.load(name1)
#        roi2 = nibabel.load(name2)        
#        
#        scan1.store_adata(key='IR_tumour_rois', data = roi1.get_data()[:,:,:,0], force = True)
#        scan2.store_adata(key='IR_tumour_rois', data = roi2.get_data()[:,:,:,0], force = True)
#
#    except:
#        print('Error saving {0}'.format(k))
#


