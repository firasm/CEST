# -*- coding: utf-8 -*-
"""
Created on Sun Apr 14 15:17:12 2013

@author: fmoosvi
"""


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
        key_list = []
        key_list.append(k) # 0
        key_list.append('NecIndex-') #1
        key_list.append('0h-IR_A') #2
        key_list.append('24h-IR_A') #3
        key_list.append('24h-IR_B') #4
#        key_list.append('IR_tumour_rois') #5
        
        data_list = []
        
        scanB = sarpy.Scan(master_sheet[k][key_list[3]][0]).pdata[0].data
        scanC = sarpy.Scan(master_sheet[k][key_list[4]][0]).pdata[0].data
        
        data_list.append(scanB)
        data_list.append(scanC)
        
        clims = sarpy.fmoosvi.getters.get_image_clims(scanB)

        sarpy.fmoosvi.wrappers.create_summary(data_list,key_list, clims)
        
    except KeyError:
        
        print('Key not found for {0}'.format(k))