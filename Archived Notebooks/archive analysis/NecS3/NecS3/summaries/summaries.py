# -*- coding: utf-8 -*-

# Testing on T1 maps from NecS3

import sarpy
import sarpy.fmoosvi.getters
import sarpy.fmoosvi.wrappers
import pylab
import numpy
import json

with open('/Volumes/Data/Dropboxes/PhD./Dropbox/studies/analysis/NecS3/NecS3.json','r') as infile:
    master_sheet = json.load(infile)
    
##################### LL T1 Maps ##############################################  
data_list = []
for k,v in master_sheet.iteritems():

    try:
        data_list = []
        key_list = []
        key_list.append(k)
        key_list.append('T1map-')
        key_list.append('0h-LL')
        key_list.append('24h-LL')
        key_list.append('T1map_LL')

        
        data1 = sarpy.Scan(master_sheet[k][key_list[2]][0]).adata[key_list[4]].data
        data2 = sarpy.Scan(master_sheet[k][key_list[3]][0]).adata[key_list[4]].data
        
        data_list.append(data1)
        data_list.append(data2)
        
        sarpy.fmoosvi.wrappers.create_summary(data_list, key_list, clims = [600, 3500] )
        
    except:
        print 'No LL for', k

##################### AUC60 Maps ##############################################       
for k,v in master_sheet.iteritems():

    try:
        try:
            data_list = []
            key_list = []
            key_list.append(k)
            key_list.append('AUC60-')
            key_list.append('0h-DCE1')
            key_list.append('24h-DCE1')
            key_list.append('AUC60')
    
            
            data1 = sarpy.Scan(master_sheet[k][key_list[2]][0]).adata[key_list[4]].data
            data2 = sarpy.Scan(master_sheet[k][key_list[3]][0]).adata[key_list[4]].data
            
            data_list.append(data1)
            data_list.append(data2)
            
            sarpy.fmoosvi.wrappers.create_summary(data_list, key_list,clims = [10,60])            
        except:
            
            data_list = []
            key_list = []
            key_list.append(k)
            key_list.append('AUC60-')
            key_list.append('24h-DCE1')
            key_list.append('AUC60')
    
            
            data2 = sarpy.Scan(master_sheet[k][key_list[2]][0]).adata[key_list[3]].data
            data1 = 0*data2
            
            data_list.append(data1)
            data_list.append(data2)
            
            sarpy.fmoosvi.wrappers.create_summary(data_list, key_list,clims = [10,60])
 
    except:
        print 'No AUC60 for',k

#################### IR-RARE ##############################################       
for k,v in master_sheet.iteritems():

    try:
        try:
            data_list = []
            key_list = []
            key_list.append(k)
            key_list.append('IR-')
            key_list.append('0h-IR_A')
            key_list.append('0h-IR_B')
            key_list.append('24h-IR_A')
            key_list.append('24h-IR_B')            
    
            
            data1 = sarpy.Scan(master_sheet[k][key_list[2]][0]).pdata[0].data
            data2 = sarpy.Scan(master_sheet[k][key_list[3]][0]).pdata[0].data
            data3 = sarpy.Scan(master_sheet[k][key_list[4]][0]).pdata[0].data
            data4 = sarpy.Scan(master_sheet[k][key_list[5]][0]).pdata[0].data
            
            data_list.append(data1)
            data_list.append(data2)
            data_list.append(data3)
            data_list.append(data4)
            
            lim = sarpy.fmoosvi.getters.get_image_clims(data1)
            
            sarpy.fmoosvi.wrappers.create_summary(data_list, key_list, clims = lim, colour_map = 'gray')
        except:
            data_list = []
            key_list = []
            key_list.append(k)
            key_list.append('IR-')
            key_list.append('0h-IR_A')
            key_list.append('24h-IR_A')
            key_list.append('24h-IR_B')            
                
            data1 = sarpy.Scan(master_sheet[k][key_list[2]][0]).pdata[0].data
            data2 = data1*1e9
            data3 = sarpy.Scan(master_sheet[k][key_list[3]][0]).pdata[0].data
            data4 = sarpy.Scan(master_sheet[k][key_list[4]][0]).pdata[0].data
            
            data_list.append(data1)
            data_list.append(data2)
            data_list.append(data3)
            data_list.append(data4)
            
            lim = sarpy.fmoosvi.getters.get_image_clims(data1)
            
            sarpy.fmoosvi.wrappers.create_summary(data_list, key_list, clims = lim, colour_map = 'gray')            
            

    except:
        print 'No IR for',k
