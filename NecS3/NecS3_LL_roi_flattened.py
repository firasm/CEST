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

        curr_T1s1 = LL_1.data * sarpy.fmoosvi.analysis.h_image_to_mask(roi_resample1)
        curr_T1s2 = LL_2.data * sarpy.fmoosvi.analysis.h_image_to_mask(roi_resample1)

        T1_vals_day1.append( curr_T1s1[numpy.isfinite(curr_T1s1)] )
        T1_vals_day2.append( curr_T1s2[numpy.isfinite(curr_T1s2)] )
#        
#        if k == 'NecS3Hs11':  #debuggimg please ignore 
#            pylab.figure()
#            bins = numpy.linspace(200,3500,50)
#            pylab.hist(T1_vals_day2[-1], bins, alpha=0.6)
#            continue
    
        
#        print k
#        pylab.figure()
#        bins = numpy.linspace(200,3500,50)
#        pylab.hist(T1_vals_day2[-1], bins, alpha=0.6)
#        pylab.title(k)

        
    except KeyError:
        #pylab.close('all')
        print('Expected error, please ignore {0}'.format(k))
        
# Now squash the lists together, and turn into an array
# use an incomprehensible list comprehension, as named by someone on SE
T1_vals_day1 = [item for sublist in T1_vals_day1 for item in sublist]
T1_vals_day1 = numpy.array(T1_vals_day1)

T1_vals_day2 = [item for sublist in T1_vals_day2 for item in sublist]
T1_vals_day2 = numpy.array(T1_vals_day2)


        
bins = numpy.linspace(200,3500,50)
pylab.hist(T1_vals_day1, bins, alpha=0.6)
pylab.hist(T1_vals_day2, bins, alpha=0.6)

#bins = numpy.linspace(200,3500,50)
#        
#pylab.figure()        
#pylab.hist(T1_vals_day1, bins, alpha=0.6)
#pylab.hist(T1_vals_day2, bins, alpha=0.6)
#        
#        
#        
#b=pylab.figure()
#a=pylab.imshow(necs3hs10.adata['T1map_LL'].data[:,:,3])
#a.set_clim([500,3500])
#pylab.colorbar()
#b.canvas.draw()
#        
