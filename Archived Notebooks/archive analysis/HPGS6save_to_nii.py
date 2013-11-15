# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 14:16:31 2013

@author: fmoosvi
"""
import sarpy
import nibabel
import os
import sarpy.fmoosvi.wrappers
import sarpy.fmoosvi.getters
import sarpy.ImageProcessing.resample_onto
import pylab
import numpy

HPGS6 = sarpy.Experiment('HPGS6')
HPSS6 = sarpy.Experiment('HPSS6')

for study in HPSS6.studies:
    HPGS6.add_study(study)

##HPGS6_IR = []
##HPGS6_IR.append(HPGS6.studies[0].scans[-4])
##HPGS6_IR.append(HPGS6.studies[1].scans[-3])
##HPGS6_IR.append(HPSS6.studies[0].scans[-7])
##HPGS6_IR.append(HPSS6.studies[1].scans[-8])
##
##
##names = ['HPGS6-ix1',
##         'HPGS6-iy1',
##         'HPGS6Ss03',
##         'HPGS6Ss04']
##
##for i in xrange(len(HPGS6_IR)):
##    
##    scan = HPGS6_IR[i]
##    scan.pdata[0].export2nii(names[i])
#
#print ('Most of this will stop working due to the changes on resample_onto, be warned!!')
#
#
##### Loading in the ROIs
#
#scan1 = sarpy.Scan("HPSS6Ss03.iy1/8")
#n1 = os.path.join('/Volumes','Data','Dropboxes','PhD.', 'Dropbox','Studies','analysis','HPGS6','roi','output1','HPGS6Ss03.nii')
#
#roi1 = nibabel.load(n1)
#scan1.store_adata(key='IR_tumour_rois', data = roi1.get_data()[:,:,:,0],force = True)
#
#scan2 = sarpy.Scan("HPSS6Ss04.iy1/11")
#n2 = os.path.join('/Volumes','Data','Dropboxes','PhD.', 'Dropbox','Studies','analysis','HPGS6','roi','output1','HPGS6Ss04.nii')
#
#roi2 = nibabel.load(n2)
#scan2.store_adata(key='IR_tumour_rois', data = roi2.get_data()[:,:,:,0], force = True)
#
#
##### Get the T1 distributions by ROI
#
#HPGS6_LLscans = HPGS6.find_scan_by_protocol('04_ubcLL2')
#
#names = ['T1distrib-HPSS6Ss03a', 'T1distrib-HPSS6Ss03b', 'T1distrib-HPGS6Ss04-after']
#rois = [roi1, roi1, roi2]
#counter = 0
#
#for scan in HPGS6_LLscans[-3:]:
#    try:
#        #T1map_LL, T1_fit_dict = sarpy.fmoosvi.wrappers.calculate_T1map(scan, protocol_name = '04_ubcLL+')
#              
##        scan.store_adata(key='T1map_LL', data = T1map_LL)
##        scan.store_adata(key='T1_fit_dict', data = T1_fit_dict)
#        
#        scan.adata['T1map_LL'].export2nii('LL1.nii')
#       
#        LLresample1 = sarpy.ImageProcessing.resample_onto.resample_onto('LL1.nii',rois[counter].get_filename())      
#        
#        name = names[counter]
#        
#        bins = numpy.linspace(500,3500,20)
#        roi = numpy.array(rois[counter].get_data()[:,:,:,0]).copy()
#
#        sarpy.fmoosvi.wrappers.roi_distribution(LLresample1, roi, bins, save_histogram = True, save_name = name)
#        
#        pylab.close('all')
#        counter = counter+1
#        
#    except:
#        
#        print scan
#        print("You're an idiot")
#        
######## Exporting T1 maps #########        
#        
#counter = 0        
#names = ['HPSS6Ss03-1', 'HPSS6Ss03-2', 'HPGS6Ss04']
#        
#data_list = []
#data1 = HPGS6_LLscans[-3].adata['T1map_LL'].data
#data2 = HPGS6_LLscans[-2].adata['T1map_LL'].data
#data3 = numpy.zeros(shape = data2.shape)
#data4 = HPGS6_LLscans[-1].adata['T1map_LL'].data
#            
#data_list.append(data1)
#data_list.append(data2)
#data_list.append(data3)
#data_list.append(data4)
#
#key_list = []
#key_list.append('HPGS6-T1Summary')
#key_list.append('T1map-')
#key_list.append('T1map_LL')
#
#sarpy.fmoosvi.wrappers.create_summary(data_list, key_list,clims = [500,3000])  
#    
######## Exporting AUC maps #########        
#
#scan1 = sarpy.Scan("HPSS6Ss03.iy1/9")
#scan2 = sarpy.Scan("HPSS6Ss04.iy1/12")
#
#AUC1 = sarpy.fmoosvi.wrappers.calculate_AUC(scan1)
#scan1.store_adata(key='AUC60', data = AUC1, force = True)
#
#AUC2 = sarpy.fmoosvi.wrappers.calculate_AUC(scan2)
#scan2.store_adata(key='AUC60', data = AUC2, force = True)
#
#key_list = []
#data_list = [scan1.adata['AUC60'].data]
#
#key_list.append('HPSS6Ss03')
#key_list.append('AUC60-')
#
#sarpy.fmoosvi.wrappers.create_summary(data_list, key_list,clims = [0,25])  
#
#key_list = []
#data_list = [scan2.adata['AUC60'].data]
#
#key_list.append('HPSS6Ss04')
#key_list.append('AUC60-')
#
#sarpy.fmoosvi.wrappers.create_summary(data_list, key_list,clims = [0,25])  

####### Exporting DCE-MRI Enhancement curves #########        

scan1 = sarpy.Scan("HPSS6Ss03.iy1/9")
scan1_roi = sarpy.Scan("HPSS6Ss03.iy1/8").adata['IR_tumour_rois']

ec_scan1 = sarpy.fmoosvi.getters.get_enhancement_curve(scan1, scan1_roi)

scan2 = sarpy.Scan("HPSS6Ss04.iy1/12")
scan2_roi = sarpy.Scan("HPSS6Ss04.iy1/11").adata['IR_tumour_rois']

ec_scan2 = sarpy.fmoosvi.getters.get_enhancement_curve(scan2, scan2_roi)

key_list = []
data_list = [ec_scan1]

key_list.append('HPSS6Ss03')
key_list.append('EC-')

sarpy.fmoosvi.wrappers.create_plot(data_list, key_list)  

key_list = []
data_list = [ec_scan2]

key_list.append('HPSS6Ss04')
key_list.append('EC-')    
    
sarpy.fmoosvi.wrappers.create_plot(data_list, key_list)  

    
    