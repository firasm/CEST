# -*- coding: utf-8 -*-

## NecS3 Analysis File

# Imports

import sarpy
import numpy
import sarpy.fmoosvi.wrappers
import time
import pylab

NecS3_exp = sarpy.Experiment('NecS3')

## AUC 60 Maps
NecS3_DCEscans = [study.find_scan_by_protocol('06')[0] for study in NecS3_exp.studies if len(study.find_scan_by_protocol('06'))>0]

start_time = time.time()
for scan in NecS3_DCEscans:
    
    try:
        AUC = sarpy.fmoosvi.wrappers.calculate_AUC(scan)
        scan.store_adata(key='AUC60', data = AUC)
    except:
        print scan

end_time = time.time()

print('This run took {0} seconds.'.format(round(end_time - start_time)))




# Look-Locker T1 maps
NecS3_LLscans = NecS3_exp.find_scan_by_protocol('04_ubcLL+')

start_time = time.time()

for scan in NecS3_LLscans:
    try:
        T1map_LL, T1_fit_dict = sarpy.fmoosvi.wrappers.calculate_T1map(scan, protocol_name = '04_ubcLL+')
        scan.store_adata(key='T1map_LL', data = T1map_LL)
        scan.store_adata(key='T1_fit_dict', data = T1_fit_dict)

    except IOError:
        print scan
    except:
        continue

end_time = time.time() 
print 'This run took {0} seconds.'.format(round(end_time - start_time))






#BSB1map_MSME = sarpy.fmoosvi.wrappers.calculate_BSB1map(NecS3_exp, protocol_name = 'BSB1map\-MSME_bas')
#BSB1map_FLASH = sarpy.fmoosvi.wrappers.calculate_BSB1map(NecS3_exp, protocol_name = '07_bSB1mapFLASH')
#NecS3_BSB1_FLASH = sarpy.Experiment('NecS3').find_scan_by_protocol('07_bSB1mapFLASH')
#NecS3_BSB1_MSME = sarpy.Experiment('NecS3').find_scan_by_protocol('BSB1map\-MSME_bas')
#    T1map_LL_MSME = sarpy.fmoosvi.wrappers.calculateT1map(NecS3_exp, protocol_name = '04_ubcLL+', flip_angle_map = BSB1map_MSME)
#    T1map_LL_FLASH = sarpy.fmoosvi.wrappers.calculateT1map(NecS3_exp, protocol_name = '04_ubcLL+', flip_angle_map = BSB1map_FLASH)

#    scan.store_adata(key='BSB1map_MSME', data = BSB1map_MSME)
#    scan.store_adata(key='BSB1map_FLASH', data = BSB1map_FLASH)
#    scan.store_adata(key='T1map_LL_MSME', data = T1map_LL_MSME)
#    scan.store_adata(key='T1map_LL_FLASH', data= T1map_LL_FLASH)        


