# -*- coding: utf-8 -*-

## NecS3 Analysis File

# Imports

import sarpy
import numpy
import sarpy.fmoosvi.wrappers
import time
import pylab

NecS3_exp = sarpy.Experiment('NecS3')

### AUC 60 Maps
#NecS3_DCEscans = [study.find_scan_by_protocol('06')[0] for study in NecS3_exp.studies if len(study.find_scan_by_protocol('06'))>0]
#
#start_time = time.time()
#for scan in NecS3_DCEscans:
#    
#    try:
#        AUC = sarpy.fmoosvi.wrappers.calculate_AUC(scan)
#        scan.store_adata(key='AUC60', data = AUC, force = True)
#    except:
#        print scan
#
#end_time = time.time()
#
#print('This run took {0} seconds.'.format(round(end_time - start_time)))


# Look-Locker T1 maps
NecS3_LLscans = NecS3_exp.find_scan_by_protocol('04_ubcLL+')

start_time = time.time()

for scan in NecS3_LLscans:
    try:
        T1map_LL, T1_fit_dict = sarpy.fmoosvi.wrappers.calculate_T1map(scan, protocol_name = '04_ubcLL+')
                
        scan.store_adata(key='T1map_LL', data = T1map_LL, force = True)
        scan.store_adata(key='T1_fit_dict', data = T1_fit_dict, force = True)

    except IOError:
        print scan
    except:
        continue

end_time = time.time() 
print 'This run took {0} seconds.'.format(round(end_time - start_time))

