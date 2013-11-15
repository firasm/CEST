# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 17:28:17 2013

@author: fmoosvi
"""

import sarpy
import sarpy.fmoosvi.wrappers
from datetime import datetime
from collections import defaultdict

startTime = datetime.now()

masterlist = 'HPGS6/HPGS6'

analyze_list = []

fFlag = False

analyze_list.append(['dce1','auc60',fFlag])
analyze_list.append(['dce2','auc60',fFlag])

analyze_list.append(['LL_before','T1map_LL',fFlag])
analyze_list.append(['LL_after','T1map_LL',fFlag])

#analyze_list.append(['gd_conc'])
#
#analyze_list.append(['dce1','augc60',fFlag])
#analyze_list.append(['dce2','augc60',fFlag])

analyze_list.append(['dce1','vtc',False])
analyze_list.append(['dce2','vtc',False])

for items in analyze_list:
    
#    if items[0] == 'gd_conc':
#        
#        sarpy.fmoosvi.wrappers.conc_from_signal(masterlist, 'dce1', 
#                                        data_label_T1map = '0h-LL',
#                                        adata_label = 'T1map_VFA')                
#        continue
#                                           
    if len(items) == 3:
                  
        sarpy.fmoosvi.wrappers.bulk_analyze(masterlist, items[0], items[1], forceVal = items[2])

    else:
        sarpy.fmoosvi.wrappers.bulk_analyze(masterlist, items[0], items[1], items[2], forceVal = items[3])
