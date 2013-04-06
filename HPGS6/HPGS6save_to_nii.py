# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 14:16:31 2013

@author: fmoosvi
"""
import sarpy

HPGS6 = sarpy.Experiment('HPGS6')
HPSS6 = sarpy.Experiment('HPSS6')

for study in HPSS6.studies:
    HPGS6.add_study(study)

HPGS6_IR = []
HPGS6_IR.append(HPGS6.studies[0].scans[-4])
HPGS6_IR.append(HPGS6.studies[1].scans[-3])
HPGS6_IR.append(HPSS6.studies[0].scans[-7])
HPGS6_IR.append(HPSS6.studies[1].scans[-8])


names = ['HPGS6-ix1',
         'HPGS6-iy1',
         'HPGS6Ss03',
         'HPGS6Ss04']

for i in xrange(len(HPGS6_IR)):
    
    scan = HPGS6_IR[i]
    scan.pdata[0].export2nii(names[i])

