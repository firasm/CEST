# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 21:19:18 2013

@author: fmoosvi
"""
import sarpy 
import json
import collections
import sarpy.fmoosvi.getters

def magic_tuple(scn):
    return (scn.shortdirname, scn.acqp.ACQ_protocol_name)  
    

pat_list = sarpy.fmoosvi.getters.get_patients_from_experiment('NecS3', verbose = True)

master_sheet = collections.OrderedDict()
for pat in pat_list:
    nm = pat.get_SUBJECT_id()
    assert len(set(nm)) ==1
    master_sheet[nm[0]]={}
    
    if pat.get_SUBJECT_id()[0] == 'NecS3Hs12':
        stdy1= pat.studies[0]
        stdy2= pat.studies[3]
    elif pat.get_SUBJECT_id()[0] == 'NecS3Hs04':
        stdy1= pat.studies[0]
        stdy2= pat.studies[2]
    elif len(pat.studies) != 2:
        print 'not 2 studies: %s ' % pat.get_SUBJECT_id()
        continue
    else:
        stdy1 = pat.studies[0]
        stdy2 = pat.studies[1]

    master_sheet[nm[0]]['0h'] = stdy1.shortdirname
    master_sheet[nm[0]]['24h'] = stdy2.shortdirname
    
    #LL
    scn = stdy1.find_scan_by_protocol('04_')
    if len(scn) > 1:
        print 'more than one LL on day 1 %s' % len(scn) + scn[0].shortdirname
    master_sheet[nm[0]]['0h-LL']= magic_tuple(scn[0])
    scn = stdy2.find_scan_by_protocol('04_')
    if len(scn) > 1:
        print 'more than one LL on day 2 %s' % len(scn) + scn[0].shortdirname
    master_sheet[nm[0]]['24h-LL']= magic_tuple(scn[0])

    #DCE 1 & 2         
    scn = stdy1.find_scan_by_protocol('06_')
    if len(scn) >= 2:
        master_sheet[nm[0]]['0h-DCE1']= magic_tuple(scn[-2])
        master_sheet[nm[0]]['0h-DCE2']= magic_tuple(scn[-1])
    elif len(scn) == 1:
        print '1 DCE on day 1 %s ' % stdy1.shortdirname
        
    scn = stdy2.find_scan_by_protocol('06_')
    if len(scn) >= 2:
        master_sheet[nm[0]]['24h-DCE1']= magic_tuple(scn[-2])
        master_sheet[nm[0]]['24h-DCE2']= magic_tuple(scn[-1])
    elif len(scn) ==1:
        print '1 DCE on day 2 %s ' % stdy1.shortdirname
      
     # RARE
    scn = stdy1.find_scan_by_protocol('05_.*RARE')
    if len(scn) == 2:
        master_sheet[nm[0]]['0h-IR_A']= magic_tuple(scn[0])
        master_sheet[nm[0]]['0h-IR_B']= magic_tuple(scn[1])
    elif len(scn) == 1:
        master_sheet[nm[0]]['0h-IR_A']= magic_tuple(scn[0])
    else:
        print 'RARE not clear, day1 %s ' % stdy1.shortdirname
        
    scn = stdy2.find_scan_by_protocol('05_.*RARE')
    if len(scn) == 2:
        master_sheet[nm[0]]['24h-IR_A']= magic_tuple(scn[0])
        master_sheet[nm[0]]['24h-IR_B']= magic_tuple(scn[1])
    elif len(scn) == 1:
        master_sheet[nm[0]]['24h-IR_A']= magic_tuple(scn[0])
    else:
        print 'RARE not clear, day2 %s ' % stdy2.shortdirname
    
with open('/Volumes/Data/Dropboxes/PhD./Dropbox/Studies/NecS3/NecS3.json','wb') as outfile:
    json.dump(master_sheet, outfile, indent=4)

#    with open('/Users/fmoosvi/NecS3.json','r') as infile:
#        x = json.load(infile)
    
    
#    for k,v in x.iteritems():
#        try:
#            print k, x[k]['24h-DCE2'][0]
#        except:
#            print k, 'too bad'
    
    
## After running  the summary thing:

#==============================================================================
# No LL for NecS3Ao 

# No AUC60 for NecS3Hs04

# No IR for NecS3Hs10 - manually fixed    
# No IR for NecS3Hs02 - manually fixed 
#==============================================================================
