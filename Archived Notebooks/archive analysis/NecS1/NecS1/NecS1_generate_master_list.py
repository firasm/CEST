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
    
pat_list = sarpy.fmoosvi.getters.get_patients_from_experiment('NecS1')
#pat_list = sarpy.fmoosvi.getters.get_patients_from_experiment('NecS1', verbose = True)

master_sheet = collections.OrderedDict()
for pat in pat_list:
    nm = pat.get_SUBJECT_id()
    assert len(set(nm)) ==1
    master_sheet[nm[0]]={}
    
# NecS1Hs01 - Day 2 imaging session not present 
          
    if pat.get_SUBJECT_id()[0] == 'NecS1Hs01':
        continue
    elif pat.get_SUBJECT_id()[0] == 'NecS1Hs04':
        continue
    elif pat.get_SUBJECT_id()[0] == 'NecS1Hs06':
        # Continuing because this mouse was poorly positiond galore!
        # the T1 maps also look like crap, likely due to poor snr. investigate
        #TODO: Investgate wat happened here
        continue
    
    elif pat.get_SUBJECT_id()[0] == 'NecS1Hs07':
        stdy1 = pat.studies[0]
        stdy2 = pat.studies[2]
        stdy3 = pat.studies[2] 
        print('Warning,the NecS3Hs07 does NOT have a third day and is filled with day 2 study!')
    
    elif len(pat.studies) != 3:
        print 'not 2 studies: %s ' % pat.get_SUBJECT_id()
    else:
        stdy1 = pat.studies[0]
        stdy2 = pat.studies[1]
        stdy3 = pat.studies[2]

    master_sheet[nm[0]]['0h'] = stdy1.shortdirname
    master_sheet[nm[0]]['24h'] = stdy2.shortdirname
    master_sheet[nm[0]]['48h'] = stdy3.shortdirname


    #LL
    scn = stdy1.find_scan_by_protocol('04_')
    if len(scn) != 1:
        if len(scn) >1:
            print 'LL weird on day 1 %s' % len(scn) + scn[0].shortdirname
        else:
            print 'LL weird on day 1 %s' % len(scn)
    else:
         master_sheet[nm[0]]['0h-LL']= magic_tuple(scn[0])
    
    scn = stdy2.find_scan_by_protocol('04_')
    if len(scn) != 1:
        if len(scn) >1:
            print 'LL weird on day 2 %s' % len(scn) + scn[0].shortdirname
        else:
            print 'LL weird on day 2 %s' % len(scn)
            stx = stdy2
    else: 
         master_sheet[nm[0]]['24h-LL']= magic_tuple(scn[0])

    scn = stdy3.find_scan_by_protocol('04_')
    if len(scn) != 1:
        if len(scn) >1:
            print 'LL weird on day 3 %s' % len(scn) + scn[0].shortdirname
        else:
            print 'LL weird on day 3 %s' % len(scn)
    else:
         master_sheet[nm[0]]['48h-LL']= magic_tuple(scn[0])
    
    
    #DCE 1 & 2     
    scn = stdy1.find_scan_by_protocol('06_')
    if len(scn) == 2:
        master_sheet[nm[0]]['0h-DCE1']= magic_tuple(scn[-2])
        master_sheet[nm[0]]['0h-DCE2']= magic_tuple(scn[-1])
    elif len(scn) == 1:
        print '1 DCE on day 1 %s ' % stdy1.shortdirname
    else:
        print 'Funky stuff going on here on Day 1 %s' % stdy1.shortdirname
        
    scn = stdy2.find_scan_by_protocol('06_')
    if len(scn) == 2:
        master_sheet[nm[0]]['24h-DCE1']= magic_tuple(scn[-2])
        master_sheet[nm[0]]['24h-DCE2']= magic_tuple(scn[-1])
    elif len(scn) ==1:
        print '1 DCE on day 2 %s ' % stdy2.shortdirname
    else:
        print 'Funky stuff going on here on Day 2 %s' % stdy2.shortdirname      
        
    scn = stdy3.find_scan_by_protocol('06_')
    if len(scn) == 2:
        master_sheet[nm[0]]['48h-DCE1']= magic_tuple(scn[-2])
        master_sheet[nm[0]]['48h-DCE2']= magic_tuple(scn[-1])
    elif len(scn) ==1:
        print '1 DCE on day 3 %s ' % stdy3.shortdirname        
    elif len(scn) ==0:
        print 'No Day 3 DCE for %s' % stdy3.shortdirname        
      
     # RARE
    scn = stdy1.find_scan_by_protocol('05_.*RARE')
    if len(scn) == 2:
        master_sheet[nm[0]]['0h-IR_A']= magic_tuple(scn[0])
        master_sheet[nm[0]]['0h-IR_B']= magic_tuple(scn[1])
    elif len(scn) == 1:
        master_sheet[nm[0]]['0h-IR_A']= magic_tuple(scn[0])        
    else:
        print 'RARE not clear, day1. Trying scan 7 %s ' % stdy1.shortdirname
        try:
            scn = stdy1.find_scan_by_protocol('07_.*RARE*')
            if len(scn) == 2:
                master_sheet[nm[0]]['0h-IR_A']= magic_tuple(scn[0])
                master_sheet[nm[0]]['0h-IR_B']= magic_tuple(scn[1])
                print '\t Try Suceeded, warning: different scan'
            else:
                raise IOError
        except IOError:
            print '\t Try failed, study likely not useful'
        
    scn = stdy2.find_scan_by_protocol('05_.*RARE')
    if len(scn) == 2:
        master_sheet[nm[0]]['24h-IR_A']= magic_tuple(scn[0])
        master_sheet[nm[0]]['24h-IR_B']= magic_tuple(scn[1])
    elif len(scn) == 1:
        master_sheet[nm[0]]['24h-IR_A']= magic_tuple(scn[0])        
    else:
        print 'RARE not clear, day2. Trying scan 7 %s ' % stdy2.shortdirname
        try:
            scn = stdy2.find_scan_by_protocol('07_.*RARE*')
            if len(scn) == 2:
                master_sheet[nm[0]]['24h-IR_A']= magic_tuple(scn[0])
                master_sheet[nm[0]]['24h-IR_B']= magic_tuple(scn[1])
                print '\t Try Suceeded, warning: different scan'
            else:
                raise IOError
        except IOError:
            print '\t Try failed, study likely not useful'
            

                
    scn = stdy3.find_scan_by_protocol('05_.*RARE')
    if len(scn) == 2:
        master_sheet[nm[0]]['48h-IR_A']= magic_tuple(scn[0])
        master_sheet[nm[0]]['48h-IR_B']= magic_tuple(scn[1])
    elif len(scn) == 1:
        master_sheet[nm[0]]['48h-IR_A']= magic_tuple(scn[0])
    else:
        print 'RARE not clear, day3. Trying scan 7 %s ' % stdy1.shortdirname
        try:
            scn = stdy3.find_scan_by_protocol('07_.*RARE*')
            if len(scn) == 2:
                master_sheet[nm[0]]['48h-IR_A']= magic_tuple(scn[0])
                master_sheet[nm[0]]['48h-IR_B']= magic_tuple(scn[1])
                print '\t Try Suceeded, warning: different scan'
            else:
                stx = stdy3
                raise IOError
        except IOError:
            print '\t Try failed, study likely not useful'
        
        
        
with open('/Volumes/Data/Dropboxes/PhD./Dropbox/studies/analysis/NecS1/NecS1.json','wb') as outfile:
    json.dump(master_sheet, outfile, indent=4)

with open('/Volumes/Data/Dropboxes/PhD./Dropbox/studies/analysis/NecS1/NecS1.json','r') as infile:
        x = json.load(infile)
       