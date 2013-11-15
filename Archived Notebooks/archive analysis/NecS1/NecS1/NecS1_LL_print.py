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

with open('/Volumes/Data/Dropboxes/PhD./Dropbox/studies/analysis/NecS1/NecS1.json','r') as infile:
    master_sheet = json.load(infile)

for k,v in master_sheet.iteritems():
    
    try:
        key_list = []
        key_list.append(k) # 0
        key_list.append('Nec-') #1
        key_list.append('24h-LL') #2

        # Get the LL data
        LL_1 = sarpy.Scan(master_sheet[k][key_list[2]][0])
        
        for slice in xrange(6):
            fig = pylab.figure(figsize = (20,20))
            G = pylab.matplotlib.gridspec.GridSpec(5,5)   
            
            for i in xrange(25):
                
                data = LL_1.pdata[0].data[:,:,slice,i]
                fig.add_subplot(G[int(numpy.floor(i/5)),i%5],frameon=False, xticks=[], yticks=[])
                clims = sarpy.fmoosvi.getters.get_image_clims(data)

                a = pylab.imshow(data, cmap = 'jet' )
                a.set_clim(clims[0],clims[1])
                
                fig.canvas.draw()

            pylab.suptitle('Slice {0}'.format(slice+1), fontsize = 14)

            # Figure spacing adjustments
            #fig.subplots_adjust(right = 0.85, wspace = 0.0001, hspace=0.0001)
            G.tight_layout(fig, h_pad = 0.1, w_pad = 0.001)
            G.update(right = 0.87)
            
            # Saving Figure    
            filename = key_list[1] + key_list[0] + '-sl-' + str(slice+1) + '.png'                
            pylab.savefig(filename, bbox_inches=0, dpi=300)
            pylab.close('all')        
            
        # Print the number of pixels used in each study
        print k
        
    except KeyError:
        print('Key error, in {0}'.format(k))

