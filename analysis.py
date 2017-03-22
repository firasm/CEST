# CEST Analysis Code

####### Imports #######

import numpy
import sarpy
import sarpy.analysis.analysis
import scipy
import collections
import pylab
import random

####### Fitting Functions #######

def h_zspectrum_N(params,freqs):
    
    arr = numpy.empty_like(freqs)*0
    shift =  params[0]
    
    for i in numpy.arange(0,len(params[1:]),3):
        
        A = params[i+1] # Amplitude of the peak
        lw =  params[i+2] # Line width (in units of freqs)
        w0 = params[i+3] # Resonant frequency of sample        

        # Construct the multiple lorentzians added together
        tmp = numpy.divide(A,(1+4*((freqs-w0)/lw)**2))
        
        arr += tmp
    return (arr+shift)

def h_zspectrum_New(paramStructuredArray,freqs):
    ''' Updated Zspectrum N function to now require a structured array of params'''
    
    arr = numpy.zeros_like(freqs)
    shift =  paramStructuredArray['offset']

    # Convert this to a regular list so that it can be added as lorentzians
    params = []
    for f in paramStructuredArray.dtype.fields:
        params.append(paramStructuredArray[f])
        
    params = numpy.array(params)
    
    for i in numpy.arange(0,len(params[1:]),3):
        
        A = params[i+1] # Amplitude of the peak
        lw =  params[i+2] # Line width (in units of freqs)
        w0 = params[i+3] # Resonant frequency of sample        

        # Construct the multiple lorentzians added together
        tmp = numpy.divide(A,(1+4*((freqs-w0)/lw)**2))
        
        arr = arr+tmp
    return (arr+shift)

def h_residual_Zspectrum_N(params, y_data, w):
    
    return numpy.abs(y_data - h_zspectrum_N(params,w))    

def fit_water_peak(data,offset_freqs,allParams = False):

    """
    Fits a Lorentzian to the data, and 
    returns the water offset frequency

    """

    # First get rid of all attempts to pass in bad data. If there are any nans in there, toss it.

    if numpy.isnan(numpy.sum(data)):
        return numpy.nan

    else:

        # Presumably the maximum signal will be at 
        # some large frequency offset, so let's just
        # use that as our baseline parameter

        # Also, normalize the data so the fit is easier

        params_passed = [numpy.max(data),-1.,0.6,offset_freqs[numpy.argmin(data)]]

        fit_params,cov,infodict,mesg,ier = scipy.optimize.leastsq(h_residual_Zspectrum_N,
            params_passed,
            args=(data, offset_freqs),
            full_output = True,
            maxfev = 200)
        
        if allParams:
            return fit_params
        else:
            return fit_params[3]

def shift_water_peak(scn_to_analyse=None, 
                     bbox = None,
                     pdata_num = 0,
                     **kwargs):

# DELETE; likely not needed anymore March 15, 2017
    
    """
    Returns a new array that should replace 
    scn.method.CEST_FreqListPPM as the x-axis
    for the CEST-spectrum.
    
    #TODO Currently only works for Single slice data

    :param str scn_to_analyse: scan.shortdirname
    :param integer pdata_num: reconstruction number, according to python numbering.
            default reconstruction is pdata_num = 0.
    :param bbox: bounding bbox
    :return: array with shifted peak positions
    """ 

    scn = sarpy.Scan(scn_to_analyse)

    # Size info
    x_size = scn.pdata[pdata_num].data.shape[0]
    y_size = scn.pdata[pdata_num].data.shape[1]  

    # Offsets
    freq_list = scn.method.CEST_FreqListPPM
    offsets = len(freq_list)

    #### Deal with bounding boxes

    if bbox is None:        
        bbox = numpy.array([0,x_size-1,0,y_size-1])    
    else:
        bbox = scn.adata['bbox'].data

    # Fit the water peak pixel by pixel
    water_shift_map = numpy.nan + numpy.empty([x_size,y_size])

    dat = scn.pdata[0].data
    for x in range(bbox[0],bbox[1]):
        for y in range(bbox[2],bbox[3]):

            fit_data = dat[x,y] / numpy.max(dat[x,y])

            water_shift_map[x,y] = fit_water_peak(fit_data,freq_list)

    return water_shift_map

####### Displaying and Plotting Functions #######    

def cest_spectrum(scn_to_analyse,
                  xval,
                  yval,
                  shift_water_peak = True,
                  normalize=True,
                  normalize_to_ppm = 200,
                  ppm_limit_min = -50,
                  ppm_limit_max = 50,
                  exclude_ppm = 66.6,
                  pdata_num = 0):
    
    '''
    Does the grunt work of processing freqs to give
    a sequential list of frequencies (rather than alternating).

    :param str scn_to_analyse: scan.shortdirname
    :param bool shift_water_peak: Interpolates intensities so the water peak is at 0 ppm (you can now average after this)
    :param float normalize_to_ppm: freq to normalize to
    :param float ppm_limit: Only return frequencies within +ppm_limit and -ppm_limit
    :param float exclude_ppm: Exclude the dummy frequencies
    :param int pdata_num: reconstruction number, according to python numbering.
            default reconstruction is pdata_num = 0.

    :return: new (frequencies in sequential order)
             tmp (data at frequencies in 'new')
    '''
    scn = sarpy.Scan(scn_to_analyse)

    # Get the Frequencies and make them sequential 
    # so that alternating points don't plot badly

    # First get the freq list
    freq_list = scn.method.CEST_FreqListPPM

    if normalize:
    # Find the frequency to normalize to, throw error if not found
        possibleNormalizations = [i for i, x in enumerate(freq_list) if numpy.abs(x - normalize_to_ppm) <1E-4]

        # By default, select the LAST instance of the freq to avoid approach to steady state issues
        normalizeTo = possibleNormalizations[-1]

    # Get only the frequencies within the ppm_limit
    new = [f for f in freq_list if f > ppm_limit_min]
    new = [f for f in new if f < ppm_limit_max]

    # Exclude the dummy frequencies at the beginning (66.6 ppm)
    new = sorted([n for n in new if n!= exclude_ppm])

    # Get the index of the good frequencies relative to the original list 
    ind = [freq_list.index(c) for c in new]  

    if normalize:
        # Get the data and normalize it to index of normalize_to_ppm
        tmp = scn.pdata[pdata_num].data[xval,yval,:] / scn.pdata[0].data[xval,yval,normalizeTo] 

    else:

        tmp = scn.pdata[pdata_num].data[xval,yval,:]

    if shift_water_peak:

        # If there is a strong CEST signal, it actually results in a bad fit. So 
        # I'm going to find the lowest frequency 
        lowest_freq = numpy.min(numpy.abs(freq_list))
        fitted_freqs = [f for f in freq_list if numpy.abs(f)<lowest_freq+1]
        lowest_freq_ind = sorted([freq_list.index(c) for c in fitted_freqs])
        shift = fit_water_peak(tmp[lowest_freq_ind],fitted_freqs)

        # Interpolation happens here

        s_shifted_back = scipy.interp(new+shift, new, tmp[ind])

        #new_shifted = [n - shift for n in new]

        return new, s_shifted_back

    else:
        # Return the x-axis (new) and the y-axis (tmp)
        return new,tmp[ind]


def h_baselineDiffs(xd, 
                    yd,
                    polynomialDegree=5,
                    removePeaksDict = None):
    
    # Define a function to evaluate the fit given an x-array and coefficients 
    def evaluate_fit(x, coeffs):
        yy = 0
        for i in range(len(coeffs)):
            yy += coeffs[i]*x**(len(coeffs)-i - 1)
        return yy    

    # Initialize xindex (xi), force to array

    xd = numpy.array(xd)    
    yd = numpy.array(yd)
    xi = numpy.ones(shape=[len(xd)])
    
    # If multiple peakes need to be removed
    if removePeaksDict:
    
        for k,v in list(removePeaksDict.items()):
            xi = xi * numpy.where(numpy.all([xd >= v[0],xd <= v[1]],
                                            axis=0),numpy.nan,1)
    
    # Fit a polynomial to xd,yd

    try:
        coeffs = numpy.polyfit(xd[numpy.isfinite(xd*xi)],
                            yd[numpy.isfinite(yd*xi)],polynomialDegree)

    except TypeError:
        coeffs = numpy.ones(shape=[polynomialDegree])*0

    # Compute the difference between the fit and the original data
    functionEval = evaluate_fit(xd, coeffs)
    diffs = functionEval - yd    
    
    # return the x and y coordinates, diffs, and the function eval
    
    return(list(xd), list(yd), list(diffs), list(functionEval))


def fitRemoveBaseline(scn_to_analyse,
                      removePeaksDict=None,
                      polynomialOrder=5,
                      target_ppm_key = '2.0',
                      consider_min=1.,
                      consider_max=4.,
                      ppm_norm=200.,
                      pdata_num = 0,
                      bbox = None,
                      water_shift_scn_name = None,
                      water_shift_adata = 'water_shift_map'):

    scn = sarpy.Scan(scn_to_analyse)

    # Size info
    x_size = scn.pdata[pdata_num].data.shape[0]
    y_size = scn.pdata[pdata_num].data.shape[1]  

    maxVal = numpy.empty(shape=[x_size,y_size])+numpy.nan
    ppmVal = numpy.empty(shape=[x_size,y_size])+numpy.nan

    if removePeaksDict is None:
        removePeaksDict = {'2.5':[1.8,2.4],
                           '3.4':[3.,3.7],
                           'max':[5,0]}
    #### Deal with bounding boxes

    if bbox is None:        
        bbox = numpy.array([0,x_size-1,0,y_size-1])    
    else:
        bbox = scn.adata['bbox'].data

    #### Shift the water peak to 0
    
    if water_shift_scn_name:

        water_shift_scn = sarpy.Scan(water_shift_scn_name)
        water_shift_map = water_shift_scn.adata[water_shift_adata].data
    else:

        water_shift_map = numpy.zeros(shape=[x_size,y_size])

    #### Loop through the bbox'd array
    for xcoord in range(bbox[0],bbox[1]):
        
        for ycoord in range(bbox[2],bbox[3]):

            # shift_water_peak is set to False because there isn't enough data around 0 
            # to successfully shift water peak. Therefore, another method must be used
            
            x,y = cest_spectrum(scn.shortdirname,
                                     xcoord,ycoord,
                                     normalize_to_ppm=ppm_norm,
                                     shift_water_peak=False)

            ## Shift the x values such that the water peak is at 0
            x = x + water_shift_map[xcoord,ycoord]
            x = numpy.array(x)

            xvals,yvals,diffs,fevals = h_baselineDiffs(x,y,polynomialOrder,removePeaksDict)
            xvals = numpy.array(xvals)
            diffs = numpy.array(diffs)
            # This fills the Peak size and location arrays.
            # If the value doesn't work, catch the error
            
            try:
                target_ppm_ind = numpy.where(numpy.all([xvals>=removePeaksDict[target_ppm_key][0],
                                                        xvals<=removePeaksDict[target_ppm_key][1],
                                                        xvals>=consider_min,
                                                        xvals<=consider_max],axis=0))

                maxVal[xcoord,ycoord] = numpy.max(diffs[target_ppm_ind])
                ppmVal[xcoord,ycoord] = xvals[list(diffs).index(numpy.max(diffs[target_ppm_ind]))]
            except ValueError:
                maxVal[xcoord,ycoord] = numpy.nan
                ppmVal[xcoord,ycoord] = numpy.nan               


    return maxVal, ppmVal


def generate_offset_list(additionalDict = None,
                         manuallyInsertedOffsets = None,
                         manuallyInsertedPositions = None,
                         alternateFreqs = True):

    if additionalDict is None:
        additionalDict = collections.OrderedDict([
                               ('dummy',[-60,60,10]),
                               ('start',[-5.1,5.1,0.1]),
                               ('baseline',[-60,60,10]),       
                               ('2.5',[2.2,2.5,0.01]),
                               ('3.4',[3.3,3.5,0.01])
                              ])

    offsetList = []

    for k,v in list(additionalDict.items()):

        offsetList.extend(numpy.round(numpy.arange(v[0],
                                                   v[1],
                                                   v[2]),3))

    # Reorder the list so it's alternating
    # (largest -> smallest -> next largest -> next smallest ->)
    # Gotten from: http://stackoverflow.com/questions/17436870/python-alternating-elements-of-a-sorted-array
    # Of course :-)

    if alternateFreqs is True:
        offsetList = list(sum(list(zip(reversed(offsetList), offsetList)), ())[:len(offsetList)])

    # Now manually insert offsets and frequency
    if manuallyInsertedOffsets is None:
        print([numpy.float("{:.3f}".format(off)) for off in offsetList])
        return numpy.round(offsetList,3)
    else:

        assert len(manuallyInsertedPositions) == len(manuallyInsertedOffsets), "List lengths not the same, please check input lists"
        
        for off,pos in zip(manuallyInsertedOffsets,manuallyInsertedPositions):

            offsetList.insert(pos,off)

        print([numpy.float("{:.3f}".format(off)) for off in offsetList])
        return numpy.round(offsetList,3)            

def cest_vtc(scn_to_analyse):

    import pylab

    scn = sarpy.Scan(scn_to_analyse)
    fig = pylab.figure()

    fig.set_size_inches(20, 20)
    G = pylab.matplotlib.gridspec.GridSpec(1,1, wspace=0.0, hspace=0.0)   
    dat = scn.pdata[0].data

    x_size = dat.shape[0]
    y_size = dat.shape[1]

    # Deal with bounding boxes
    try:
        bbox = scn.adata['bbox'].data
    except KeyError:       
        bbox = numpy.array([0,x_size-1,0,y_size-1])   

    imgdata = numpy.mean(dat,axis=2)
    vtcdata = scn.adata['vtc'].data

    reps = vtcdata.shape[-1] / imgdata.shape[0]

    axs = fig.add_subplot(G[0, 0])
    aspect= (1.0*scn.method.PVM_FovCm[0]/scn.method.PVM_Matrix[0])/ \
            (1.0*scn.method.PVM_FovCm[1]/scn.method.PVM_Matrix[1])

    axs.imshow(imgdata[bbox[0]:bbox[1],\
                       bbox[2]:bbox[3]],\
                       cmap='gray', 
                       interpolation='None',
                       alpha=1.0,
                       zorder=0,
                       aspect=aspect,
                       vmin=0,
                       vmax=1)
    axs.set_axis_off()
    fig.canvas.draw()
    pylab.axis('off')

    box = axs.get_position().bounds
    height = box[3] / (bbox[1]-bbox[0])

    for ht,i in enumerate(range(bbox[0], bbox[1])):

        #, simply use the add_axes() method which takes a list of 
        # [left, bottom, width, height] values in 0-1 relative figure coordinates

        #The add_axes method takes a list of four values, which are 
        #xmin, ymin, dx, and dy for the subplot, where xmin and ymin 
        #are the coordinates of the lower left corner of the subplot, 
        #and dx and dy are the width and height of the subplot, 
        #with all values specified in relative units 
        #(where 0 is left/bottom and 1 is top/right)

        tmpax = fig.add_axes([box[0], box[1]+ht*height,
                             box[2], height])

        tmpax.plot(vtcdata[i,((bbox[2])*reps):((bbox[3])*reps)],
                           color='r', 
                           linewidth=1.5,
                           zorder=1)
        tmpax.set_axis_off()
        pylab.ylim([0.2,1.1])
        pylab.xlim([0,((bbox[3])*reps)-(bbox[2])*reps])

    pylab.savefig('{0}.png'.format(scan_object.shortdirname.split('/')[0]),dpi=600)

    #pylab.close(fig)
    fig = pylab.figure(figsize=(12,8))

    G = pylab.matplotlib.gridspec.GridSpec(1,1, wspace=0.0, hspace=0.0)
    #axs=fig.add_subplot(G[0, 0])
    axs.imshow(numpy.flipud(imgdata[bbox[0]:bbox[1],\
                           bbox[2]:bbox[3]]),\
                           cmap='gray', 
                           interpolation='None',
                           alpha=.5,
                           zorder=0,
                           aspect=aspect)  
    axs.set_axis_off()        

def fit_5_peaks_cest(scn_to_analyse, fitrounds = 1):

    def get_neighbours_starting(fit_arr,i,j):
        #ToFIX: Make this a random game so that any of the pixels are used rather than always the same one.
        # that way if there's a problem, we can make sure it's fixed in the second iteration

        randint = random.randint(1,4)

        if numpy.isfinite(numpy.sum(fit_arr[i,j])): #first try the same px in case this is a refit
            return fit_arr[i,j]
        elif numpy.isfinite(numpy.sum(fit_arr[i,j-1])) and randint ==1:
            return fit_arr[i,j-1]
        elif numpy.isfinite(numpy.sum(fit_arr[i-1,j])) and randint ==2:
            return fit_arr[i-1,j]  
        elif numpy.isfinite(numpy.sum(fit_arr[i-1,j-1])) and randint ==3: 
            return fit_arr[i-1,j-1]
        else:
            peak1dict = {'A':-0.09,
                         'lw':1.3,
                         'w0':2.2}
            peak2dict = {'A':-0.08,
                         'lw':1.0,
                         'w0':3.5}
            peak3dict = {'A':-0.13,
                         'lw':3.5,
                         'w0': -3.3}
            peak4dict = {'A':-0.06,
                         'lw':1.2,
                         'w0': -3.}    
            peak5dict = {'A':-0.9,
                         'lw':1.3,
                         'w0':0.01}

            shift = 1
            paramDict = collections.OrderedDict()
            paramDict['peak1'] = peak1dict
            paramDict['peak2'] = peak2dict
            paramDict['peak3'] = peak3dict
            paramDict['peak4'] = peak4dict
            paramDict['peak5'] = peak5dict
            paramDict['shift'] = shift

            return makeParamArray(paramDict,fixw0=False)  

    def zspectrum_N(params,freqs):

        arr = numpy.empty_like(freqs)*0
        shift =  params[0]

        for i in numpy.arange(0,len(params[1:]),3):

            A = params[i+1]
            lw =  params[i+2]
            w0 = params[i+3]

            tmp = numpy.divide(A,(1+4*((freqs-w0)/lw)**2))
            
            for ind,v in enumerate(tmp):
                if numpy.isnan(v):
                    if numpy.isnan(tmp[ind-1]):
                        tmp[ind] = tmp[ind-2]
                        if numpy.isnan(tmp[ind-2]):

                            print('fail')
                    else:
                        tmp[ind] = tmp[ind-1]
                    
            arr = arr+tmp

        #nanlocs = numpy.isnan(arr)
        #arr[nanlocs] = arr[nanlocs-1]

        if numpy.isnan(numpy.sum(arr)):
            print('hallelujah')#, params
        return (arr+shift)

    def penaltyfn(x,centre=0., scale=1., trough_width=1., steepness=2.):

        return scale*((centre-x)/trough_width)**(2*steepness)

    def h_residual_Zspectrum_N(params, y_data, w):
        penalty = 0

        for i in numpy.arange(1,len(params[1:]),3):
        #return scale*((centre-x)/trough_width)**(2*steepness)

            if i==1: #pk1  
                # lw penalty
                penalty += penaltyfn(params[i+1], centre=0, scale=1E-3, trough_width=.5) # trough 0.3/1/1 ??  
                # w0 penalty
                penalty += penaltyfn(params[i+2], centre=2.2, scale=1E-3, trough_width=.2)
                #print('\t {0}: {1}'.format(i,penalty))
            elif i==4: #pk2
                # lw penalty 
                penalty += penaltyfn(params[i+1], centre=1.0, scale=1E-3, trough_width=.3)
                # w0 penalty
                penalty += penaltyfn(params[i+2], centre=3.6, scale=1E-3, trough_width=.1)    
                #print('\t {0}: {1}'.format(i,penalty))            
            elif i==7: #pk3
                
                # Amplitude < 0 penalty 
                #if params[i] > 0:
                #    penalty += 1E-3          
                # lw penalty 
                penalty += penaltyfn(params[i+1], centre=3.5, scale=1E-2,trough_width=.5)
                # w0 penalty
                penalty += penaltyfn(params[i+2], centre=-3.3, scale=1E-3, trough_width=.3)         
                #print('\t {0}: {1}'.format(i,penalty))

            elif i==10: #pk4
                params[i+0] = 0
                params[i+1] = 0
                params[i+2] = 0

                # Amplitude < 0 penalty 
                #if params[i] < 0:
                #    penalty += 1.
                # lw penalty
                #penalty += penaltyfn(params[i+1], centre=1.2, scale=1E-3, trough_width=.1)
                # w0 penalty
                #penalty += penaltyfn(params[i+2], centre=-3, scale=1E-3, trough_width=.3)
                #print('\t {0}: {1}'.format(i,penalty))

            elif i==13: #pk 5
                # lw penalty
                penalty += penaltyfn(params[i+1], centre=1.5, scale=1e-3, trough_width=0.3)
                # w0 penalty
                penalty += penaltyfn(params[i+2], centre=0.01, scale=1e-3, trough_width=.4)
                #print('\t {0}: {1}'.format(i,penalty))
            #print scipy.nansum(numpy.abs(y_data - zspectrum_N(params,w))), penalty
        return numpy.abs(y_data - zspectrum_N(params,w)) +penalty 

    def makeParamArray(paramDict, fixw0 = False): 
        params = []   
        params.append(paramDict['shift'])

        if fixw0:
            for paramKey,p, in list(paramDict.items()):
                if paramKey != 'shift':
                    params.append(p['A'])
                    params.append(p['lw'])
                else:
                    continue
            return params

        else:
            for paramKey,p, in list(paramDict.items()):
                if paramKey != 'shift':
                    params.append(p['A'])
                    params.append(p['lw'])    
                    params.append(p['w0'])
                else:
                    continue                   
            return params    

    scn = sarpy.Scan(scn_to_analyse)
    pdata_num = 0
    x_size = scn.pdata[pdata_num].data.shape[0]
    y_size = scn.pdata[pdata_num].data.shape[1]  
    try:
        roi = scn.adata['roi'].data
    except KeyError:
        roi = scn.pdata[0].data[:,:,0]*0+1

    # Get the bbox so that the whole image isn't fit
    try:
        bbox = scn.adata['bbox'].data
    except KeyError:       
        bbox = numpy.array([0,x_size-1,0,y_size-1])   

    datashape = roi.shape
    roi_reshaped = numpy.reshape(roi,[roi.shape[0],roi.shape[1],1])
    cestscan_roi = scn.pdata[0].data * roi_reshaped       

    # Fit multiple peaks, need some empty arrays
    offst = numpy.empty_like(roi) + numpy.nan
    pk1_amp = numpy.empty_like(roi) + numpy.nan
    pk1_pos = numpy.empty_like(roi) + numpy.nan
    pk1_width = numpy.empty_like(roi) + numpy.nan

    pk2_amp = numpy.empty_like(roi) + numpy.nan
    pk2_pos = numpy.empty_like(roi) + numpy.nan
    pk2_width = numpy.empty_like(roi) + numpy.nan

    pk3_amp = numpy.empty_like(roi) + numpy.nan
    pk3_pos = numpy.empty_like(roi) + numpy.nan
    pk3_width = numpy.empty_like(roi) + numpy.nan

    pk4_amp = numpy.empty_like(roi) + numpy.nan
    pk4_pos = numpy.empty_like(roi) + numpy.nan
    pk4_width = numpy.empty_like(roi) + numpy.nan

    pk5_amp = numpy.empty_like(roi) + numpy.nan
    pk5_pos = numpy.empty_like(roi) + numpy.nan
    pk5_width = numpy.empty_like(roi) + numpy.nan

    fit_quality = numpy.empty_like(roi) + numpy.nan
    fit_params_arr = numpy.empty_like(roi, dtype=object)
    fit_params_arr[:] = [numpy.nan]
    ppm_corrected_arr = numpy.empty_like(roi, dtype=object)

    # Defining parameters

    freq_list = scn.method.CEST_FreqListPPM

    ppm_limit_min = -50
    ppm_limit_max = 50
    exclude_ppm = 200
    normalize_to_ppm = 66.6
    

    possibleNormalizations = [i for i, x in enumerate(freq_list) if numpy.abs(x - normalize_to_ppm) <1E-4]
    normalizeTo = possibleNormalizations[-1]

    # Get only the frequencies within the ppm_limit
    ppm_filtered = [f for f in freq_list if ppm_limit_max > f > ppm_limit_min]

    # Exclude the dummy frequencies at the beginning (66.6 ppm)
    ppm_filtered = sorted([n for n in ppm_filtered if n!= exclude_ppm])

    # Get the index of the good frequencies relative to the original list 
    ppm_filtered_ind = [freq_list.index(c) for c in ppm_filtered]  

    # get the freqs that'll be used for water fit
    water_fit_freqs = [f for f in ppm_filtered if (numpy.abs(f)< 3.)]
    water_fit_freqs_ind = sorted([ppm_filtered.index(c) for c in water_fit_freqs])

    # Create some empty arrays
    water_shifts = numpy.empty_like(roi) + numpy.nan
    new_shifted = numpy.empty(shape=(water_shifts.shape[0], water_shifts.shape[0], len(ppm_filtered))) + numpy.nan

    # Fit count, this counts the number of rounds the data has been fit
    fitcount = 0

    while fitcount < fitrounds:
        for xval in range(bbox[0],bbox[1]):    
            for yval in range(bbox[2],bbox[3]):
                # Get the data and normalize it to index of normalize_to_ppm
                tmp = cestscan_roi[xval,yval][ppm_filtered_ind] / scn.pdata[0].data[xval,yval,normalizeTo]           

                # Check to make sure I'm inside the ROI
                if numpy.isfinite(numpy.sum(tmp)):            
                    # First do the water fit and shift the data so water is at 0  
                    shiftParams = fit_water_peak(tmp[water_fit_freqs_ind],water_fit_freqs,allParams=True)
                    shift = shiftParams[3]
                    water_shifts[xval,yval] = shift

                    # Interpolation happens here
                    if numpy.isfinite(shift):
                        s_shifted_back = scipy.interp(ppm_filtered+shift, ppm_filtered, tmp)
                        new_shifted[xval,yval,:] = s_shifted_back       
                    else:
                        print(shift)
                        pass            

                    testParams = get_neighbours_starting(fit_params_arr,xval,yval)

                    fit_params,cov,infodict,mesg,ier = scipy.optimize.leastsq(
                                                                h_residual_Zspectrum_N,
                                                                testParams,
                                                                args=(new_shifted[xval,yval], ppm_filtered), 
                                                                full_output = True,
                                                                maxfev = 900,
                                                                ftol =1E-9)
                    w = numpy.arange(-20.,20.,0.01)

                    # Specify paramsets for peaks:
                    #TOFIX: why is the offset applied to each peak
                    pk1 = [fit_params[0]]+list(fit_params[1:4])
                    pk2 = [fit_params[0]]+list(fit_params[4:7])
                    pk3 = [fit_params[0]]+list(fit_params[7:10])
                    pk4 = [fit_params[0]]+list(fit_params[10:13])
                    pk5 = [fit_params[0]]+list(fit_params[13:16]) 

                    offst[xval,yval] = fit_params[0]
                    pk1_amp[xval,yval] = fit_params[1]
                    pk1_width[xval,yval] = fit_params[2]
                    pk1_pos[xval,yval] = fit_params[3]

                    pk2_amp[xval,yval] = fit_params[4]
                    pk2_width[xval,yval] = fit_params[5]
                    pk2_pos[xval,yval] = fit_params[6]

                    pk3_amp[xval,yval] = fit_params[7]
                    pk3_width[xval,yval] = fit_params[8]
                    pk3_pos[xval,yval] = fit_params[9]

                    pk4_amp[xval,yval] = fit_params[10]
                    pk4_width[xval,yval] = fit_params[11]            
                    pk4_pos[xval,yval] = fit_params[12]

                    pk5_amp[xval,yval] = fit_params[13]
                    pk5_width[xval,yval] = fit_params[14]            
                    pk5_pos[xval,yval] = fit_params[15]                
                  
                    fit_quality[xval,yval] = scipy.nansum(numpy.abs(new_shifted - zspectrum_N(fit_params,ppm_filtered)))
                    fit_params_arr[xval,yval] = fit_params
                    ppm_corrected_arr[xval,yval] = ppm_filtered

        fitcount+=1 # increment fitcounter
    
    # Save the data as a structured array
    newstruct = numpy.empty(roi.shape, dtype=[('offset', 'float32'),
       ('A1', 'float32'),('w1', 'float32'),('p1', 'float32'),
       ('A2', 'float32'),('w2', 'float32'),('p2', 'float32'),
       ('A3', 'float32'),('w3', 'float32'),('p3', 'float32'),
       ('A4', 'float32'),('w4', 'float32'),('p4', 'float32'),
       ('A5', 'float32'),('w5', 'float32'),('p5', 'float32')])

    # Nan the array so there are no zeroes anywhere
    newstruct[:] = numpy.nan

    newstruct['offset'] = offst
    newstruct['A1'] = pk1_amp
    newstruct['w1'] = pk1_width
    newstruct['p1'] = pk1_pos
    newstruct['A2'] = pk2_amp
    newstruct['w2'] = pk2_width
    newstruct['p2'] = pk2_pos
    newstruct['A3'] = pk3_amp
    newstruct['w3'] = pk3_width
    newstruct['p3'] = pk3_pos
    newstruct['A4'] = pk4_amp
    newstruct['w4'] = pk4_width
    newstruct['p4'] = pk4_pos
    newstruct['A5'] = pk5_amp
    newstruct['w5'] = pk5_width
    newstruct['p5'] = pk5_pos

    return {'':newstruct,'fit_quality':fit_quality}


