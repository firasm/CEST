# CEST Analysis Code

####### Imports #######


from __future__ import division
import numpy
import sarpy
import sarpy.fmoosvi.analysis
import scipy


####### Fitting Functions #######

def h_zspectrum_N(params,freqs):
    
    arr = numpy.empty_like(freqs)*0
    shift =  params[0]
    
    for i in numpy.arange(0,len(params[1:]),3):
        
        A = params[i+1] # Amplitude of the peak
        w0 = params[i+2] # Resonant frequency of sample
        lw =  params[i+3] # Line width (in units of freqs)

        # Construct the multiple lorentzians added together
        tmp = numpy.divide(A,(1+4*((freqs-w0)/lw)**2))
        
        arr = arr+tmp
    return (arr+shift)

def h_residual_Zspectrum_N(params, y_data, w):
    
    return numpy.abs(y_data - h_zspectrum_N(params,w))    

def fit_water_peak(data,offset_freqs):

    """
    Fits a Lorentzian to the data, and 
    returns the water offset frequency

    """

    # Presumably the maximum signal will be at 
    # some large frequency offset, so let's just
    # use that as our baseline parameter

    # Also, normalize the data so the fit is easier

    params_passed = [numpy.max(data),-1.,0.1,0.6]

    fit_params,cov,infodict,mesg,ier = scipy.optimize.leastsq(h_residual_Zspectrum_N,
        params_passed,
        args=(data, offset_freqs),
        full_output = True,
        maxfev = 200)

    return fit_params[2]

def shift_water_peak(scn_to_analyse=None, 
                     bbox = None,
                     pdata_num = 0,
                     **kwargs):
    
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
        bbox = sarpy.fmoosvi.getters.convert_bbox(scn_to_analyse,bbox) 

    # Shape the bbox properly to account for offsets
    if bbox.shape == (4,):            
    
        bbox_mask = numpy.empty([x_size,y_size])
        bbox_mask[:] = numpy.nan        
        bbox_mask[bbox[0]:bbox[1],bbox[2]:bbox[3]] = 1
    
    # Tile for Offsets
    bbox_mask = numpy.tile(bbox_mask.reshape(x_size,y_size,1),offsets) 

    # Fit the water peak pixel by pixel
    water_shift_map = numpy.nan + numpy.empty([x_size,y_size])

    dat = scn.pdata[0].data
    for x in xrange(bbox[0],bbox[1]):
        for y in xrange(bbox[2],bbox[3]):

            fit_data = dat[x,y] / numpy.max(dat[x,y])

            water_shift_map[x,y] = fit_water_peak(fit_data,freq_list)

    return water_shift_map

####### Displaying and Plotting Functions #######    

def cest_spectrum(scn_to_analyse,
                  xval,
                  yval,
                  shift_water_peak = True,
                  normalize_to_ppm = 200,
                  ppm_limit = 50,
                  exclude_ppm = 66.6,
                  pdata_num = 0):
    
    """
    Does the grunt work of processing freqs to give
    a sequential list of frequencies (rather than alternating).

    :param str scn_to_analyse: scan.shortdirname
    :param bool shift_water_peak: Fits the water peak and shifts the freqs so that 
            water is at 0 ppm (you can now average after this)
    :param float normalize_to_ppm: freq to normalize to
    :param float ppm_limit: Only return frequencies within +ppm_limit and -ppm_limit
    :param float exclude_ppm: Exclude the dummy frequencies
    :param int pdata_num: reconstruction number, according to python numbering.
            default reconstruction is pdata_num = 0.

    :return: new (frequencies in sequential order)
             tmp (data at frequencies in 'new')
    """ 
    scn = sarpy.Scan(scn_to_analyse)

    # Get the Frequencies and make them sequential 
    # so that alternating points don't plot badly

    # First get the freq list
    freq_list = scn.method.CEST_FreqListPPM

    # Find the frequency to normalize to, throw error if not found
    normalizeTo = freq_list.index(normalize_to_ppm)

    # Get only the frequencies within the ppm_limit
    new = [f for f in freq_list if numpy.abs(f)<ppm_limit]

    # Exclude the dummy frequencies at the beginning (66.6 ppm)
    new = sorted([n for n in new if n!= exclude_ppm])

    # Get the index of the good frequencies relative to the original list 
    ind = [freq_list.index(c) for c in new]  

    # Get the data and normalize it to index of normalize_to_ppm
    tmp = scn.pdata[pdata_num].data[xval,yval,ind] / scn.pdata[0].data[xval,yval,normalizeTo] 

    if shift_water_peak:
        shift = fit_water_peak(tmp,new)
        new_shifted = [n - shift if n>=0 else n+shift for n in new]

        return new_shifted, tmp

    else:
        # Return the x-axis (new) and the y-axis (tmp)
        return new,tmp








