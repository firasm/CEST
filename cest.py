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
                  normalize=True,
                  normalize_to_ppm = 200,
                  ppm_limit_min = -50,
                  ppm_limit_max = 50,
                  exclude_ppm = 66.6,
                  pdata_num = 0):
    
    """
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
    """ 
    scn = sarpy.Scan(scn_to_analyse)

    # Get the Frequencies and make them sequential 
    # so that alternating points don't plot badly

    # First get the freq list
    freq_list = scn.method.CEST_FreqListPPM

    if normalize:
    # Find the frequency to normalize to, throw error if not found
        possibleNormalizations = [i for i, x in enumerate(freq_list) if x == normalize_to_ppm]

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

def plotBaselineDiffs(xrawdata, 
                      yrawdata,
                      polynomialDegree=5,
                      removePeaksDict = None):
    
    # Define a function to evaluate the fit given an x-array and coefficients 
    def evaluate_fit(x, coeffs):
        yy = 0
        for i in xrange(len(coeffs)):
            yy += coeffs[i]*x**(len(coeffs)-i - 1)
        return yy    

    # Initialize xd (xdata) and xindex (xi), force to array
    xrawdata = numpy.array(xrawdata)
    xd = xrawdata
    xi = numpy.ones(shape=[len(xd)])
    
    # If multiple peakes need to be removed
    if removePeaksDict:
    
        for k,v in removePeaksDict.iteritems():
            xi = xi * numpy.where(numpy.all([xd >= v[0],xd <= v[1]],
                                            axis=0),numpy.nan,1)

    # Get the corresponding values of y,x after it has been cleaned, force to array
    yd = numpy.array(yrawdata)
    xd = numpy.array(xrawdata)
    
    # Fit a polynomial of degree 5 to xd,yd

    try:
        coeffs = numpy.polyfit(xd[numpy.isfinite(xd*xi)],
                            yd[numpy.isfinite(yd*xi)],polynomialDegree)

    except TypeError:
        coeffs = numpy.ones(shape=[polynomialDegree])*0

    # Compute the difference between the fit and the original data
    diffs = evaluate_fit(xrawdata, coeffs) - yd

    functionEval = evaluate_fit(xrawdata, coeffs)
    
    # return the x and y coordinates, and the function eval
    
    return(list(xd), list(yd), list(diffs), list(functionEval))






