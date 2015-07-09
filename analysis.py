# CEST Analysis Code

####### Imports #######


from __future__ import division
import numpy
import sarpy
import sarpy.fmoosvi.analysis
import scipy
import collections

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
        bbox = scn.adata['bbox'].data

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


def h_baselineDiffs(xd, 
                    yd,
                    polynomialDegree=5,
                    removePeaksDict = None):
    
    # Define a function to evaluate the fit given an x-array and coefficients 
    def evaluate_fit(x, coeffs):
        yy = 0
        for i in xrange(len(coeffs)):
            yy += coeffs[i]*x**(len(coeffs)-i - 1)
        return yy    

    # Initialize xindex (xi), force to array

    xd = numpy.array(xd)    
    yd = numpy.array(yd)
    xi = numpy.ones(shape=[len(xd)])
    
    # If multiple peakes need to be removed
    if removePeaksDict:
    
        for k,v in removePeaksDict.iteritems():
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
    for xcoord in xrange(bbox[0],bbox[1]):
        
        for ycoord in xrange(bbox[2],bbox[3]):

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

    for k,v in additionalDict.iteritems():

        offsetList.extend(numpy.round(numpy.arange(v[0],
                                                   v[1],
                                                   v[2]),3))

    # Reorder the list so it's alternating
    # (largest -> smallest -> next largest -> next smallest ->)
    # Gotten from: http://stackoverflow.com/questions/17436870/python-alternating-elements-of-a-sorted-array
    # Of course :-)

    if alternateFreqs is True:
        offsetList = sum(zip(reversed(offsetList), offsetList), ())[:len(offsetList)]

    # Now manually insert offsets and frequency
    if manuallyInsertedOffsets is None:
        print [numpy.float("{:.3f}".format(off)) for off in offsetList]
        return numpy.round(offsetList,3)
    else:

        assert len(manuallyInsertedPositions) == len(manuallyInsertedOffsets), "List lengths not the same, please check input lists"
        
        for off,pos in zip(manuallyInsertedOffsets,manuallyInsertedPositions):

            offsetList.insert(pos,off)

        print [numpy.float("{:.3f}".format(off)) for off in offsetList]
        return numpy.round(offsetList,3)            

def cest_vtc(scn_to_analyse):

    import pylab

    scan_object = sarpy.Scan(scn_to_analyse)
    fig = pylab.figure()

    fig.set_size_inches(20, 20)
    G = pylab.matplotlib.gridspec.GridSpec(1,1, wspace=0.0, hspace=0.0)   
    dat = scan_object.pdata[0].data

    x_size = dat.shape[0]
    y_size = dat.shape[1]

    # Deal with bounding boxes
    try:
        bbox = scan_object.adata['bbox'].data
    except KeyError:       
        bbox = numpy.array([0,x_size-1,0,y_size-1])   

    imgdata = numpy.mean(dat,axis=2)
    vtcdata = scan_object.adata['vtc'].data

    reps = vtcdata.shape[-1] / imgdata.shape[0]

    axs = fig.add_subplot(G[0, 0])
    aspect= (1.0*scan_object.method.PVM_FovCm[0]/scan_object.method.PVM_Matrix[0])/ \
            (1.0*scan_object.method.PVM_FovCm[1]/scan_object.method.PVM_Matrix[1])

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

    for ht,i in enumerate(xrange(bbox[0], bbox[1])):

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
    fig = pylab.figure()

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

