# CEST Analysis Code

####### Imports #######

import numpy
import sarpy
import sarpy.analysis.analysis
import scipy
import collections
import pylab
import random
import cest

####### Fitting Functions #######

def h_convertBetweenStructArrays(params, toType = 'struct'):

    if toType =='array': # probably starting with a struct

        if len(params.dtype) < 1: # check if already an array
            return params

        paramAsArray = []
        shift =  params['offset']
        paramAsArray.extend(shift)

        # Now get the other peaks as regular lorentzians
        for i in numpy.arange(1,int(len(params.dtype.descr)/3)):

            # Genius of @drSar, did not know you could access structuredArrays this way
            A = params['A{}'.format(i)]
            w = params['w{}'.format(i)]
            p = params['p{}'.format(i)]

            paramAsArray.extend(A)
            paramAsArray.extend(w)
            paramAsArray.extend(p)

        # Water parameters
        A = params['water_A']
        w = params['water_w']
        p = params['water_p']

        paramAsArray.extend(A)
        paramAsArray.extend(w)
        paramAsArray.extend(p)

        return numpy.array(paramAsArray)

    if toType =='struct': # probably starting with a struct

        if len(params.dtype) > 1: # check if already a struct
            return params

        newstruct = numpy.zeros((1), dtype=[('offset', 'float64'),
           ('A1', 'float64'),('w1', 'float64'),('p1', 'float64'),
           ('A2', 'float64'),('w2', 'float64'),('p2', 'float64'),
           ('A3', 'float64'),('w3', 'float64'),('p3', 'float64'),
           ('A4', 'float64'),('w4', 'float64'),('p4', 'float64'),
           ('water_A', 'float64'),('water_w', 'float64'),('water_p', 'float64')])

        if len(params) < 5: # water fit
            newstruct['offset'] = params[0]
            newstruct['water_A'] = params[1]
            newstruct['water_w']= params[2]
            newstruct['water_p'] = params[3]

        else: # all other fits
            newstruct['offset'] = params[0]
            newstruct['A1']= params[1]
            newstruct['w1']= params[2]
            newstruct['p1']= params[3]
            newstruct['A2']= params[4]
            newstruct['w2']= params[5]
            newstruct['p2']= params[6]
            newstruct['A3']= params[7]
            newstruct['w3']= params[8]
            newstruct['p3'] = params[9]
            newstruct['A4']= params[10]
            newstruct['w4']= params[11]
            newstruct['p4']= params[12]
            newstruct['water_A'] = params[13]
            newstruct['water_w']= params[14]
            newstruct['water_p'] = params[15]
        return newstruct

def h_lorentzian(A,w,p,freqs):
    return numpy.divide(-A,(1+4*((freqs-p)/w)**2))

def h_superlorentzian(A,w,p,freqs):
    theta=numpy.arange(0,1.01,0.01) * numpy.pi/2
    s=0;
    for i in range(len(numpy.arange(0,1.01,0.01))):
        trig=numpy.abs(3*numpy.cos(theta[i])**2-1) # Trig Identity with cos^2
        temp=numpy.sin(theta[i])*numpy.sqrt(2/numpy.pi)*1/(w*trig)
        temp*=numpy.exp(-2*((2*numpy.pi*(freqs+p)*(1/w))/trig)**2)
        s+=temp
    return A*s/numpy.max(s)

def h_peak_N(params,freqs,peak = 1):
    ''' Zspectrum N function to only return one peak'''

    paramStructuredArray = h_convertBetweenStructArrays(params,toType = 'struct')

    # First get the water peak as a super lorentzian
    shift =  paramStructuredArray['offset']

    # ^ needs to change to -1 and additional things added for water peak (super lorentzian)
    # Genius of @drSar, did not know you could access structuredArrays this way
    A = paramStructuredArray['A{}'.format(peak)][0]+1E-5
    w = paramStructuredArray['w{}'.format(peak)][0]+1E-6
    p = paramStructuredArray['p{}'.format(peak)][0]+1E-7

    return h_lorentzian(A,w,p,freqs) + shift    

def h_zspectrum_N(params,freqs):
    ''' Updated Zspectrum N function to now require a structured array of params'''

    paramStructuredArray = h_convertBetweenStructArrays(params,toType = 'struct')

    arr = numpy.zeros_like(freqs)

    # First get the water peak as a super lorentzian
    shift =  paramStructuredArray['offset']

    # Now get the other peaks as regular lorentzians
    for i in numpy.arange(1,int(len(paramStructuredArray.dtype.descr)/3)):
        # ^ needs to change to -1 and additional things added for water peak (super lorentzian)
        # Genius of @drSar, did not know you could access structuredArrays this way
        A = paramStructuredArray['A{}'.format(i)][0]+1E-5
        w = paramStructuredArray['w{}'.format(i)][0]+1E-6
        p = paramStructuredArray['p{}'.format(i)][0]+1E-7

        newpeak = h_lorentzian(A,w,p,freqs)

        if numpy.isnan(numpy.sum(newpeak)):
            print(A,w,p)
            print('struct: {}'.format(paramStructuredArray['A{}'.format(i)],paramStructuredArray['w{}'.format(i)],paramStructuredArray['p{}'.format(i)]))
        arr += newpeak

    A = paramStructuredArray['water_A'][0]+1E-8
    w = paramStructuredArray['water_w'][0]+1E-9
    p = paramStructuredArray['water_p'][0]+1E-10

    arr += h_lorentzian(A,w,p,freqs)
    return arr + shift

def h_residual_Zspectrum_N(params, y_data, w):

    # Define the penalty function:
    #def old penaltyfn(x,centre=0., scale=1., trough_width=1., steepness=2.):
    #    return scale*((centre-x)/trough_width)**(2*steepness)

    params = h_convertBetweenStructArrays(params,toType = 'struct')

    def penaltyfn(x, centre, trough = 0.1, scale=1E-4,hardlimit = -50):

        #centre = numpy.mean([rightEdge,leftEdge])
        pen = scale*((centre-x)/trough)**2*2
        #print(pen)
        return 0
        return pen

        #return numpy.where(numpy.logical_or(x < leftEdge,x > rightEdge),
        #                   penaltyScale*(numpy.abs(x-center))**2,
        #                  1E-2*(numpy.abs(x-center))**2)

    penaltyStruct = numpy.zeros((2), dtype=[('w1', 'float64'),('p1', 'float64'),
                                             ('w2', 'float64'),('p2', 'float64'),
                                             ('w3', 'float64'),('p3', 'float64'),
                                             ('w4', 'float64'),('p4', 'float64'),
                                             ('water_w', 'float64'),('water_p', 'float64')])
    # Theoretical centers of peaks
    penaltyStruct['p1'] = [2.2,0.2]
    penaltyStruct['p2'] = [3.6,0.1]
    penaltyStruct['p3'] = [-3.3,0.3]
    penaltyStruct['p4'] = [-3,0.3]
    penaltyStruct['water_p'] = [0.01,0.2]

    # Theoretical widths of peaks
    penaltyStruct['w1'] = [1.3,0.1]
    penaltyStruct['w2'] = [1.0,0.3]
    penaltyStruct['w3'] = [3.5,0.5]
    penaltyStruct['w4'] = [5.,0.4]
    penaltyStruct['water_w'] = [1.5,0.3]

    penalty = 0
    # Penalties for other peaks
    for i in numpy.arange(1,int(len(params.dtype.descr)/3)):

        penalty += penaltyfn(params['w{}'.format(i)],
                             centre = penaltyStruct['w{}'.format(i)][0],
                             trough = penaltyStruct['w{}'.format(i)][1],
                             hardlimit = 0)

        penalty += penaltyfn(params['p{}'.format(i)],
                             centre = penaltyStruct['w{}'.format(i)][0],
                             trough = penaltyStruct['w{}'.format(i)][1])

    # Penalties for water peak fitting 

    penalty += penaltyfn(params['water_w'],
                             centre = 1.5,
                             trough = 0.3,
                             hardlimit = 0) 

    penalty += penaltyfn(params['water_p'],
                         centre = 0.01,
                         trough = 0.4)

    params = h_convertBetweenStructArrays(params,toType = 'array')

    return numpy.abs(y_data - h_zspectrum_N(params,w)) + penalty

def fit_water_peak(data,offset_freqs,allParams = False):

    """
    Fits a Super Lorentzian to the data, and 
    returns the water offset frequency

    """
    params = lmfit.Parameters()

    # add with tuples: (NAME VALUE VARY MIN  MAX  EXPR  BRUTE_STEP)    
    params.add_many(('DCoffset', 1, True, None, None, None),
                    ('A_0', -0.8, True, None, None, None),
                    ('w_0', 1.3, True, None, None, None),
                    ('p_0', 0.01, True, None, None, None),
    )        

    try:
        out = lmfit.minimize(h_zspectrum_N_residuals, params, args=(numpy.array(offset_freqs), data))
        return numpy.array(out.params['p_0'].value)
    except ValueError: # when the data given is nonsense
        return numpy.array(0)

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


def plotCestPeaks(scn_to_analyse, cest_adata):
    ''' This function takes in a full fit and returns a plot of the individual peaks plotted on the flipped axis'''

    scn = sarpy.Scan(scn_to_analyse)
    result = scn.adata[cest_adata].data
    freqs = numpy.arange(-20,20,0.1)

    pylab.figure(figsize=(10,10))
    pylab.subplot(221)
    pylab.plot(ppm_filtered,new_shifted,'.',label='raw')
    pylab.plot(freqs,cest.analysis.h_zspectrum_N(fit_params,freqs),label='fit')
    pylab.xlim(10,-10)

    fps = cest.analysis.h_convertBetweenStructArrays(fit_params,toType='struct')
    for i in numpy.arange(1,int(len(fps.dtype.descr)/3)):
        pylab.plot(freqs,1-cest.analysis.h_peak_N(fit_params,freqs,i))

    pylab.axvline(2.2,ymin=0,label='2.2 amine',color='y', alpha=0.4)
    pylab.axvline(3.5,ymin=0,label='3.5 amide',color='r', alpha=0.4)
    pylab.axvline(-3.25,ymin=0,label='-3.25 aliphatic',color='g', alpha=0.4)
    pylab.axvline(-3.0,ymin=0,label='-3.0 ppm-amine')    
    #pylab.axvline(1.5,ymin=0.6,label='1.5 OH',color='b', alpha=0.4)

    pylab.subplot(222)
    pylab.plot(ppm_filtered,new_shifted,'.',label='raw')
    pylab.plot(freqs,cest.analysis.h_zspectrum_N(fit_params,freqs),label='fit')

    pylab.xlim(1,-1)
    #pylab.ylim(0,1)
    pylab.legend()











    ###############################

    shift =  paramStructuredArray['offset']

    w = numpy.arange(-20.,20.,0.01)

    pk1 = [fit_params[0]]+list(fit_params[1:4])
    pk2 = [fit_params[0]]+list(fit_params[4:7])
    pk3 = [fit_params[0]]+list(fit_params[7:10])
    pk4 = [fit_params[0]]+list(fit_params[10:13])
    pk5 = [fit_params[0]]+list(fit_params[13:16]) 

    # Separate this out into a different function - this is absurd.
    pylab.figure(figsize=(12,8))                         

    pylab.plot(w,zspectrum_N(pk1,w),'-',label='w0 = {0}, lw = {1}, A={2}'.format(numpy.round(pk1[-1],2),numpy.round(pk1[2],2),numpy.round(pk1[1],2)),color='g') # peak1
    pylab.plot(w,zspectrum_N(pk2,w),'-',label='w0 = {0}, lw = {1}, A={2}'.format(numpy.round(pk2[-1],2),numpy.round(pk2[2],2),numpy.round(pk2[1],2)),color='r') # peak2
    pylab.plot(w,zspectrum_N(pk3,w),'-',label='w0 = {0}, lw = {1}, A={2}'.format(numpy.round(pk3[-1],2),numpy.round(pk3[2],2),numpy.round(pk3[1],2)),color='c') # peak3
    pylab.plot(w,zspectrum_N(pk4,w),'-',label='w0 = {0}, lw = {1}, A={2}'.format(numpy.round(pk4[-1],2),numpy.round(pk4[2],2),numpy.round(pk4[1],2)),color='y') # peak4
    pylab.plot(w,zspectrum_N(pk5,w),'-',label='w0 = {0}, lw = {1}, A={2}'.format(numpy.round(pk5[-1],2),numpy.round(pk5[2],2),numpy.round(pk5[1],2)),color='pink') # peak5

    # Draw vertical lines at peak positions
    pylab.axvline(pk1[-1],color='g')
    pylab.axvline(pk2[-1],color='r')
    pylab.axvline(pk3[-1],color='c')
    pylab.axvline(pk4[-1],color='y')
    #axvline(pk5[-1],color='purple')                

    pylab.plot(ppm_filtered,data_watersupp,'o-',color='pink',label='raw')
    pylab.xlim(15,-15)
    #pylab.ylim(-0.5,0.3)
    pylab.title('Fit for pixel {0},{1} \n Residual: {2}'.format(xval,yval,fit_quality[xval,yval]))
    pylab.legend(loc='lower right')

    pylab.savefig('pixelBypixel/{0}_{1}-{2},{3}.png'.format(scn.patientname,scn.studyname,xval,yval,dpi=400))
    pylab.clf()


################################################

####### Fitting CEST data Functions #######    

################################################

def get_neighbours_starting(fit_arr=None,i=1,j=1):
    #ToFIX: Make this a random game so that any of the pixels are used rather than always the same one.
    # that way if there's a problem, we can make sure it's fixed in the second iteration

    if fit_arr is None: # in case there is nothing specified for fit_arr
        fit_arr = numpy.array([[numpy.nan]*4, [numpy.nan]*4])

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
        newstruct = numpy.zeros((1), dtype=[('offset', 'float64'),
           ('A1', 'float64'),('w1', 'float64'),('p1', 'float64'),
           ('A2', 'float64'),('w2', 'float64'),('p2', 'float64'),
           ('A3', 'float64'),('w3', 'float64'),('p3', 'float64'),
           ('A4', 'float64'),('w4', 'float64'),('p4', 'float64'),
           ('water_A', 'float64'),('water_w', 'float64'),('water_p', 'float64')])

        # Nan the array so there are no zeroes anywhere
        newstruct[:] = numpy.float64(numpy.nan)

        newstruct['offset'] = 1
        newstruct['A1'] = 0.09
        newstruct['w1'] = 1.3
        newstruct['p1'] = 2.2
        newstruct['A2'] = 0.08
        newstruct['w2'] = 1.0
        newstruct['p2'] = 3.5
        newstruct['A3'] = 0.13
        newstruct['w3'] = 3.5
        newstruct['p3'] = -3.3
        newstruct['A4'] = 0.06
        newstruct['w4'] = 1.2
        newstruct['p4'] = -3.
        newstruct['water_A'] = 0.9
        newstruct['water_w'] = 1.3
        newstruct['water_p'] = 0.01

        return newstruct   

def fit_px_cest(scn_to_analyse, xval, yval, fitrounds = 1):

    scn = sarpy.Scan(scn_to_analyse)
    pdata_num = 0 

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
    water_shifts = numpy.zeros(shape=(1))  + numpy.nan
    new_shifted = numpy.zeros(shape=(1)) + numpy.nan

    newstruct = numpy.zeros((1), dtype=[('offset', 'float64'),
       ('A1', 'float64'),('w1', 'float64'),('p1', 'float64'),
       ('A2', 'float64'),('w2', 'float64'),('p2', 'float64'),
       ('A3', 'float64'),('w3', 'float64'),('p3', 'float64'),
       ('A4', 'float64'),('w4', 'float64'),('p4', 'float64'),
       ('water_A', 'float64'),('water_w', 'float64'),('water_p', 'float64')])

    # Fit count, this counts the number of rounds the data has been fit
    fitcount = 0

    while fitcount < fitrounds:
        # Get the data and normalize it to index of normalize_to_ppm
        tmp = scn.pdata[0].data[xval,yval,:][ppm_filtered_ind] / scn.pdata[0].data[xval,yval,normalizeTo]           
   
        # First do the water fit and shift the data so water is at 0  
        shiftParams = fit_water_peak(tmp[water_fit_freqs_ind],water_fit_freqs,allParams=True)
        shift = shiftParams[3]

        # Interpolating the Y-data so that it gets shifted to the acquired offsets!
        if numpy.isfinite(shift):
            s_shifted_back = scipy.interp(ppm_filtered, ppm_filtered+shift/2, tmp)
            new_shifted = s_shifted_back       
        else:
            print(shift)
            pass            

        if fitcount>0: # Use parameters from last fit 
            testParams = h_convertBetweenStructArrays(newstruct,toType='array')

        else: # Get initial starting parameters
            testParams = get_neighbours_starting()
            testParams = h_convertBetweenStructArrays(testParams,toType = 'array')

        fit_params,cov,infodict,mesg,ier = scipy.optimize.leastsq(
                                                    h_residual_Zspectrum_N,
                                                    testParams,
                                                    args=(new_shifted, ppm_filtered), 
                                                    full_output = True,
                                                    maxfev = 900,
                                                    ftol =1E-9)
        newstruct['offset'] = fit_params[0]
        newstruct['A1'] = fit_params[1]
        newstruct['w1'] = fit_params[2]
        newstruct['p1'] = fit_params[3]
        newstruct['A2'] = fit_params[4]
        newstruct['w2'] = fit_params[5]
        newstruct['p2'] = fit_params[6]
        newstruct['A3'] = fit_params[7]
        newstruct['w3'] = fit_params[8]
        newstruct['p3'] = fit_params[9]
        newstruct['A4'] = fit_params[10]
        newstruct['w4'] = fit_params[11]
        newstruct['p4'] = fit_params[12]
        newstruct['water_A'] = fit_params[13]
        newstruct['water_w'] = fit_params[14]
        newstruct['water_p'] = fit_params[15]
      
        fitcount+=1
    freqs = numpy.arange(-10,10,0.1)
    fit_quality = scipy.nansum(numpy.abs(new_shifted - h_zspectrum_N(fit_params,ppm_filtered)))

    return {'fit_params':newstruct, 
            'data': [ppm_filtered, new_shifted],
            'fitdata': [freqs, h_zspectrum_N(newstruct, freqs)],
            'fit_quality': fit_quality}

def fit_5_peaks_cest(scn_to_analyse, fitrounds = 1):

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

    water_amp = numpy.empty_like(roi) + numpy.nan
    water_pos = numpy.empty_like(roi) + numpy.nan
    water_width = numpy.empty_like(roi) + numpy.nan

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

    # Create temporary structs to store paramArrays
    tempstruct = numpy.zeros((1), dtype=[('offset', 'float64'),
       ('A1', 'float64'),('w1', 'float64'),('p1', 'float64'),
       ('A2', 'float64'),('w2', 'float64'),('p2', 'float64'),
       ('A3', 'float64'),('w3', 'float64'),('p3', 'float64'),
       ('A4', 'float64'),('w4', 'float64'),('p4', 'float64'),
       ('water_A', 'float64'),('water_w', 'float64'),('water_p', 'float64')])

    newstruct = numpy.zeros(roi.shape, dtype=[('offset', 'float64'),
       ('A1', 'float64'),('w1', 'float64'),('p1', 'float64'),
       ('A2', 'float64'),('w2', 'float64'),('p2', 'float64'),
       ('A3', 'float64'),('w3', 'float64'),('p3', 'float64'),
       ('A4', 'float64'),('w4', 'float64'),('p4', 'float64'),
       ('water_A', 'float64'),('water_w', 'float64'),('water_p', 'float64')])

    # Nan the array so there are no zeroes anywhere
    newstruct[:] = numpy.nan
    #tempstruct[:] = numpy.nan

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
                        s_shifted_back = scipy.interp(ppm_filtered, ppm_filtered+shift/2, tmp)
                        new_shifted[xval,yval,:] = s_shifted_back       
                    else:
                        print(shift)
                        pass            

                    testParams = get_neighbours_starting(fit_params_arr,xval,yval)
                    testParams = h_convertBetweenStructArrays(testParams,toType = 'array')

                    fit_params,cov,infodict,mesg,ier = scipy.optimize.leastsq(
                                                                h_residual_Zspectrum_N,
                                                                testParams,
                                                                args=(new_shifted[xval,yval], ppm_filtered), 
                                                                full_output = True,
                                                                maxfev = 900,
                                                                ftol =1E-9)
                    # Specify paramsets for peaks:
                    #TOFIX: why is the offset applied to each peak
                    #pk1 = [fit_params[0]]+list(fit_params[1:4])
                    #pk2 = [fit_params[0]]+list(fit_params[4:7])
                    #pk3 = [fit_params[0]]+list(fit_params[7:10])
                    #pk4 = [fit_params[0]]+list(fit_params[10:13])
                    #waterpk = [fit_params[0]]+list(fit_params[13:16]) 

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

                    water_amp[xval,yval] = fit_params[13]
                    water_width[xval,yval] = fit_params[14]            
                    water_pos[xval,yval] = fit_params[15]                
                  
                    fit_quality[xval,yval] = scipy.nansum(numpy.abs(new_shifted - h_zspectrum_N(fit_params,ppm_filtered-shift)))
                    fit_params_arr[xval,yval] = fit_params
                    ppm_corrected_arr[xval,yval] = ppm_filtered

        fitcount+=1 # increment fitcounter
    
    # Save the data as a structured array

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
    newstruct['water_A'] = water_amp
    newstruct['water_w'] = water_width
    newstruct['water_p'] = water_pos

    return {'':newstruct,'fit_quality':fit_quality}