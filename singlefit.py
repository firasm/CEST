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

def h_zspectrum_N(params,freqs):
    ''' Updated Zspectrum N function to now require a structured array of params'''

    paramStructuredArray = cest.analysis.h_convertBetweenStructArrays(params,toType = 'struct')

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

        newpeak = cest.analysis.h_lorentzian(A,w,p,freqs)

        if numpy.isnan(numpy.sum(newpeak)):
            print(A,w,p)
            print('struct: {}'.format(paramStructuredArray['A{}'.format(i)],paramStructuredArray['w{}'.format(i)],paramStructuredArray['p{}'.format(i)]))
        arr += newpeak

    A = paramStructuredArray['water_A'][0]+1E-8
    w = paramStructuredArray['water_w'][0]+1E-9
    p = paramStructuredArray['water_p'][0]+1E-10

    arr += cest.analysis.h_lorentzian(A,w,p,freqs)
    return arr + shift

def h_residual_Zspectrum_N(params, y_data, w):

    # Define the penalty function:
    #def old penaltyfn(x,centre=0., scale=1., trough_width=1., steepness=2.):
    #    return scale*((centre-x)/trough_width)**(2*steepness)

    params = cest.analysis.h_convertBetweenStructArrays(params,toType = 'struct')

    def penaltyfn(x, centre, trough = 0.1, scale=1E-4,hardlimit = -50):

        pen = scale*((centre-x)/trough)**2*2
        #print(pen)
        return 0
        return pen

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

    params = cest.analysis.h_convertBetweenStructArrays(params,toType = 'array')

    return numpy.abs(y_data - cest.analysis.h_zspectrum_N(params,w)) + penalty

def fit_water_peak(data,offset_freqs,allParams = False):

    """
    Fits a Super Lorentzian to the data, and 
    returns the water offset frequency

    """

    # First get rid of all attempts to pass in bad data. If there are any nans in there, toss it.

    if numpy.isnan(numpy.sum(data)):
        return numpy.nan

    else:

        # Presumably the maximum signal will be at 
        # some large frequency offset, so let's just
        # use that as our baseline parameter

        params = numpy.zeros((1), dtype=[('offset', 'float64'),('water_A', 'float64'), 
                                         ('water_w', 'float64'),('water_p', 'float64')])

        # Also, normalize the data so the fit is easier

        params['offset'] = numpy.max(data)
        params['water_A'] = -1.
        params['water_w'] = 0.6
        params['water_p'] = offset_freqs[numpy.argmin(data)]

        params = cest.analysis.h_convertBetweenStructArrays(params, toType = 'array')

        fit_params,cov,infodict,mesg,ier = scipy.optimize.leastsq(h_residual_Zspectrum_N,
                                                                    params,
                                                                    args=(data, offset_freqs),
                                                                    full_output = True,
                                                                    maxfev = 200)
        
        if allParams:
            return fit_params
        else:
            return cest.analysis.h_convertBetweenStructArrays(params, toType = 'struct')['water_p']


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

################################# Start of fitting portion

# Function inputs
scn_to_analyse = scn.shortdirname
xval =27
yval = 40
fitrounds = 1
###################################################
#Function start

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
        if shift <0: 
            s_shifted_back = scipy.interp(ppm_filtered, ppm_filtered+shift/2, tmp)
        else:
            s_shifted_back = scipy.interp(ppm_filtered, ppm_filtered-shift, tmp)


        new_shifted = s_shifted_back       
    else:
        print(shift)
        pass            

    if fitcount>0: # Use parameters from last fit 
        testParams = cest.analysis.h_convertBetweenStructArrays(newstruct,toType='array')
    else: # Get initial starting parameters
        testParams = cest.analysis.get_neighbours_starting()
        originaltestParams = cest.analysis.h_convertBetweenStructArrays(testParams,toType='struct').copy()       
        testParams = cest.analysis.h_convertBetweenStructArrays(testParams,toType = 'array')

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

fitresults = {'fit_params':newstruct, 
                'data': [ppm_filtered, new_shifted],
                'fitdata': [freqs, h_zspectrum_N(newstruct, freqs)],
                'fit_quality': fit_quality}

print('Fit Quality is: {0}'.format(fit_quality))