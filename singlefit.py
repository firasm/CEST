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