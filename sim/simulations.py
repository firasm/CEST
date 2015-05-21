##Collection of code to simulate CEST sequences

import numpy

def testThis():

    return 2*100

def xrot(phi):
    return numpy.array([[1,0,0,0,0,0,0],
                        [0,1,0,0,0,0,0],
                        [0,0,cos(phi),0,sin(phi),0,0],
                        [0,0,0,cos(phi),0,sin(phi),0],
                        [0,0,-sin(phi),0,cos(phi),0,0],
                        [0,0,0,-sin(phi),0,cos(phi),0],
                        [0,0,0,0,0,0,1],])

def yrot(phi):
    return numpy.array([[cos(phi),0,0,0,-sin(phi),0,0],
                        [0,cos(phi),0,0,0,-sin(phi),0],
                        [0,0,1,0,0,0,0],
                        [0,0,0,1,0,0,0],
                        [sin(phi),0,0,0,cos(phi),0,0],
                        [0,sin(phi),0,0,0,cos(phi),0],
                        [0,0,0,0,0,0,1]])

def zrot(phi):
    return numpy.array([[numpy.cos(phi), 0, numpy.sin(phi), 0, 0, 0, 0],
                        [0, numpy.cos(phi), 0, numpy.sin(phi), 0, 0, 0],
                        [-numpy.sin(phi), 0, numpy.cos(phi), 0, 0, 0, 0],
                        [0, -numpy.sin(phi), 0, numpy.cos(phi), 1, 0, 0],
                        [0, 0, 0, 0, 1, 0, 0],
                        [0, 0, 0, 0, 0, 1, 0],
                        [0, 0, 0, 0, 0, 0, 1],])

def freeprecessTwoPool(dt, relaxationTimes = [T1a, T2a, T1b, T2b], poolMagnitudes = [1000.0, 1.0], domega=domega):
    ''' return the A matrix and B vector for the dM/dt magnetization evolution '''
    phi = domega*dt	 # Resonant precession, radians.
    E1a = numpy.exp(-dt/T1a)
    E1b = numpy.exp(-dt/T1b)
    E2a = numpy.exp(-dt/T2a)
    E2b = numpy.exp(-dt/T2b)
    
    B = numpy.array([0, 0, 0, 0, M0a*(1. - E1a), M0b*(1. - E1b), 0])

    A = numpy.array([[E2a, 0, 0, 0, 0, 0, 0],
                     [0, E2b, 0, 0, 0, 0, 0],
                     [0, 0, E2a, 0, 0, 0, 0],
                     [0, 0, 0, E2b, 0, 0, 0],
                     [0, 0, 0, 0, E1a, 0, 0],
                     [0, 0, 0, 0, 0, E1b, 0],
                     [0, 0, 0, 0, 0, 0, 0 ]])
    return numpy.dot(A, zrot(phi)),B

## Takes a starting magnetization, returns a history of the magnetization over the course of a CEST FLASH sequence

def cestSequence(Mstart, physicsVariables, sequenceParams):
    
    [satDur, ti, tacq, tpresat, accFactor, tinterfreq, hardTheta, m, dt, delta] = sequenceParams
    [B0, omega1, domegaSpecies,  M0a, M0b, T1a, T2a, T1b, T2b, ka, kb] = physicsVariables
    
    
    ################    SATURATION PULSE    ##################################    
    def dMdtTwoPool(t, M_vec, M0a = M0a, M0b = M0b, T1a = T1a, T2a = T2a, T1b = T1b, T2b = T2b,
                ka = ka, kb = kb, domegaa=-delta, domegab=-delta-domegaSpecies, omega1=omega1):

        A = numpy.array([[-1./T2a - ka, kb,  domegaa, 0, 0, 0, 0],
                         [ka, -1./T2b - kb, 0, domegab, 0, 0 , 0],
                         [-domegaa, 0, -1./T2a-ka, kb, omega1, 0, 0],
                         [0, -domegab, ka, -1./T2b-kb, 0, omega1, 0],
                         [0, 0, -omega1, 0, -1./T1a - ka, kb, M0a/T1a],
                         [0, 0, 0, -omega1, ka, -1./T1b - kb, M0b/T1b],
                         [0, 0, 0, 0, 0, 0,0]])
        
        return numpy.dot(A,M_vec)
    
    
    tresult = numpy.empty(satDur)
    Mresult = numpy.empty((satDur,7))
    Mresult[0,:] = Mstart
    
    r = scipy.integrate.ode(dMdtTwoPool)
    r = r.set_integrator('dopri5')
    r = r.set_initial_value(Mstart, t=0)
    
    idx = 1
    while r.successful() and r.t < satDur and idx<satDur:
        Mresult[idx,:] = r.integrate(r.t+dt)
        tresult[idx]=r.t
        idx += 1
        
    #signal = numpy.sqrt(Mresult[-1,0]**2 + Mresult[-1,2]**2)
    #signal = numpy.sqrt(Mresult[-1,4]**2)

    ##################    END OF SATURATION PULSE  #####################################
   
    ##################   IMAGING SEQUENCE     ##########################################
    A_fp, B_fp = ABtwoPool(dt, T1a=T1a, T2a=T2a, T1b = T1b, T2b = T2b, M0a=M0a, M0b = M0b, domega=0)
    
    Mresult = freePrecessTwoPool(Mresult, ti, A_fp, B_fp)## between sat pulse and aquisition pulse

    for mFactor in range(m):  ##  Acquisition, including hard pulse 
        for accFactor in range(1, accFactor+1):
            Mresult = numpy.concatenate((Mresult, [np.dot(yrot(hardTheta), Mresult[-1])]))
            signal = numpy.sqrt(Mresult[-1,0]**2 + Mresult[-1,2]**2)
            Mresult = freePrecessTwoPool(Mresult, tacq, A_fp, B_fp)
            Mresult[-1][0:4] = [0,0,0,0] ## Spoiler Gradient
    
    
    Mresult = freePrecessTwoPool(Mresult, tpresat, A_fp, B_fp)     ## tPresat - Does this exist? 
    
    Mresult = freePrecessTwoPool(Mresult, tinterfreq, A_fp, B_fp)     ## after acquisition, before the next frequency offset
    
    
    #################     END OF IMAGING SEQUENCE     ####################################
    
    return Mresult, signal
