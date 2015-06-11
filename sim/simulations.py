##Collection of code to simulate CEST sequences

import numpy
from scipy.signal import argrelextrema
import scipy.optimize
from scipy.integrate import ode

def setCESTdefaults(B1):
    satDur = 4000
    ti = 10  #time between sat pulse and acquisition
    tacq =10
    tpresat = 0
    accFactor = 72
    tinterfreq = 7000
    hardTheta = np.pi/9
    PElines =72
    m = PElines/accFactor
    dt = 1e-3
    delta = -1000.0 #freqency offset of the saturation pulse
    sequenceParams = [satDur, ti, tacq, tpresat, accFactor, tinterfreq, hardTheta, m, dt, delta]
    
    gamma = 2*pi*42.6e6 # rad/(s T)
    B0 = 7.0 #Tesla
    omega0 = gamma * B0
    omega1 = gamma * B1
    omegaWater = gamma * B0
    #omega1 = omegaWater + delta #frequency of sat pulse
    domegaSpecies = 4000 #rad/s
    omegaSpecies = omegaWater + domegaSpecies #chemical resonance frequency


    #Pool a is water, b is the species
    M0a = 1000.0
    M0b = M0a*1e-3*10
    T1a = 2.5 # seconds
    T2a = 0.6
    T1b = 1.0
    T2b = 0.2
    kb = 40.0
    ka = M0b/M0a*kb #Transfer rate of pool a to pool b, s^-1. We want a-->b = b-->a
    physicsVariables = [B0, omega1, domegaSpecies,  M0a, M0b, T1a, T2a, T1b, T2b, ka, kb]
    Mstart = array([0,0,0,0,M0a,M0b,1.])
    
    print 'sequenceParams = [satDur, ti, tacq, tpresat, accFactor, tinterfreq, hardTheta, m, dt, delta]'
    print 'physicsVariables = [B0, omega1, domegaSpecies,  M0a, M0b, T1a, T2a, T1b, T2b, ka, kb]'
    return sequenceParams, physicsVariables, Mstart


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
                        [0, -numpy.sin(phi), 0, numpy.cos(phi), 0, 0, 0],
                        [0, 0, 0, 0, 1, 0, 0],
                        [0, 0, 0, 0, 0, 1, 0],
                        [0, 0, 0, 0, 0, 0, 1],])


def ABtwoPool(dt, T1a=numpy.inf, T2a=numpy.inf, T1b = numpy.inf, T2b = numpy.inf, M0a=1000, M0b = 1.0, domega=0):
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
                     [0, 0, 0, 0, 0, 0, 1 ]])
    return numpy.dot(A, zrot(phi)),B

def freePrecessTwoPool(Mresult, t, A_fp, B_fp):
    if t > 0:
        Mresult_fp = numpy.empty((t+1,7))
        Mresult_fp[0,:] = np.array(Mresult)[-1,:]
        for i in range(1, t+1):
            Mresult_fp[i,:] = numpy.dot(A_fp, Mresult_fp[i-1,:]) + B_fp
        return numpy.concatenate((Mresult, Mresult_fp[1:-1]), 0)
    else:
        return Mresult



def cestSequence(Mstart, physicsVariables, sequenceParams):
    ## Takes a starting magnetization, returns a history of the magnetization over the course of a CEST FLASH sequence
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
    
    Mresult = np.empty((int(satDur),7))
    Mresult[0,:], tresult = Mstart, [0]
 
    r = scipy.integrate.ode(dMdtTwoPool)
    r = r.set_integrator('dopri5')
    r = r.set_initial_value(Mstart, t=0)

    t = 0.0
    idx = 1
    while r.successful() and idx < satDur:
        Mresult[idx,:] = r.integrate(r.t + dt)
        t+= dt
        idx += 1
        
    if delta > 0.:
        Mresult[-1,0] = -Mresult[-1,0]
        Mresult[-1,2] = -Mresult[-1,2]
    
    ##################    END OF SATURATION PULSE  #####################################

    
    ##################   IMAGING SEQUENCE     ##########################################
    signals = []
    dt = 0.001
    A_fp, B_fp = ABtwoPool(dt, T1a=T1a, T2a=T2a, T1b = T1b, T2b = T2b, M0a=M0a, M0b = M0b, domega=0)
    Mresult = freePrecessTwoPool(Mresult, ti, A_fp, B_fp)## between sat pulse and aquisition pulse
    #for mFactor in range(m):  ##  Acquisition, including hard flip
    for i in range(accFactor):
        Mresult = numpy.concatenate((Mresult, [np.dot(yrot(hardTheta), Mresult[-1])]))
        signals.append(numpy.sqrt(Mresult[-1,0]**2 + Mresult[-1,2]**2))
        Mresult = freePrecessTwoPool(Mresult, tacq, A_fp, B_fp)
        Mresult[-1][0:4] = [0,0,0,0] ## Spoiler Gradient
    
    
    Mresult = freePrecessTwoPool(Mresult, tpresat, A_fp, B_fp)     ## tPresat - Does this exist? 
    
    Mresult = freePrecessTwoPool(Mresult, tinterfreq, A_fp, B_fp)     ## after acquisition, before the next frequency offset
    #################     END OF IMAGING SEQUENCE     ####################################
    
    return Mresult, signals

def Zspectrum(freqs, Mstart):
    signals = []
    Mresults = []
    for freq in freqs:
        print(freq)
        sequenceParams[-1] = freq
        Mresult, signal = cestSequence(Mstart, physicsVariables, sequenceParams)
        Mstart = Mresult[-1,:]
        Mresults.append(Mresult)
        signals.append(signal)
    return Mresult, signals
