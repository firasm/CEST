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
