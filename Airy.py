#!/usr/bin/env python

"""Airy.py implements fine granularity Airy function lookup."""

import string
from scipy import arange, pi, set_printoptions
from scipy.special import j1
from optparse import OptionParser

class Airy(object):
    """
    The Airy class maintains an abstract fine-granularity Airy lookup table.
    http://en.wikipedia.org/wiki/Angular_resolution
    http://en.wikipedia.org/wiki/Airy_disk
    """
    u0 = 3.8317059702075125 # u0: First zero of j1(u)/u
    r0 = 1.2196698912665045 # r0: Resolution limit u0/r0 == pi, u0/pi == r0
    point, zeros, peaks, spike = [], [], [], [] # Initialized singletons

    def element(this,i,r):
        u = r * Airy.u0 * 1e+4
        sqrtI = 1.0 if u == 0.0 else 2.0*j1(u)/u
        return {'index': i, 'radius': r, 'parameter': u, 'amplitude': sqrtI}

    def __init__(this, **kwargs):
        if not Airy.point:
            delta    = kwargs.get(  'delta', 1e-7)
            epsilon  = kwargs.get('epsilon', 1e-3)
            valn, val0 = this.element(-2,2*delta), this.element(-1,1*delta)
            terminal, avoid_neighbors = False, False
            for index, radius in enumerate(arange(0.0, 1e-3, delta)):
                valp = this.element(index,radius)
                Airy.point += [valp,]
                vn, v0, vp = (val['amplitude'] for val in (valn,val0,valp))
                if vn * vp < 0.0:
                    if avoid_neighbors:
                        avoid_neighbors = False
                    else:
                        Airy.zeros += [val0,]
                        Airy.spike += [(val0['parameter'],abs(Airy.peaks[-1]['amplitude'])),]
                        avoid_neighbors = True
                    "Find first zero past below epsilon peak"
                    if terminal: break
                elif abs(vn) <= abs(v0) >= abs(vp):
                    Airy.peaks += [val0,]
                    "Find a below epsilon peak"
                    if abs(val0['amplitude']) < epsilon: terminal = True
                valn, val0 = val0, valp

Airy() # Generate singleton

if __name__ == "__main__":

    def test1(): print len(Airy.zeros), "zeros"
    def test2(): print Airy.zeros[0]['parameter'], Airy.point[1]['amplitude']
    def test3(): print string.join(["%1.16e, %1.16e" % (u,am) for (u,am) in Airy.spike], '\n')
    def test4(): print Airy.zeros[0]['parameter']/Airy.r0

    def runtest(key, val): print '\t%s' %(key); val()

    tests = {"zerocount": test1, "firstzero": test2, "spikes": test3, "u0/r0": test4}

    parser = OptionParser()
    parser.add_option("-d", "--delta",      default=1e-7,  help="radial step size")
    parser.add_option("-e", "--epsilon",    default=1e-3,  help="lower cutoff for kernel")
    parser.add_option("-l", "--linewidth",  default=1000,  help="print width")
    parser.add_option("-p", "--precision",  default=2,     help="print precision")
    (opts, args) = parser.parse_args()
    kwargs = vars(opts)

    set_printoptions(precision=kwargs.get('precision',2), linewidth=kwargs.get('linewidth', 1000))
    [runtest(key,val) for (key,val) in tests.iteritems()]
