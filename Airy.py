#!/usr/bin/env python


"""
Airy.py

Implements fine granularity Airy function lookup.

Copyright(c) 2013 Jonathan D. Lettvin, All Rights Reserved"

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

__module__     = "Airy.py"
__author__     = "Jonathan D. Lettvin"
__copyright__  = """\
Copyright(C) 2011 Jonathan D. Lettvin, All Rights Reserved"""
__credits__    = [ "Jonathan D. Lettvin" ]
__license__    = "GPLv3"
__version__    = "0.0.1"
__maintainer__ = "Jonathan D. Lettvin"
__email__      = "jlettvin@gmail.com"
__contact__    = "jlettvin@gmail.com"
__status__     = "Demonstration"
__date__       = "20111027"

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
