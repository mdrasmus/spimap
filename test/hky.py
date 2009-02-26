"""

   test HKY matrix

"""


import sys, unittest, ctypes
from pprint import pprint
from math import *

sys.path.append("python")
import spidir

from test import *

from rasmus.common import *


class TestHKY (unittest.TestCase):
    
    def test_hky(self):
        """general test"""
        
        bgfreq = [.3, .2, .3, .2]
        kappa = 2.0
        t = 0.2

        pprint(spidir.make_hky_matrix(bgfreq, kappa, t))


    def test_JC(self):
        """test equivalence to JC"""
        
        bgfreq = [.25, .25, .25, .25]
        kappa = 1.0

        for t in frange(0, 1.0, .1):
            mat = spidir.make_hky_matrix(bgfreq, kappa, t)

            a = 1/3.0
            r = (1/4.0)*(1 + 3*exp(-4*a*t))
            s = (1/4.0)*(1 - exp(-4*a*t))

            mat2 = [[r, s, s, s],
                    [s, r, s, s],
                    [s, s, r, s],
                    [s, s, s, r]]

            for i in xrange(4):
                for j in xrange(4):
                    fequal(mat[i][j], mat2[i][j])




if __name__ == "__main__":
    unittest.main(testRunner=TestRunner())
