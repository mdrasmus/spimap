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

rplot_set_viewer("display")


class TestHKY (unittest.TestCase):
    
    def test_hky(self):
        """general test"""
        
        bgfreq = [.3, .2, .3, .2]
        kappa = 2.0
        t = 0.2

        pprint(spidir.make_hky_matrix(bgfreq, kappa, t))

    def test_hky_deriv(self):
        """general test"""
        
        bgfreq = [.3, .2, .3, .2]
        kappa = 2.0
        i = random.randint(0, 3)
        j = random.randint(0, 3)
        
        x = list(frange(0, 1.0, .01))
        y = [spidir.make_hky_matrix(bgfreq, kappa, t)[i][j]
             for t in x]
        dy = [spidir.make_hky_deriv_matrix(bgfreq, kappa, t)[i][j]
              for t in x]
        dy2 = [(spidir.make_hky_matrix(bgfreq, kappa, t+.01)[i][j] -
                spidir.make_hky_matrix(bgfreq, kappa, t)[i][j]) / .01
               for t in x]

        prep_dir("test/output/hky")
                
        rplot_start("test/output/hky/hky_deriv.pdf")
        rplot("plot", x, y, t="l", ylim=[min(dy + y), max(dy + y)])
        rp.lines(x, dy, col="red")
        rp.lines(x, dy2, col="blue")
        rplot_end(True)


    def test_hky_deriv2(self):
        """general test"""
        
        bgfreq = [.3, .2, .3, .2]
        kappa = 2.0
        i = random.randint(0, 3)
        j = random.randint(0, 3)
        
        x = list(frange(0, 1.0, .01))
        y = [spidir.make_hky_deriv_matrix(bgfreq, kappa, t)[i][j]
             for t in x]
        dy = [spidir.make_hky_deriv2_matrix(bgfreq, kappa, t)[i][j]
              for t in x]
        dy2 = [(spidir.make_hky_deriv_matrix(bgfreq, kappa, t+.01)[i][j] -
                spidir.make_hky_deriv_matrix(bgfreq, kappa, t)[i][j]) / .01
               for t in x]

        prep_dir("test/output/hky2")
                
        rplot_start("test/output/hky2/hky_deriv2.pdf")
        rplot("plot", x, y, t="l", ylim=[min(dy2 + dy + y), max(dy2 + dy + y)])
        rp.lines(x, dy, col="red")
        rp.lines(x, dy2, col="blue")
        rplot_end(True)


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
