"""

    SPIDIR

    Test sequence likelihood functions

"""

import sys, unittest, ctypes


sys.path.append("python")
import spidir

from test import *

from rasmus.common import *

class TestSeqLikelihood (unittest.TestCase):

    def setUp(self):
        pass        
        

    def test_calc_hky_seq_likelihood(self):

        bgfreq = [.258,.267,.266,.209]
        kappa = 1.59
        tree = treelib.parseNewick("((A:.1,B:.1):.1,((C:.1,D:.1):.2,E:.3):.1);")
        align = {"A": "CGCAGACAACTCCCCCGACCACACATAGTACGAAATCCTCAGCCGCTGCCGACTCCGACGCGCGGACTGTCCGGGTTCAGCGAGGCTTAAGAGAACGGCC",
                 "B": "CCCAAACAACTCCCCCGACCAGACATAGTACGAGATCCTCAGCCACTGGCGACTCGGACGCGCAGAGTGTCCCGCTTAAGCGAGGCTGCAGAGAACGGCC",
                 "C": "GGCCAGCAATTCCTCCGACCACGCATAGTACGAGATCGTCTGCCTCCTGCGAATCGGACGCGCAGAGTGTTCCGGTTAAGGGAGACTTCAGAGACCTGGC",
                 "D": "CGCTAACAATTCCCCCGACCACACTGAGTACGAGATACTCGGACTCCGGCGATCTCTACTCGCAGAGAGTCCCACTTAAGCGAGACTGACGAGCACGGGC",
                 "E": "ATTCTTCCACACCTGCGTGTTCGTCACGTATCAAATGCGGAGCCCACGTCCAATGGCACACGAACAGTCGGCCACGGAATCGCAGACTCGTTGACCAACG"}

        drawTree(tree)

        l = spidir.calc_seq_likelihood_hky(tree,
                                           mget(align, tree.leafNames()),
                                           bgfreq, kappa)

        l2 = spidir.find_ml_branch_lengths_hky(tree,
                                              mget(align, tree.leafNames()),
                                              bgfreq, kappa)

        drawTree(tree)
        
        print "log lk", l, l2
        self.assert_(l2 > l)


        

if __name__ == "__main__":
    unittest.main(testRunner=TestRunner())


