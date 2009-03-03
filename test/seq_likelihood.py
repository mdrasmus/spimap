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
        

    def _test_calc_hky_seq_likelihood(self):

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


    def test_branch_likelihood_hky(self):
        """Test likelihood function"""

        # params
        bgfreq = [.258,.267,.266,.209]
        kappa = 1.59
        seqlen = 1000
        div = .1
        t0 = 0
        t1 = .5
        step = .001

        # prep probabilities
        probs1 = []
        probs2 = []
        for i in xrange(seqlen):
            if random.random() < div:
                for j in xrange(4):
                    probs1.append(.01)
                    probs2.append(.01)
                k = random.randint(1, 4)
                probs1[-k] = 1.0
                k = (k % 4) + 1
                probs2[-k] = 1.0
            else:
                for j in xrange(4):
                    probs1.append(.2)
                    probs2.append(.2)
                k = random.randint(1, 4)
                probs1[-k] = 1.0
                probs2[-k] = 1.0

        # estimate MLE
        mle = spidir.mle_distance_hky(probs1, probs2, seqlen, bgfreq, kappa,
                                      t0, t1, step)

        x = list(frange(0, .5, .01))
        y = [spidir.branch_likelihood_hky(probs1, probs2, seqlen,
                                           bgfreq, kappa, t)
             for t in x]

        top = x[argmax(y)]

        print div, top, mle


        prep_dir("test/output/branch_likelihood/")
        
        rplot_start("test/output/branch_likelihood/branch_function.pdf")
        rplot("plot", x, y, t="l",
              xlab="distance",
              ylab="likelihood")
        #rp.lines([div, div], [-1e300, 0], col="green")
        rp.lines([top, top], [-1e300, 0], col="blue")        
        rp.lines([mle, mle], [-1e300, 0], col="red")
        rplot_end(True)


        x = list(frange(0, .5, .01))
        dy = [spidir.derive_branch_likelihood_hky(probs1, probs2, seqlen,
                                                  bgfreq, kappa, t)
              for t in x]
        dy2 = [(spidir.branch_likelihood_hky(probs1, probs2, seqlen,
                                             bgfreq, kappa, t+.01)  -
                spidir.branch_likelihood_hky(probs1, probs2, seqlen,
                                             bgfreq, kappa, t))
               for t in x]

        rplot_start("test/output/branch_likelihood/deriv_branch_function.pdf")
        rplot("plot", x, dy2, t="l",
              xlab="distance",
              ylab="d/dt likelihood")
        rp.lines(x, dy, col="grey")
        #rp.lines([div, div], [-1e300, 0], col="green")
        rp.lines([top, top], [-1e300, 0], col="blue")        
        rp.lines([mle, mle], [-1e300, 0], col="red")
        rplot_end(True)
        

if __name__ == "__main__":
    unittest.main(testRunner=TestRunner())


