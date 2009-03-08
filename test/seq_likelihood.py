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


    def test_branch_likelihood_hky(self):
        """Test likelihood function"""

        # params
        bgfreq = [.2,.3,.3,.2]
        kappa = 1.59
        seqlen = 100
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
                                      t0, t1)

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

        #================================
        # 1st derivative

        x = list(frange(0, .5, .01))
        dy = [spidir.branch_likelihood_hky_deriv(probs1, probs2, seqlen,
                                                 bgfreq, kappa, t)
              for t in x]
        dy2 = [(spidir.branch_likelihood_hky(probs1, probs2, seqlen,
                                             bgfreq, kappa, t+.01)  -
                spidir.branch_likelihood_hky(probs1, probs2, seqlen,
                                             bgfreq, kappa, t)) / .01
               for t in x]
        
        
        rplot_start("test/output/branch_likelihood/deriv_branch_function.pdf")
        rplot("plot", x, dy2, t="l",
              xlab="distance",
              ylab="d/dt likelihood")
        rp.lines([min(x), max(x)], [0,0], col="black")
        rp.lines(x, dy, col="grey")
        #rp.lines([div, div], [-1e300, 0], col="green")
        rp.lines([top, top], [-1e300, 0], col="blue")        
        rp.lines([mle, mle], [-1e300, 0], col="red")
        rplot_end(True)


        #=============================
        # 2nd derivative

        x = list(frange(0, .5, .01))
        d2y = [spidir.branch_likelihood_hky_deriv2(probs1, probs2, seqlen,
                                                   bgfreq, kappa, t)
               for t in x]
        d2y2 = [(spidir.branch_likelihood_hky_deriv(probs1, probs2, seqlen,
                                                    bgfreq, kappa, t+.01)  -
                 spidir.branch_likelihood_hky_deriv(probs1, probs2, seqlen,
                                                    bgfreq, kappa, t)) / .01
                for t in x]
        
        
        rplot_start("test/output/branch_likelihood/deriv2_branch_function.pdf")
        rplot("plot", x, d2y2, t="l",
              xlab="distance",
              ylab="d^2/dt^2 likelihood")
        rp.lines([min(x), max(x)], [0,0], col="black")
        rp.lines(x, d2y, col="grey")
        #rp.lines([div, div], [-1e300, 0], col="green")
        #rp.lines([top, top], [-1e300, 0], col="blue")        
        #rp.lines([mle, mle], [-1e300, 0], col="red")
        rplot_end(True)


    def _test_calc_lktable_row(self):
        """test the function CalcLktbaleRow"""

        def branchlk(probs1, probs2, seqlen, bgfreq, kappa, t):

            model1 = spidir.make_hky_matrix(bgfreq, kappa, t)
            model2 = spidir.make_hky_matrix(bgfreq, kappa, 0)
            
            logl = 0.0
            for j in xrange(seqlen):
                s = sum(bgfreq[k] *
                        sum(model1[k][x] * probs1[4*j+x] for x in xrange(4)) *
                        sum(model2[k][y] * probs2[4*j+y] for y in xrange(4))
                        for k in xrange(4))
                logl += safelog(s, e)
                

            return logl


        def dbranchlk(probs1, probs2, seqlen, bgfreq, kappa, t):

            model1 = spidir.make_hky_matrix(bgfreq, kappa, t)
            model2 = spidir.make_hky_matrix(bgfreq, kappa, 0.0)

            dmodel1 = spidir.make_hky_deriv_matrix(bgfreq, kappa, t)
            dmodel2 = spidir.make_hky_deriv_matrix(bgfreq, kappa, 0.0)

            
            logl = 0.0
            for j in xrange(seqlen):
                ds = sum(bgfreq[k] *
                        sum(dmodel1[k][x] * probs1[4*j+x] for x in xrange(4)) *
                        sum(model2[k][y] * probs2[4*j+y] for y in xrange(4))
                        for k in xrange(4))
                
                s = sum(bgfreq[k] *
                        sum(model1[k][x] * probs1[4*j+x] for x in xrange(4)) *
                        sum(model2[k][y] * probs2[4*j+y] for y in xrange(4))
                        for k in xrange(4))

                logl += safediv(ds, s, INF)                
            return logl
        
        def d2branchlk(probs1, probs2, seqlen, bgfreq, kappa, t):

            model1 = spidir.make_hky_matrix(bgfreq, kappa, t)
            model2 = spidir.make_hky_matrix(bgfreq, kappa, 0.0)

            dmodel1 = spidir.make_hky_deriv_matrix(bgfreq, kappa, t)
            dmodel2 = spidir.make_hky_deriv_matrix(bgfreq, kappa, 0.0)

            d2model1 = spidir.make_hky_deriv2_matrix(bgfreq, kappa, t)
            d2model2 = spidir.make_hky_deriv2_matrix(bgfreq, kappa, 0.0)
            
            
            logl = 0.0
            for j in xrange(seqlen):
                g = sum(bgfreq[k] *
                        sum(model1[k][x] * probs1[4*j+x] for x in xrange(4)) *
                        sum(model2[k][y] * probs2[4*j+y] for y in xrange(4))
                        for k in xrange(4))

                dg = sum(bgfreq[k] *
                        sum(dmodel1[k][x] * probs1[4*j+x] for x in xrange(4)) *
                        sum(model2[k][y] * probs2[4*j+y] for y in xrange(4))
                        for k in xrange(4))                

                d2g = sum(bgfreq[k] *
                        sum(d2model1[k][x] * probs1[4*j+x] for x in xrange(4)) *
                        sum(model2[k][y] * probs2[4*j+y] for y in xrange(4))
                        for k in xrange(4))                


                logl += - safediv(dg*dg, g*g, INF) + \
                        safediv(d2g, g, INF)
            return logl


        bgfreq = [.25,.25,.25,.25]
        kappa = 1.59
        seqlen = 100

        # prep probabilities
        div = .1
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

        #probs1 = [0.0, 0.0, 1.0, 0.0] + \
        #         [1.0, 0.0, 0.0, 0.0] * 5
        #probs2 = [0.0, 0.0, 0.0, 1.0] + \
        #         [1.0, 0.0, 0.0, 0.0] * 5
        #seqlen = len(probs1) / 4

        x = list(frange(0, 1.0, .01))

        y = [spidir.branch_likelihood_hky(probs1, probs2, seqlen,
                                          bgfreq, kappa, t)
             for t in x]
        y2 = [branchlk(probs1, probs2, seqlen,
                       bgfreq, kappa, t)
             for t in x]

        prep_dir("test/output/branch_likelihood_simple/")
        
        rplot_start("test/output/branch_likelihood_simple/cmp_c_py.pdf")
        rplot("plot", x, y, t="l")
        rp.lines(x, y2, t="l", col="red")
        rplot_end(True)


        x = list(frange(0, 1.0, .01))

        y = [spidir.branch_likelihood_hky_deriv(probs1, probs2, seqlen,
                                                bgfreq, kappa, t)
             for t in x]
        y2 = [dbranchlk(probs1, probs2, seqlen,
                        bgfreq, kappa, t)
              for t in x]

        rplot_start("test/output/branch_likelihood_simple/cmp_c_py_deriv.pdf")
        rplot("plot", x, y, t="l")
        rp.lines(x, y2, t="l", col="red")
        rplot_end(True)


        x = list(frange(0, 1.0, .01))

        y = [spidir.branch_likelihood_hky_deriv2(probs1, probs2, seqlen,
                                                 bgfreq, kappa, t)
             for t in x]
        #y = [(spidir.branch_likelihood_hky_deriv(probs1, probs2, seqlen,
        #                                         bgfreq, kappa, t+.01) -
        #      spidir.branch_likelihood_hky_deriv(probs1, probs2, seqlen,
        #                                         bgfreq, kappa, t)) / .01
        #     for t in x]
        y2 = [d2branchlk(probs1, probs2, seqlen,
                          bgfreq, kappa, t)
              for t in x]

        rplot_start("test/output/branch_likelihood_simple/cmp_c_py_deriv2.pdf")
        rplot("plot", x, y, t="l")
        rp.lines(x, y2, t="l", col="red")
        rplot_end(True)

if __name__ == "__main__":
    unittest.main(testRunner=TestRunner())


