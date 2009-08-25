"""

    SPIDIR

    Test birth death functions - estimating lambda and mu

"""

import sys, unittest, ctypes


sys.path.append("python")
import spidir
from spidir import calc_birth_death_prior as c_calcBirthDeathPrior

from test import *

from rasmus.common import *
from rasmus import stats
from rasmus.bio import phylo, birthdeath


if os.system("which xpdf 2> /dev/null") != 0:
    rplot_set_viewer("display")

def prob_birth_death_same(a, b, t, lam):

    al = lam*t / (1 + lam*t)
    
    return sum(choose(a, j) * choose(a+b-j-1, a-1) *
               al**(a+b-2*j) * (1-2*al)**j
               for j in range(min(a,b)+1))


#=============================================================================


class TestBirthDeath (unittest.TestCase):

    def setUp(self):
        pass

    def test_birth_death_large(self):
        
        t = 0.0392320007
        birth = 0.400400013
        death = 0.400000006

        print "large"
        print birthdeath.prob_birth_death(60, 71, t, birth, death)
        print spidir.birthDeathCounts2(60, 71, t, birth, death)
        print spidir.birthDeathCounts(60, 71, t, birth, death)


        t = 0.067008
        birth = 0.010100
        death = 0.010000 
        print spidir.birthDeathCounts(51, 51, t, birth, death)
        print spidir.birthDeathCounts2(51, 51, t, birth, death)
        print spidir.birthDeathCounts(51, 51, t, birth, death)

    def test_birth_death_same(self):

        counts = list(range(1, 10))
        birth = .2
        death = .1
        lam = birth
        t = 10.4

        for a in range(1, 10):
            for b in range(0, 10):                
                l = spidir.birthDeathCounts(a, b, t, birth, death)
                #l2 = spidir.birthDeathCounts2(a, b, t, birth, death*.999)
                #l2 = prob_birth_death_same(a, b, t, lam)
                l2 = l
                l3 = birthdeath.prob_birth_death(a, b, t, birth, death)
                
                print a, b, l, l2, l3
                #fequal(l, l2, .01)
                fequal(l, l3, .01)

            #fequal(sum(spidir.birthDeathCounts(a, b, t, birth, death)
            #           for b in xrange(0, 20)),
            #       1.0, .001)


    def test_birth_death(self):
        
        counts = list(range(1, 10))
        birth = 0.01
        death = 0.02
        t = 1.0

        for a in range(1, 10):
            for b in range(0, 10):
                l = spidir.birthDeathCounts(a, b, t, birth, death)
                #l2 = spidir.birthDeathCounts2(a, b, t, birth, death)
                l2 = l
                l3 = birthdeath.prob_birth_death(a, b, t, birth, death)
                
                print a, b, l, l2, l3
                #fequal(l, l2, .1)
                fequal(l, l3, .01)

            fequal(sum(spidir.birthDeathCounts(a, b, t, birth, death)
                       for b in xrange(0, 20)),
                   1.0, .001)


    def test_sample_birth_times(self):

        stree = read_tree("test/data/flies.stree")
        
        birth = 8.1
        death = 2.1 * 1.0001

        print spidir.birthDeathCount(4, .005, birth, death)
        print spidir.birthDeathCounts(1, 4, .005, birth, death)
        

        
        maxgenes = [5, 10, 20, 40, 50]
        tic("time")
        counts = [2, 2, 3, 3, 3, 2, 1, 1, 1, 3, 3, 2]
        for maxgene in maxgenes:
            logl = spidir.birth_death_tree_counts(stree, counts,
                                                     birth, death,
                                                     maxgene=maxgene,
                                                     rootgene=1)
            print maxgene, logl
        toc()

        counts = [1] * 12
        for maxgene in maxgenes:
            logl = spidir.birth_death_tree_counts(stree, counts,
                                                     birth, death,
                                                     maxgene=maxgene,
                                                     rootgene=1)
            print maxgene, logl

                                          
        counts = [4] * 12        
        for maxgene in maxgenes:
            logl = spidir.birth_death_tree_counts(stree, counts,
                                                     birth, death,
                                                     maxgene=maxgene,
                                                     rootgene=1)
            print maxgene, logl


    def test_sample_birth_rates(self):

        prep_dir("test/output/birth_death_est")

        stree = read_tree("test/data/flies.stree")
    
        maxgene = 20
        counts = [2, 2, 3, 3, 3, 2, 0, 0, 0, 3, 3, 2]
        #counts = [5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
        
        rates = list(frange(.1, 20, 1))
        
        top = -INF

        x = []
        y = []
        z = []
        
        for birth in rates:
            for death in rates:
                l = spidir.birth_death_tree_counts(stree, counts,
                                                   birth, death,
                                                   maxgene=maxgene,
                                                   rootgene=1)
                x.append(birth)
                y.append(death)
                z.append(l)
                
                if l > top:
                    top = l
                    b = birth
                    d = death
                print birth, death, l

        print "best:", b, d, top

        write_delim("test/output/birth_death_est/rates.txt", zip(x, y, z))

        mat = [mget(z, range(ind, ind+len(rates)))
               for ind in xrange(0, len(z), len(rates))]
        #mat = map2(log, mat)

        heatmap(mat, rlabels=map(str, rates), clabels=map(str, rates),
                display=True, xmargin=100, ymargin=100)


    def test_sample_birth_rates_ml(self):

        stree = read_tree("test/data/flies.stree")
    
        maxgene = 20
        counts = [[2, 2, 3, 3, 3, 2, 0, 0, 0, 3, 3, 2]]
        birth0 = .1
        death0 = .1
        step = 1.0

        opt = spidir.birth_death_counts_ml_alloc(stree, counts,
                                                 birth0, death0, step,
                                                 maxgene=maxgene)

        for i in range(300):
            status, size, (b,d) = spidir.birth_death_counts_ml_iter(opt)
            print status, size, (b, d)
            if status == 0:
                break

        print "ml b,d =", (b,d)
    
        spidir.birth_death_counts_ml_free(opt)


    def _test_sample_birth_rates_forest(self):

        prep_dir("test/output/birth_death_est")

        stree = read_tree("test/data/flies.stree")
    
        maxgene = 20
        counts = [[2, 2, 3, 3, 3, 2, 0, 0, 0, 3, 3, 2]] * 20
                  
        
        rates = list(frange(.1, 20, 2))
        
        top = -INF

        x = []
        y = []
        z = []
        
        for birth in rates:
            for death in rates:
                l = spidir.birth_death_forest_counts(stree, counts,
                                                     birth, death,
                                                     maxgene=maxgene,
                                                     rootgene=1)
                x.append(birth)
                y.append(death)
                z.append(l)
                
                if l > top:
                    top = l
                    b = birth
                    d = death
                print birth, death, l

        print "best:", b, d, top

        write_delim("test/output/birth_death_est/rates2.txt", zip(x, y, z))

        mat = [mget(z, range(ind, ind+len(rates)))
               for ind in xrange(0, len(z), len(rates))]
        #mat = map2(log, mat)

        heatmap(mat, rlabels=map(str, rates), clabels=map(str, rates),
                display=True, xmargin=100, ymargin=100)


        

if __name__ == "__main__":
    unittest.main(testRunner=TestRunner())
