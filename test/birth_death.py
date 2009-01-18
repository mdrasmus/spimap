"""

    SPIDIR

    Test birth death functions

"""

import sys, unittest


sys.path.append("python")
import spidir


from rasmus.common import *



def birthDeathCount(ngenes, time, birthRate, deathRate):
    """
    Returns the probability that one lineage leaves 'ngenes' genes
    after time 'time'
    """
    
    l = birthRate
    u = deathRate
    r = l - u
    a = u / l

    ut = (1.0 - exp(-r*time)) / (1.0 - a * exp(-r*time))
    p0 = a*ut
    
    if ngenes == 0:
        return p0
    
    # (1.0 - p0)*(1.0 - ut) * ut^{ngenes-1}
    return (1.0 - p0)*(1.0 - ut) * (ut**(ngenes-1))



class TestBirthDeath (unittest.TestCase):

    def setUp(self):
        pass


    def test_birthDeathCount(self):
        l = 3
        u = .5

        for t in frange(0, 5, .01):
            for s in range(0, 20):
                p1 = spidir.birthDeathCount(s, t, l, u)
                p2 = birthDeathCount(s, t, l, u)
                #print t, s, p1, p2
                self.assert_(abs(p1 - p2) < .0001)


    def test_numHistories(self):
        for i in range(1, 20):
            n = spidir.numHistories(i)
            n2 = int(factorial(i) * factorial(i-1) / 2**(i-1))
            self.assert_(abs(n - n2) / abs(n + n2) < .00001)
        

    def test_make_trees(self):
        tree = treelib.parseNewick("((a,b),(c,d));")
        ctree = spidir.tree2ctree(tree)
        spidir.deleteTree(ctree)

        ctree = spidir.ptree2ctree([3,3,4,4,-1])
        spidir.deleteTree(ctree)

        ctree = spidir.makeTree(5, spidir.c_list(spidir.c_int, [3,3,4,4,-1]))
        spidir.deleteTree(ctree)


    def test_numTopologyHistories(self):

        tree = spidir.tree2ctree(treelib.parseNewick("((a,b),(c,d));"))
        self.assertEqual(spidir.numTopologyHistories(tree), 2)
        spidir.deleteTree(tree)

        
        tree = spidir.tree2ctree(treelib.parseNewick("(((a,b),(c,d)),(e,f));"))
        self.assertEqual(spidir.numTopologyHistories(tree), 8)
        spidir.deleteTree(tree)

        tree = spidir.tree2ctree(treelib.parseNewick("((((a,b),(c,d)),(e,f)),g);"))
        self.assertEqual(spidir.numTopologyHistories(tree), 8)
        spidir.deleteTree(tree)

        tree = spidir.tree2ctree(treelib.parseNewick("((((a,b),(c,d)),(e,f)),((g,h),i));"))
        self.assertEqual(spidir.numTopologyHistories(tree), 7*6/2*8)
        spidir.deleteTree(tree)



if __name__ == "__main__":
    unittest.main()
