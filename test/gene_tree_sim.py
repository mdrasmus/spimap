"""

    Test simulation program

"""





import sys, unittest, ctypes, shutil
from pprint import pprint
from math import *

sys.path.append("python")
import spidir

from test import *

from rasmus.common import *
from rasmus import treelib
from rasmus.bio import phylo, birthdeath




class TestSim (unittest.TestCase):

    def setUp(self):
        pass

    def _test_sim(self):
        """simple test of simulation program"""

        prep_path("test/output/sim")

        ret = os.system(
            "bin/gene-tree-sim "
            "-p test/data/flies.nt.param "
            "-s test/data/flies.stree "
            "-D 3 -L 1 "
            "-n 10 "
            "-O test/output/sim")

        self.assertEqual(ret, 0)



    def _test_expected_births(self):
        """ensure proper number of births are occurring"""

        T = 1.5
        birth = 2.0
        death = 0.0
        ntrees = 1000

        birth_counts = []  # number of births along a random path

        for i in xrange(ntrees):
            tree = birthdeath.sample_birth_death_tree(T, birth, death)
            #treelib.drawTree(tree)

            node = tree.root
            i = 0
            while not node.isLeaf():
                node = node.children[random.randint(0, 1)]
                i += 1
            birth_counts.append(i)

        fequal(mean(birth_counts), birth * T, .05)


    def _test_survivor_distrib(self):
        """ensure proper number of survivors"""

        T = 1.2
        birth = 2.0
        death = 1.0
        ntrees = 10000

        survivors = []     # number of survivors in a tree

        for i in xrange(ntrees):
            tree = birthdeath.sample_birth_death_tree(T, birth, death)
            if len(tree.leaves()) == 1:
                if tree.leaves()[0].dist >= T:
                    survivors.append(1)
                else:
                    survivors.append(0)
            else:
                survivors.append(len(tree.leaves()))

        doomprob = birthdeath.prob_birth_death(1, 0, T, birth, death)
        
        h = [(1.0 - doomprob) * x/float(ntrees)
             for x in util.hist_int(survivors)]
        h[0] = doomprob
        h2 = [birthdeath.prob_birth_death(1, x, T, birth, death)
              for x in range(0, 10)]

        diff = util.vsub(h[:10], h2[:10])
        
        #print h[:10]
        #print h2[:10]
        self.assert_(max(map(abs, diff)) < .01)


    def test_birth_wait_time(self):
        """ensure reconstructed birth wait time has right distribution"""

        T = 1.2
        birth = 2.0
        death = 1.2

        times = [birthdeath.sample_birth_wait_time(1, T, birth, death)
                 for i in xrange(10000)]

        hx, hy = util.distrib(times, 30)

        cond = 1.0 - birthdeath.prob_no_birth(1, T, birth, death)
        x = list(frange(0, T, T/40.0))        
        y = [birthdeath.birth_wait_time(t, 1, T, birth, death) / cond
             for t in x]

        print sum(y) * T / 40.0

        prep_dir("test/output/sim-wait")
        rplot_start("test/output/sim-wait/wait_time.pdf")
        rplot("plot", hx, hy, t="l", ylim=[0, max(hy)])
        rp.lines(x, y, col="red")
        rplot_end(False)
                 


if __name__ == "__main__":
    unittest.main(testRunner=TestRunner())
