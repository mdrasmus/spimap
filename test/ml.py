"""

   test ML code

"""


import sys, unittest, ctypes
from pprint import pprint
from math import *

sys.path.append("python")
import spidir

from test import *

from rasmus import util, treelib
from rasmus.bio import phylo, fasta

class TestHKY (unittest.TestCase):
    
    def test_ml(self):
        """Test ML code"""

        # params
        bgfreq = [.258,.267,.266,.209]
        kappa = 1.59

        # data
        tree = treelib.readTree("test/data/0.nt.tree")
        align = fasta.readFasta("test/data/1.nt.align")


        likes = []
        dists = []

        nodes = sorted(tree.nodes.values(), key=lambda x: x.dist)

        util.tic("find ML")
        for i in range(40):
            l = spidir.find_ml_branch_lengths_hky(
                tree,
                util.mget(align, tree.leafNames()),
                bgfreq, kappa,
                parsinit=(i == 0),
                maxiter=1)
            
            dists.append([n.dist for n in nodes])
            likes.append(l)
        util.toc()

        print likes

        prep_dir("test/output/ml/")

        # distances plot
        util.rplot_start("test/output/ml/ml_branches.pdf")
        util.rplot("plot", util.cget(dists, -1),
                   ylim=[0, max(dists[0])], t="l",
                   main="branch length convergence",
                   xlab="iterations",
                   ylab="branch lengths (sub/site)")
        for d in zip(* dists):
            util.rplot("lines", d)
        util.rplot_end(True)

        # likelihood plot
        util.rplot_start("test/output/ml/ml_likelihood.pdf")
        util.rplot("plot", likes, t="l",
                   xlab="iterations",
                   ylab="log likelihood",
                   main="likelihood convergence")
        util.rplot_end(True)


if __name__ == "__main__":
    unittest.main(testRunner=TestRunner())
