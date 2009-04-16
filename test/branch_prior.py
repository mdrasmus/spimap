
import sys, os

import pygsl
import pygsl.sf

while "python" not in os.listdir("."):
    os.chdir("..")

sys.path.append("python")
import spidir

from rasmus.common import *
from rasmus.bio import phylo
from test import *

if os.system("which xpdf") != 0:
    rplot_set_viewer("display")


class TestBranchPrior (unittest.TestCase):

    def test_branch_prior(self):
        """Test branch prior"""

        prep_dir("test/output/branch_prior")
        
        tree = readTree("test/data/sample.tree")
        tree2 = readTree("test/data/sample2.tree")
        stree = readTree("test/data/sample.stree")
        gene2species = genomeutil.readGene2species("test/data/sample.smap")
        params = spidir.read_params("test/data/sample.param")
        birth = .4
        death = .39
        
        
        recon = phylo.reconcile(tree, stree, gene2species)
        events = phylo.labelEvents(tree, recon)
        p = [spidir.branch_prior(tree, stree, recon, events,
                                 params, birth, death)
             for i in xrange(30)]
        print mean(p), sdev(p)
        
        recon2 = phylo.reconcile(tree2, stree, gene2species)
        events2 = phylo.labelEvents(tree2, recon2)
        p = [spidir.branch_prior(tree2, stree, recon2, events2,
                                 params, birth, death)
             for i in xrange(30)]
        
        print mean(p), sdev(p)
        
        
if __name__ == "__main__":
    unittest.main(testRunner=TestRunner())
