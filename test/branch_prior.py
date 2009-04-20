
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

if os.system("which xpdf 2>/dev/null") != 0:
    rplot_set_viewer("display")


class TestBranchPrior (unittest.TestCase):

    def test_branch_prior_approx(self):
        """Test branch prior"""

        prep_dir("test/output/branch_prior")
        #out = open("test/output/branch_prior/flies.nt.approx.txt", "w")
        out = sys.stderr

        treeids = os.listdir("test/data/flies.nt")[:1]
        treeids = ["0"]

        for treeid in treeids:
        
            tree = readTree("test/data/flies.nt/%s/%s.tree" % (treeid, treeid))

            drawTree(tree)
            
            stree = readTree("test/data/flies.norm.stree")
            gene2species = genomeutil.readGene2species("test/data/flies.smap")
            params = spidir.read_params("test/data/flies.nt.param")
            birth = .4
            death = .39
            nsamples = 1000
        
            recon = phylo.reconcile(tree, stree, gene2species)
            events = phylo.labelEvents(tree, recon)
            p = [spidir.branch_prior(tree, stree, recon, events,
                                     params, birth, death,
                                     nsamples, True)
                 for i in xrange(30)]
            p2 = [spidir.branch_prior(tree, stree, recon, events,
                                      params, birth, death,
                                      nsamples, False)
                 for i in xrange(30)]
            print >>out, "\t".join(map(str, [treeid, mean(p), sdev(p),
                                             mean(p2), sdev(p2)]))

        #out.close()

        
    def test_branch_prior_predup(self):
        """Test branch prior"""

        prep_dir("test/output/branch_prior_predup")
        #out = open("test/output/branch_prior/flies.nt.approx.txt", "w")
        out = sys.stderr
        treeid = "predup"

        tree = readTree("test/data/flies.predup.tree")
        drawTree(tree)
            
        stree = readTree("test/data/flies.norm.stree")
        gene2species = genomeutil.readGene2species("test/data/flies.smap")
        params = spidir.read_params("test/data/flies.nt.param")
        birth = .4
        death = .39
        nsamples = 1000
        
        recon = phylo.reconcile(tree, stree, gene2species)
        events = phylo.labelEvents(tree, recon)
        p = [spidir.branch_prior(tree, stree, recon, events,
                                 params, birth, death,
                                 nsamples, True)
             for i in xrange(30)]
        p2 = [spidir.branch_prior(tree, stree, recon, events,
                                  params, birth, death,
                                  nsamples, False)
              for i in xrange(30)]
        print >>out, "\t".join(map(str, [treeid, mean(p), sdev(p),
                                         mean(p2), sdev(p2)]))

        #out.close()
        

        
        
        
if __name__ == "__main__":
    unittest.main(testRunner=TestRunner())
