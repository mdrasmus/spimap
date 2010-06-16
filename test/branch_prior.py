
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


def exc_default(func, val, exc=Exception):
    """Specify a default value for when an exception occurs"""
    try:
        return func()
    except exc:
        return val


class TestBranchPrior (unittest.TestCase):

    def test_branch_prior_samples(self):
        """Test branch prior"""

        prep_dir("test/output/branch_prior")
        out = open("test/output/branch_prior/flies.approx.txt", "w")
        out = sys.stderr

        treeids = os.listdir("test/data/flies")
        treeids = ["3"]

        for treeid in treeids:
        
            tree = read_tree("test/data/flies-duploss/%s/%s.tree" % (treeid, treeid))

            print treeid
            draw_tree(tree)
            
            stree = read_tree("test/data/flies.stree")
            gene2species = phylo.read_gene2species("test/data/flies.smap")
            params = spidir.read_params("test/data/flies.param")
            birth = .0012
            death = .0013
            pretime = 1.0
            nsamples = 100
        
            recon = phylo.reconcile(tree, stree, gene2species)
            events = phylo.label_events(tree, recon)

            p = [spidir.branch_prior(tree, stree, recon, events,
                                     params, birth, death,
                                     nsamples=nsamples, approx=True)
                 for i in xrange(30)]

            #row = [treeid,
            #       mean(p), exc_default(lambda: sdev(p), INF)]
            print treeid, p

            #print >>out, "\t".join(map(str, row))
            #self.assert_(INF not in row and -INF not in row)

        out.close()


    def _test_branch_prior_approx(self):
        """Test branch prior"""

        prep_dir("test/output/branch_prior")
        out = open("test/output/branch_prior/flies.approx.txt", "w")
        out = sys.stderr

        treeids = os.listdir("test/data/flies")

        for treeid in treeids:
        
            tree = read_tree("test/data/flies-duploss/%s/%s.nt.tree" % (treeid, treeid))

            print treeid
            draw_tree(tree)
            
            stree = read_tree("test/data/flies.stree")
            gene2species = phylo.read_gene2species("test/data/flies.smap")
            params = spidir.read_params("test/data/flies.param")
            birth = .0012
            death = .0013
            pretime = 1.0
            nsamples = 100
        
            recon = phylo.reconcile(tree, stree, gene2species)
            events = phylo.label_events(tree, recon)
            p = [spidir.branch_prior(tree, stree, recon, events,
                                     params, birth, death,
                                     nsamples=nsamples, approx=False)
                 for i in xrange(30)]
            p2 = [spidir.branch_prior(tree, stree, recon, events,
                                     params, birth, death,
                                     nsamples=nsamples, approx=True)
                 for i in xrange(30)]


            row = [treeid,
                   mean(p), exc_default(lambda: sdev(p), INF),
                   mean(p2),exc_default(lambda: sdev(p2), INF)]

            print >>out, "\t".join(map(str, row))
            self.assert_(INF not in row and -INF not in row)

        out.close()

        
    def _test_branch_prior_predup(self):
        """Test branch prior"""

        prep_dir("test/output/branch_prior_predup")
        #out = open("test/output/branch_prior/flies.nt.approx.txt", "w")
        out = sys.stderr
        treeid = "predup"

        tree = read_tree("test/data/flies.predup.tree")
        drawTree(tree)
            
        stree = read_tree("test/data/flies.stree")
        gene2species = phylo.read_gene2species("test/data/flies.smap")
        params = spidir.read_params("test/data/flies.param")
        birth = .4
        death = .39
        pretime = 1.0
        nsamples = 100
        
        recon = phylo.reconcile(tree, stree, gene2species)
        events = phylo.label_events(tree, recon)
        p = [spidir.branch_prior(tree, stree, recon, events,
                                 params, birth, death, pretime,
                                 nsamples, True)
             for i in xrange(30)]
        p2 = [spidir.branch_prior(tree, stree, recon, events,
                                  params, birth, death, pretime,
                                  nsamples, False)
              for i in xrange(30)]
        print >>out, "\t".join(map(str, [treeid, mean(p), sdev(p),
                                         mean(p2), sdev(p2)]))

        #out.close()
        

        
        
        
if __name__ == "__main__":
    unittest.main(testRunner=TestRunner())
