
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


class TestAllTerms (unittest.TestCase):


    def _test_all_terms(self):
        """Test all terms"""

        prep_dir("test/output/all_terms")
        out = open("test/output/all_terms/flies.nt.txt", "w")
        #out = sys.stderr

        treeids = os.listdir("test/data/flies.nt")[:100]
        #treeids = ["0"]

        for treeid in treeids:
        
            tree = readTree("test/data/flies.nt/%s/%s.tree" % (treeid, treeid))
            align = readFasta("test/data/flies.nt/%s/%s.align" % (treeid, treeid))

            print >>out, treeid
            drawTree(tree, out=out)
            
            stree = readTree("test/data/flies.norm.stree")
            gene2species = genomeutil.readGene2species("test/data/flies.smap")
            params = spidir.read_params("test/data/flies.nt.param")
            birth = .4
            death = .39
            pretime = 1.0
            nsamples = 100
            maxdoom = 20
            bgfreq = [.258,.267,.266,.209]
            kappa = 1.59
        
            recon = phylo.reconcile(tree, stree, gene2species)
            events = phylo.labelEvents(tree, recon)

            branchp, topp, seqlk = spidir.calc_joint_prob(
                align, tree, stree, recon, events, params,
                birth, death, pretime,
                bgfreq, kappa, maxdoom=maxdoom, terms=True)
            joint = topp + branchp + seqlk
            
            
            print >>out, "topp   ", topp
            print >>out, "branchp", branchp
            print >>out, "seqlk  ", seqlk
            print >>out, "joint  ", joint


        out.close()

        
    def test_search(self):
        """Test all terms"""

        prep_dir("test/output/all_terms_search")
        out = open("test/output/all_terms_search/flies.nt.txt", "w")
        #out = sys.stderr

        treeids = os.listdir("test/data/flies.nt")
        #treeids = ["3"]

        for treeid in treeids:
        
            tree_correct = readTree("test/data/flies.nt/%s/%s.tree" %
                                    (treeid, treeid))
            align = readFasta("test/data/flies.nt/%s/%s.align" %
                              (treeid, treeid))

            phylo.hashOrderTree(tree_correct)

            print >>out, treeid
            print >>out, "correct"
            drawTree(tree_correct, out=out)
            
            stree = readTree("test/data/flies.norm.stree")
            gene2species = genomeutil.readGene2species("test/data/flies.smap")
            params = spidir.read_params("test/data/flies.nt.param")
            birth = .4
            death = .39
            pretime = 1.0
            maxdoom = 20
            bgfreq = [.258,.267,.266,.209]
            kappa = 1.59

            genes = align.keys()
            seqs = align.values()
            
            tree = spidir.search_climb(genes, seqs,
                                       stree, gene2species,
                                       params, birth, death, pretime,
                                       bgfreq, kappa,
                                       maxdoom=maxdoom,
                                       niter=50, quickiter=100,
                                       nsamples=100, branch_approx=True)

            phylo.hashOrderTree(tree)
            
            

            print >>out, "constructed"
            drawTree(tree, out=out)
            

            print >>out, "is_correct:", (phylo.hashTree(tree) ==
                                         phylo.hashTree(tree_correct))
            

        out.close()
        
if __name__ == "__main__":
    unittest.main(testRunner=TestRunner())
