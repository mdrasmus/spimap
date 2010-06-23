
import sys, os

import pygsl
import pygsl.sf

while "python" not in os.listdir("."):
    os.chdir("..")

sys.path.append("python")
import spidir

from rasmus.common import *
from rasmus.bio import phylo, birthdeath
from test import *

try:
    from ctypes import cdll
    c = cdll.LoadLibrary("libc.so.6")
    c.srand(int(time.time()))
except:
    pass



if os.system("which xpdf 2>/dev/null") != 0:
    rplot_set_viewer("display")


def exc_default(func, val, exc=Exception):
    """Specify a default value for when an exception occurs"""
    try:
        return func()
    except exc:
        return val


class TestBranchPrior (unittest.TestCase):


    def test_branch_prior_simple1(self):
        """Test branch prior"""
        
        tree = treelib.parse_newick("((a1:1, b1:1):2, c1:3);")
        stree = treelib.parse_newick("((A:2, B:2):1, C:3);")
        
        gene2species = lambda x: x[0].upper()

        params = {"A": (1.0, 1.0),
                  "B": (3.0, 3.0),
                  "C": (4, 3.5),
                  2: (2.0, 2.0),
                  1: (1.0, 1.0),
                  "baserate": (11.0, 10.0)}
        birth = .01
        death = .02
        pretime = 1.0
        nsamples = 1
        
        recon = phylo.reconcile(tree, stree, gene2species)
        events = phylo.label_events(tree, recon)
        #pd(mapdict(recon, key=lambda x: x.name, val=lambda x: x.name))
        #pd(mapdict(events, key=lambda x: x.name))

        p = spidir.branch_prior(tree, stree, recon, events,
                                params, birth, death,
                                nsamples=nsamples, approx=False,
                                generate=1)


        tot = 0.0
        gs = list(frange(.0001, 4, .01))
        gs = list(frange(1, 1.01, .01))
        for g in gs:
            pg = invgammaPdf(g, params["baserate"])
            pa = gammaPdf(tree.nodes["a1"].dist,
                          [params["A"][0],
                           params["A"][1] / (g * stree.nodes["A"].dist)])

            pb = gammaPdf(tree.nodes["b1"].dist,
                          [params["B"][0],
                           params["B"][1] / (g * stree.nodes["B"].dist)])

            pc = spidir.gammaSumPdf(
                tree.nodes["c1"].dist + tree.nodes[2].dist, 2,
                [params["C"][0],
                 params[2][0]],
                [params["C"][1] / (g * stree.nodes["C"].dist),
                 params[2][1] / (g * stree.nodes[2].dist)], .001)

            print g, pg, pa, pb, pc
            tot += pg * pa * pb * pc
        tot /= len(gs)


        print (tree.nodes["c1"].dist + tree.nodes[2].dist,
               [params["C"][0], params[2][0]],
               [params["C"][1], params[2][1]])
        
        print "C", p
        print "P", log(tot)


    def test_branch_prior_simple2(self):
        """Test branch prior 2"""

        tree = treelib.parse_newick("((a1:2, a2:3):.4, b1:2);")        
        stree = treelib.parse_newick("(A:2, B:2);")
        
        gene2species = lambda x: x[0].upper()

        params = {"A": (1.0, 1.0),
                  "B": (3.0, 3.0),
                  1: (1.0, 1.0),
                  "baserate": (11.0, 10.0)}
        birth = .01
        death = .02
        pretime = 1.0
        nsamples = 100
        
        recon = phylo.reconcile(tree, stree, gene2species)
        events = phylo.label_events(tree, recon)
        #pd(mapdict(recon, key=lambda x: x.name, val=lambda x: x.name))
        #pd(mapdict(events, key=lambda x: x.name))

        p = spidir.branch_prior(tree, stree, recon, events,
                                params, birth, death,
                                nsamples=nsamples, approx=False)


        tot = 0.0

        gstart = 0.01
        gend = 3.0
        step = (gend - gstart) / 20.0;
        s2 = step / 2.0
        gs = list(frange(gstart+s2, gend+s2, step))
        for g in gs:
            pg = invgammaPdf(g, params["baserate"])

            pa = 0.0
            
            for i in range(nsamples):

                t = birthdeath.sample_birth_wait_time(1, stree.nodes["A"].dist,
                                                      birth, death)
                #print t
                
                t2 = stree.nodes["A"].dist - t

                pa1 = gammaPdf(tree.nodes["a1"].dist,
                              [params["A"][0],
                               params["A"][1] / (g * t2)])


                pa2 = gammaPdf(tree.nodes["a2"].dist,
                               [params["A"][0],
                                params["A"][1] / (g * t2)])            


                pb = spidir.gammaSumPdf(
                    tree.nodes["b1"].dist + tree.nodes[2].dist, 2,
                    [params["B"][0],
                     params["A"][0]],
                    [params["B"][1] / (g * stree.nodes["B"].dist),
                     params["A"][1] / (g * t)], .001)

                if "nan" not in map(str, [pa1, pa2, pb]):                    
                    pa += pa1 * pa2 * pb / nsamples
            
            tot += pg * pa * step
        #tot /= len(gs)

        print "unfold", (tree.nodes["b1"].dist + tree.nodes[2].dist,
               [params["B"][0],
                 params["A"][0]],
                [params["B"][1] / (g * stree.nodes["B"].dist),
                 params["A"][1] / (g * t)])

        print "C", p
        print "P", log(tot)


    def _test_branch_prior_samples(self):
        """Test branch prior"""

        prep_dir("test/output/branch_prior")

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
