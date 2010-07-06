"""

    SPIDIR

    Test coal functions

"""

import sys, unittest, ctypes, os


sys.path.append("python")
import spidir

from spidir.topology_prior import calc_doom_table
calcDoomTable = calc_doom_table


from test import *

from rasmus.common import *
from compbio import phylo, birthdeath, coal




#=============================================================================


class TestCoalSim (unittest.TestCase):

    def setUp(self):
        pass



    def calc_fit(self, prefix, hist, probs, tabsize=100, pvalue=.01):
        """Calculate the goodness of fit between a histogram and an expected
           density distribution"""

        N = sum(hist.cget("count"))
        aerrs = []
        rerrs = []
        folds = []
        
        for i, row in enumerate(hist):
            aerrs.append(abs(probs[i] - row["percent"]))
            rerrs.append(abs(probs[i] - row["percent"]) / row["percent"])
            folds.append(probs[i] / row["percent"])

        hist.add_col("prior", data=probs)


        expected = [N * p for p in hist.cget("prior")]

        hist.add_col("expect", data=expected)
        hist.add_col("abs_error", data=aerrs)
        hist.add_col("rel_err", data=rerrs)
        hist.add_col("fold", data=folds)

        # display tables
        hist2 = hist.get()
        if isinstance(hist2[0]["item"], str):
            for row in hist2:
                row["item"] = row["item"][:60] # truncate colsize
        hist2[:tabsize].write_pretty()
        hist.write(prefix + "_compare.tab")

        # plot prior versus frequency
        x, y = hist.cget("percent", "prior")
        rplot_start(prefix + "_plot.pdf")
        rplot("plot", x, y, log="xy",
              xlim=[.0001, 1], ylim=[.0001, 1],
              xlab="simulated percent",
              ylab="computed prior")
        rp.lines([.0001, 1], [.0001, 1], col="red")
        rplot_end(False)
        
        # test goodness of fit
        mincount = 10
        chisq = sum((o - e)**2 / e
                    for o, e in zip(* hist.cget("count", "expect"))
                    if e > mincount)
        df = count(lambda e: e>mincount, hist.cget("expect")) - 1
        pval = 1.0 - rp.pchisq(chisq, df)
        print "chisq =", chisq, "df =", df
        print "pvalue(chisq) =", pval


        # assert chi square fits
        self.assert_(pval > pvalue)


    def do_test_coal_sim(self, stree, gene2species, n,
                         ntrees=10000, tabsize=30):
        """Perform a coal gene tree simulation test"""

        tops = []
        lookup = {}

        util.tic("simulating %d trees" % ntrees)
        for i in xrange(ntrees):
            tree, recon = coal.sample_multicoal_tree(stree, n,
                                                     namefunc=lambda x: x)
            tops.append(phylo.hash_tree(tree))
            lookup[tops[-1]] = (tree, recon)
        util.toc()
        
        hist = histtab(tops)

        probs = []
        for row in hist:
            tree, recon= lookup[row["item"]]
            try:
                #treelib.draw_tree_names(tree, maxlen=5)
                treelib.remove_single_children(tree)
                nodes = set(tree.postorder())
                for node, snode in recon.items():
                    if node not in nodes:
                        del recon[node]
                p = coal.prob_coal_recon_topology(tree, recon, stree, n)
            except:
                draw_tree(tree, maxlen=5, minlen=5)
                raise
            probs.append(exp(p))

        return hist, probs


    def test_coal_sim(self):
        """test coal prior against simulation"""

        def gene2species(gene):
            return gene[:1].upper()

        n = 1000
        stree = treelib.parse_newick("((A:1000,B:1000):1000,C:2000);")
        hist, probs = self.do_test_coal_sim(
            stree, gene2species, n=n, ntrees=1000)
        outdir = "test/output/coal_sim"
        prep_dir(outdir)
        self.calc_fit(outdir + "/sim_prior", hist, probs)


    def test_coal_sim_flies(self):
        """test coal prior against simulation"""

        def gene2species(gene):
            return gene[:1].upper()

        n = 1
        stree = treelib.read_tree("test/data/flies.stree")
        treelib.draw_tree_names(stree, maxlen=5)
        hist, probs = self.do_test_coal_sim(
            stree, gene2species, n=n, ntrees=10000)
        outdir = "test/output/coal_sim_flies"
        prep_dir(outdir)
        self.calc_fit(outdir + "/sim_prior_flies", hist, probs)



        

if __name__ == "__main__":
    unittest.main(testRunner=TestRunner())



