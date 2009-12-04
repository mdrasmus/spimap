"""

    SPIDIR

    Test birth death functions

"""

import sys, unittest, ctypes, os


sys.path.append("python")
import spidir

from test import *
from test.birthdeath import \
     calcDoomTable, calcBirthDeathPrior, c_calcBirthDeathPrior, \
     numRedunantTopology, birthDeathCount, numTopologyHistories


from rasmus.common import *
from rasmus.bio import phylo, birthdeath


# NOTE: reactivate tests by removing leading "_"


#=============================================================================


def rename_leaves(tree, stree, gene2species):

    spnames = dict.fromkeys(stree.leafNames(), 1)
    leaves = tree.leaves()
    random.shuffle(leaves)

    for node in leaves:
        sp = gene2species(node.name)
        tree.rename(node.name, sp + "%d" % spnames[sp])
        spnames[sp] += 1



class TestBirthDeathSim (unittest.TestCase):

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
        hist[:tabsize].write_pretty()
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



    def do_test_doomed(self, stree, gene2species,
                                     duprate, lossrate,
                                     ntrees=10000, maxdoom=20):
        """Perform a doom table test"""

        doomtable = calcDoomTable(stree, duprate, lossrate, maxdoom)

        doomed = 0

        util.tic("simulating %d trees" % ntrees)
        for i in xrange(ntrees):
            tree, recon, events = birthdeath.sample_birth_death_gene_tree(
                stree, duprate, lossrate, removeloss=True)            
            if len(tree.nodes) == 0:
                doomed += 1
        util.toc()

        return doomed / float(ntrees), exp(doomtable[-1])


    def do_test_birth_death_gene_sim(self, stree, gene2species,
                                     duprate, lossrate,
                                     ntrees=10000, maxdoom=20, tabsize=30):
        """Perform a birth death gene tree simulation test"""

        doomtable = calcDoomTable(stree, duprate, lossrate, maxdoom)
        
        tops = []
        lookup = {}


        def rename_tree(tree, gene2species):
            if len(tree.nodes) == 0:
                return
            spcounts = util.hist_dict(map(gene2species, tree.leaf_names()))
            names = {}
            for sp, c in spcounts.items():
                names[sp] = range(1, c+1)
                random.shuffle(names[sp])

            for node in tree.leaves():
                sp = gene2species(node.name)
                tree.rename(node.name, sp + "." + str(names[sp].pop()))
            

        util.tic("simulating %d trees" % ntrees)
        for i in xrange(ntrees):
            tree, recon, events = birthdeath.sample_birth_death_gene_tree(
                stree, duprate, lossrate, 
                removeloss=True)
            rename_tree(tree, gene2species)
            phylo.add_implied_spec_nodes(tree, stree, recon, events)
            
            if len(tree.nodes) == 0:
                tops.append("()")
                lookup["()"] = (None, None, None)
            else:
                rename_leaves(tree, stree, gene2species)
                tops.append(phylo.hash_tree(tree))
                #tops.append(phylo.hash_tree(tree, gene2species))
                lookup[tops[-1]] = (tree, recon, events)
        util.toc()
        
        hist = histtab(tops)

        probs = []
        for row in hist:
            tree, recon, events = lookup[row["item"]]

            if tree is None:
                probs.append(exp(doomtable[-1]))
            else:                
                p = c_calcBirthDeathPrior(tree, stree, recon,
                                          duprate, lossrate,
                                          maxdoom, events=events)
                p2 = calcBirthDeathPrior(tree, stree, recon,
                                         duprate, lossrate,
                                         maxdoom, events=events)
                #try:
                #    fequal(p, p2)
                #except:
                # 
                #    recon2 = phylo.reconcile(tree, stree, gene2species)
                #    if recon != recon2:
                #        print "recon is not parsimonious"
                #    tree.write_newick(oneline=True)
                #    raise
                probs.append(exp(p2))

        return hist, probs




    def test_birth_death_single_sim(self):
        """test the single branch prior"""
        
        duprate = 2
        lossrate = .5
        ntrees = 1000
        tabsize = 100
        maxdoom = 20
        T = 1.0

        tops = []
        survivors = []
        lookup = {}

        stree = treelib.parse_newick("(A:1);")
        def gene2species(gene):
            return gene[:1].upper()

        util.tic("simulating %d trees" % ntrees)
        for i in xrange(ntrees):
            tree = birthdeath.sample_birth_death_tree(T, duprate, lossrate)
            rename_leaves(tree, stree, lambda x: "A")
            tops.append(phylo.hash_tree(tree, gene2species))
            survivors.append(len(tree.leaves()))
            lookup[tops[-1]] = tree
        util.toc()

        outdir = "test/output/birthdeath_sim_simple"
        prep_dir(outdir)
        
        hist_tops = histtab(tops)
        hist = histtab(survivors)

        no_extinct = 1.0 - birthDeathCount(0, T, duprate, lossrate)

        probs = []
        for row in hist:
            ngenes = row["item"]
            probs.append(birthDeathCount(ngenes, T, duprate, lossrate) /
                         no_extinct)

        probs_tops = []
        for row in hist_tops:
            tree = lookup[row["item"]]
            nhist = numTopologyHistories(tree.root) * \
                    numRedunantTopology(tree.root, gene2species, tree.leaves())
            s = len(tree.leaves())
            thist = factorial(s) * factorial(s-1) / 2**(s-1)
            
            #probs_tops.append(nhist / float(thist) *
            #                  birthDeathCount(s, T, duprate, lossrate) /
            #                  no_extinct)

            s = stree.leaves()[0]
            recon = dict((n, s) for n in tree)
            p = calcBirthDeathPrior(tree, stree, recon,
                                    duprate, lossrate,
                                    maxdoom)
            probs_tops.append(exp(p))

        self.calc_fit(outdir + "/sim_prior_ngenes", hist, probs)
        self.calc_fit(outdir + "/sim_prior_top", hist_tops, probs_tops)




    def test_birth_death_single2_sim(self):
        """test the single branch prior"""

        duprate = 2.0
        lossrate = .5
        T = 1.0
        maxdoom = 20

        
        stree = treelib.parse_newick("(A:1,B:1);")
        def gene2species(gene):
            return gene[:1].upper()
        s = stree.leaves()[0]

        b = birthDeathCount(1, T, duprate, lossrate)

        # 1
        tree = treelib.parseNewick("(a,b);")
        recon = phylo.reconcile(tree, stree, gene2species)
        p = birthDeathCount(1, T, duprate, lossrate) * b
        p2 = exp(calcBirthDeathPrior(tree, stree, recon, duprate, lossrate,
                                     maxdoom))
        p2 *= numRedunantTopology(tree.root, gene2species)
        print p, p2
        fequal(p, p2)
        
        # 2
        tree = treelib.parseNewick("((a,a),b);")
        recon = phylo.reconcile(tree, stree, gene2species)
        p = birthDeathCount(2, T, duprate, lossrate) * b
        p2 = exp(calcBirthDeathPrior(tree, stree, recon, duprate, lossrate,
                                     maxdoom))
        p2 *= numRedunantTopology(tree.root, gene2species)
        print p, p2
        fequal(p, p2)

        # 3
        tree = treelib.parseNewick("(((a,a),a),b);")
        recon = phylo.reconcile(tree, stree, gene2species)
        p = birthDeathCount(3, T, duprate, lossrate) * b
        p2 = exp(calcBirthDeathPrior(tree, stree, recon, duprate, lossrate,
                                     maxdoom))
        p2 *= numRedunantTopology(tree.root, gene2species)
        print p, p2
        fequal(p, p2)

        # 4
        tree = treelib.parseNewick("(((a,a),(a,a)),b);")
        recon = phylo.reconcile(tree, stree, gene2species)
        p = birthDeathCount(4, T, duprate, lossrate) * b / 3.0
        p2 = exp(calcBirthDeathPrior(tree, stree, recon, duprate, lossrate,
                                     maxdoom))
        p2 *= numRedunantTopology(tree.root, gene2species)
        print p, p2
        fequal(p, p2)


    def test_birth_death_single3_sim(self):
        """test the single branch prior"""

        duprate = 2.0
        lossrate = .5
        T = 1.0
        maxdoom = 20

        
        stree = treelib.parse_newick("(A:1,B:1);")
        def gene2species(gene):
            return gene[:1].upper()
        s = stree.leaves()[0]

        b = birthDeathCount(1, T, duprate, lossrate)

        # 1
        tree = treelib.parseNewick("(a,b);")
        recon = phylo.reconcile(tree, stree, gene2species)
        p = birthDeathCount(1, T, duprate, lossrate) * b
        p2 = exp(calcBirthDeathPrior(tree, stree, recon, duprate, lossrate,
                                     maxdoom))
        print p, p2
        fequal(p, p2)
        
        # 2
        tree = treelib.parseNewick("((a1,a2),b);")
        recon = phylo.reconcile(tree, stree, gene2species)
        p = birthDeathCount(2, T, duprate, lossrate) * b
        p2 = exp(calcBirthDeathPrior(tree, stree, recon, duprate, lossrate,
                                     maxdoom))        
        print p, p2
        fequal(p, p2)

        # 3
        tree = treelib.parseNewick("(((a1,a2),a3),b);")
        recon = phylo.reconcile(tree, stree, gene2species)
        p = birthDeathCount(3, T, duprate, lossrate) * b / 3.0
        p2 = exp(calcBirthDeathPrior(tree, stree, recon, duprate, lossrate,
                                     maxdoom))
        print p, p2
        fequal(p, p2)

        # 4
        tree = treelib.parseNewick("(((a1,a2),(a3,a4)),b);")
        recon = phylo.reconcile(tree, stree, gene2species)
        p = birthDeathCount(4, T, duprate, lossrate) * b / 3.0 / 3.0
        p2 = exp(calcBirthDeathPrior(tree, stree, recon, duprate, lossrate,
                                     maxdoom))
        print p, p2
        fequal(p, p2)


    def test_birth_death_sim(self):
        """test birth death prior against simulation"""

        def gene2species(gene):
            return gene[:1].upper()

        stree = treelib.parse_newick("((A:1,B:1):1,C:2);")
        hist, probs = self.do_test_birth_death_gene_sim(
            stree, gene2species,
            duprate=.2, lossrate=.1, ntrees=10000)
        outdir = "test/output/birthdeath_sim"
        prep_dir(outdir)
        self.calc_fit(outdir + "/sim_prior", hist, probs)
        
    def test_birth_death_sim2(self):
        """test birth death prior against simulation"""

        def gene2species(gene):
            return gene[:1].upper()
        
        stree = treelib.parse_newick("((A:1,B:1):1,(C:1.5,D:1.5):0.5);")
        hist, probs = self.do_test_birth_death_gene_sim(
            stree, gene2species,
            duprate=.5, lossrate=.4,
            ntrees=10000)
        outdir = "test/output/birthdeath_sim2"
        prep_dir(outdir)
        self.calc_fit(outdir + "/sim_prior", hist, probs)
        

    def test_birth_death_sim3(self):
        """test birth death prior against simulation"""

        def gene2species(gene):
            return gene[:1].upper()

        stree = treelib.parse_newick("((A:1,B:1):1,C:2);")
        hist, probs = self.do_test_birth_death_gene_sim(
            stree, gene2species,
            duprate=2, lossrate=1.5,
            ntrees=1000)
        outdir = "test/output/birthdeath_sim3"
        prep_dir(outdir)
        self.calc_fit(outdir + "/sim_prior", hist, probs)

    def test_birth_death_sim4(self):
        """test birth death prior against simulation"""

        def gene2species(gene):
            return gene[:1].upper()

        stree = treelib.parse_newick(
            "(((A:1,B:1):1,(C:1.5,D:1.5):0.5):.5,((E:.2,F:.2):6):1.9);")
        hist, probs = self.do_test_birth_death_gene_sim(
            stree, gene2species,
            duprate=.2, lossrate=.1,
            ntrees=10000)
        outdir = "test/output/birthdeath_sim4"
        prep_dir(outdir)
        self.calc_fit(outdir + "/sim_prior", hist, probs)
        

    def test_doom_table(self):

        def gene2species(gene):
            return gene[:1].upper()
        
        stree = treelib.parseNewick("((A:1,B:1):1,(C:1.5,D:1.5):0.5);")

        sim_doom, prior_doom = self.do_test_doomed(
            stree, gene2species,
            duprate=2, lossrate=1.5,
            ntrees=10000, maxdoom=10)
        print sim_doom, prior_doom
        
        

if __name__ == "__main__":
    unittest.main(testRunner=TestRunner())



