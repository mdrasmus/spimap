"""

    SPIDIR

    Test birth death functions

"""

import sys, unittest, ctypes


sys.path.append("python")
import spidir
from spidir import calc_birth_death_prior as c_calcBirthDeathPrior

from test import *

from rasmus.common import *
from rasmus import stats
from rasmus.bio import phylo, birthdeath


if os.system("which xpdf 2> /dev/null") != 0:
    rplot_set_viewer("display")


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


def calcDoomTable(tree, birthRate, deathRate, maxdoom):
    l = birthRate
    u = deathRate
    
    ptree, nodes, nodelookup = spidir.make_ptree(tree)

    doomtable = [0] * len(tree.nodes)
    
    def walk(node):
        if node.isLeaf():
            doomtable[nodelookup[node]] = -INF
        else:
            for child in node.children:
                walk(child)
                        
            i = nodelookup[node]
            p = 1.0
            for child in node:
                p *= sum(birthDeathCount(d, child.dist, l, u) *
                         exp(doomtable[nodelookup[child]]) ** d
                         for d in range(0, maxdoom+1))
            doomtable[i] = safelog(p, e, -INF)
    walk(tree.root)
    
    return doomtable


def numTopologyHistories(node, leaves=None):

    if leaves is None:
        leaves = node.leaves()
    leaves = set(leaves)

    prod = [1.0]

    def walk(node):        
        if node in leaves:
            return 0
        else:
            internals = map(walk, node.children)
            prod[0] *= choose(sum(internals), internals[0])
            return 1 + sum(internals)
    walk(node)

    return prod[0]


def numRedunantTopology(node, gene2species, leaves=None):

    if leaves is None:
        leaves = node.leaves()
    leaves = set(leaves)
    colors = {}
    nmirrors = [0]

    def walk(node):
        if node in leaves:
            colors[node] = phylo.hashTree(node, gene2species)
        else:
            # recurse
            for child in node.children:
                walk(child)
            
            childHashes = mget(colors, node.children)
            if len(childHashes) > 1 and util.equal(* childHashes):
                nmirrors[0] += 1
            
            childHashes.sort()
            colors[node] = phylo.hashTreeCompose(childHashes)
    walk(node)

    colorsizes = util.hist_dict(util.mget(colors, leaves)).values()
    
    val = 1
    for s in colorsizes:
        if s > 1:
            val *= stats.factorial(s)
    #print "py val=", val, "nmirrors=", nmirrors[0]
    return val / (2**nmirrors[0])


def numRedunantTopology2(node, gene2species, leaves=None):

    if leaves is None:
        leaves = node.leaves()
    colors = []

    for node in leaves:
        colors.append(phylo.hashTree(node, gene2species))
        
    colorsizes = util.hist_dict(colors).values()
    return stats.multinomial(colorsizes)


def getSubTree(node, snode, recon, events):
    leaves = []

    def walk(node):
        if recon[node] != snode:
            return
        if events[node] != "dup":
            leaves.append(node)
        else:
            for child in node.children:
                walk(child)
    walk(node)

    return leaves
    

def calcBirthDeathPrior_old(tree, stree, recon, birth, death, maxdoom):

    events = phylo.labelEvents(tree, recon)
    phylo.addImpliedSpecNodes(tree, stree, recon, events)
    
    pstree, snodes, snodelookup = spidir.make_ptree(stree)

    # get doomtable
    doomtable = calcDoomTable(stree, birth, death, maxdoom)

    prod = 0.0
    for node in tree:
        if events[node] == "spec":
            for schild in recon[node].children:
                nodes2 = [x for x in node.children if recon[x] == schild]
                if len(nodes2) > 0:
                    node2 = nodes2[0]
                    subleaves = getSubTree(node2, schild, recon, events)
                    nhist = numTopologyHistories(node2, subleaves)
                    s = len(subleaves)
                    thist = factorial(s) * factorial(s-1) / 2**(s-1)
                else:
                    nhist = 1.0
                    thist = 1.0
                    s = 0

                t = sum(birthDeathCount(s + i, schild.dist, birth, death) *
                        exp(doomtable[snodelookup[schild]]) ** i
                        for i in range(maxdoom+1))
                
                prod += log(nhist) - log(thist) + log(t)


    #phylo.removeImpliedSpecNodes(tree, recon, events)
    treelib.removeSingleChildren(tree)

    return prod



def calcBirthDeathPrior(tree, stree, recon, birth, death, maxdoom,
                        events=None):

    if events is None:
        events = phylo.labelEvents(tree, recon)
    phylo.addImpliedSpecNodes(tree, stree, recon, events)
    
    
    pstree, snodes, snodelookup = spidir.make_ptree(stree)
    leaves = set(tree.leaves())

    def gene2species(gene):
        return recon[tree.nodes[gene]].name

    # get doomtable
    doomtable = calcDoomTable(stree, birth, death, maxdoom)

    prod = 0.0
    for node in tree:
        if events[node] == "spec":
            for schild in recon[node].children:
                nodes2 = [x for x in node.children if recon[x] == schild]
                if len(nodes2) > 0:
                    node2 = nodes2[0]
                    subleaves = getSubTree(node2, schild, recon, events)
                    nhist = numTopologyHistories(node2, subleaves)  * \
                            numRedunantTopology2(node2, gene2species, subleaves)
                    
                    if len(set(subleaves) & leaves) > 0:
                        nhist *= numRedunantTopology(node2, gene2species, subleaves)
                    s = len(subleaves)
                    thist = factorial(s) * factorial(s-1) / 2**(s-1)
                else:
                    nhist = 1.0
                    thist = 1.0
                    s = 0

                t = sum(birthDeathCount(s + i, schild.dist, birth, death) *
                        exp(doomtable[snodelookup[schild]]) ** i
                        for i in range(maxdoom+1))

                prod += log(nhist) - log(thist) + safelog(t)
                #print "py prod2", t, s

    # correct for renumbering
    nt = numRedunantTopology(tree.root, gene2species)
    #print "py x", nt, prod
    prod -= log(nt)

    
    treelib.removeSingleChildren(tree)

    return prod



#=============================================================================


class TestBirthDeath2 (unittest.TestCase):

    def setUp(self):
        pass


    def test_sample_birth(self):

        prep_dir("test/output/birth_death_sample")

        n = 1
        T = 1.0
        birth = 2.0
        death = 0.5

        tic("sampling")
        samples = [birthdeath.sample_birth_wait_time(n, T, birth, death)
                   for i in xrange(10000)]        
        toc()

        tic("samplingC")
        samples2 = [spidir.sampleBirthWaitTime1(T, birth, death)
                   for i in xrange(10000)]        
        toc()

        cond = 1.0 - birthdeath.prob_no_birth(n, T, birth, death)

        x, y2 = distrib(samples, 20)
        x2, y3 = distrib(samples2, 20)
        y = [birthdeath.birth_wait_time(i, n, T, birth, death) / cond
              for i in x]

        rplot_start("test/output/birth_death_sample/sample.pdf")
        rplot("plot", x, y, t="l")
        rp.lines(x, y2, col="red")
        rp.lines(x2, y3, col="blue")
        rplot_end(True)
              
        


class TestBirthDeath (unittest.TestCase):

    def setUp(self):
        pass


    def test_birthDeathCount(self):
        """birthDeathCount"""
        l = 3
        u = .5

        for t in frange(0, 5, .01):
            for s in range(0, 20):
                p1 = spidir.birthDeathCount(s, t, l, u)
                p2 = birthDeathCount(s, t, l, u)
                #print t, s, p1, p2
                fequal(p1, p2, .01)

        #s = 0
        #1 = 2
        #u = .5
        #p1 = spidir.birthDeathCount(s, t, l, u)
        #p2 = birthDeathCount(s, t, l, u)
        #fequal(sp1, p2)


    def test_numHistories(self):
        """numHistories"""
        for i in range(1, 20):
            n = spidir.numHistories(i)
            n2 = int(factorial(i) * factorial(i-1) / 2**(i-1))
            self.assert_(abs(n - n2) / abs(n + n2) < .00001)
        

    def test_make_trees(self):
        """makeTree"""
        tree = treelib.parseNewick("((a,b),(c,d));")
        ctree = spidir.tree2ctree(tree)
        spidir.deleteTree(ctree)

        ctree = spidir.ptree2ctree([3,3,4,4,-1])
        spidir.deleteTree(ctree)

        ctree = spidir.makeTree(5, spidir.c_list(spidir.c_int, [3,3,4,4,-1]))
        spidir.deleteTree(ctree)


    def test_numTopologyHistories(self):
        """numTopologyHistories"""
        
        tree = treelib.parseNewick("((a,b),(c,d));")
        ctree = spidir.tree2ctree(tree)
        self.assertEqual(spidir.numTopologyHistories(ctree), 2)
        self.assertEqual(numTopologyHistories(tree.root), 2)
        spidir.deleteTree(ctree)

        tree = treelib.parseNewick("(((a,b),(c,d)),(e,f));")
        ctree = spidir.tree2ctree(tree)
        self.assertEqual(spidir.numTopologyHistories(ctree), 8)
        self.assertEqual(numTopologyHistories(tree.root), 8)
        spidir.deleteTree(ctree)

        tree = treelib.parseNewick("((((a,b),(c,d)),(e,f)),g);")
        ctree = spidir.tree2ctree(tree)
        self.assertEqual(spidir.numTopologyHistories(ctree), 8)
        self.assertEqual(numTopologyHistories(tree.root), 8)
        spidir.deleteTree(ctree)

        tree = treelib.parseNewick("((((a,b),(c,d)),(e,f)),((g,h),i));")
        ctree = spidir.tree2ctree(tree)
        self.assertEqual(spidir.numTopologyHistories(ctree), 7*6/2*8)
        self.assertEqual(numTopologyHistories(tree.root), 7*6/2*8)
        spidir.deleteTree(ctree)


    def assert_doom_table(self, tree, l, u, maxdoom):
        ctree = spidir.tree2ctree(tree)
        ptree, nodes, nodelookup = spidir.make_ptree(tree)        

        treelib.drawTree(tree, scale=10)

        doomtable = [0] * len(tree.nodes)
        doomtable = spidir.c_list(ctypes.c_double, doomtable)
        spidir.calcDoomTable(ctree, l, u, maxdoom, doomtable)
        spidir.deleteTree(ctree)
        
        doomtable2 = calcDoomTable(tree, l, u, maxdoom)
        
        print list(doomtable)        
        print doomtable2

        for i, j in zip(doomtable, doomtable2):            
            fequal(i, j)
        

    def test_calcDoomTable(self):
        """test doom table calculation"""
        l = 2
        u = .5
        maxdoom = 10

        tree = treelib.parseNewick("(((a:1,b:1):1,(c:1,d:1):2):1,(e:1,f:1):1);")
        self.assert_doom_table(tree, l, u, maxdoom)

        tree = treelib.parseNewick("(((a:1,b:2):.3,c:5):1,(e:1,f:1):1);")
        self.assert_doom_table(tree, l, u, maxdoom)

        tree = treelib.parseNewick("((A:1,B:1):1,(C:1.5,D:1.5):0.5);")
        self.assert_doom_table(tree, l, u, maxdoom)

        # TODO: simulate random trees and test their doom tables


    def test_birthDeathPrior(self):
        """test birth death prior (simple)"""
        
        l = 2
        u = .5
        maxdoom = 10

        def gene2species(gene):
            return gene[:1].upper()


        
        stree = treelib.parseNewick("((A:1,B:1):1,((C:1,D:1):2,E:3):1);")
        
        tree = treelib.parseNewick("((((a1,a2),(a3,a4)),(b1,b2)),(((c1,d1),(c2,d2)),e1));")
        recon = phylo.reconcile(tree, stree, gene2species)
        p = c_calcBirthDeathPrior(tree, stree, recon, l, u, maxdoom)
        p2 = calcBirthDeathPrior(tree, stree, recon, l, u, maxdoom)
        print "prior", p, p2
        fequal(p, p2)

        # test gene reconciling within species tree
        tree = treelib.parseNewick("((((a1,a2),(a3,a4)),(b1,b2)),((c1,d1),(c2,c3)));")
        recon = phylo.reconcile(tree, stree, gene2species)
        p = c_calcBirthDeathPrior(tree, stree, recon, l, u, maxdoom)
        p2 = calcBirthDeathPrior(tree, stree, recon, l, u, maxdoom)
        print "prior", p, p2
        fequal(p, p2)

        # test gene reconciling within species tree
        tree = treelib.parseNewick("((a1,b1),c1);")
        recon = phylo.reconcile(tree, stree, gene2species)
        p = c_calcBirthDeathPrior(tree, stree, recon, l, u, maxdoom)
        p2 = calcBirthDeathPrior(tree, stree, recon, l, u, maxdoom)
        print "prior", p, p2
        fequal(p, p2)


        # test case that occurred during simulation
        # non parsimonious reconciliation
        stree = treelib.parseNewick("((A:1,B:1):1,C:2);")
        tree = treelib.parseNewick("((a1,a2));")
        recon = {tree.nodes["a1"]: stree.nodes["A"],
                 tree.nodes["a2"]: stree.nodes["A"],
                 tree.nodes["a1"].parent: stree.nodes["A"].parent,
                 tree.root: stree.root}
        events = {tree.nodes["a1"]: "gene",
                 tree.nodes["a2"]: "gene",
                 tree.nodes["a1"].parent: "dup",
                 tree.root: "spec"}
        p = c_calcBirthDeathPrior(tree, stree, recon, l, u, maxdoom,
                                  events=events)
        p2 = calcBirthDeathPrior(tree, stree, recon, l, u, maxdoom,
                                 events=events)
        tree.writeNewick(oneline=True)
        print "\nprior", p, p2
        fequal(p, p2)

        # complicated case
        stree = treelib.parseNewick("((A:1,B:1):1,C:2);")
        tree = treelib.parseNewick("((((B2:1.072961,B8:1.072961):0.106756,((((A1:0.427377,(((A3:0.150067,A11:0.150067):0.038521,A2:0.188588):0.121082,A5:0.309671):0.117706):0.352590,A9:0.779967):0.113269,(A8:0.266488,A7:0.266488):0.626747):0.236597,(((B9:0.160640,B7:0.160640):0.098506,B4:0.259146):0.429865,B5:0.689011):0.440822):0.049885):0.714463,(B13:1.086980,((A10:1.000000,((B10:0.408524,(((B3:0.143778,B1:0.143778):0.023788,B6:0.167566):0.058639,B12:0.226204):0.182319):0.232105,B11:0.640629):0.359371):0.082149,(A6:0.277757,A4:0.277757):0.804392):0.004830):0.807201):0.105819,(C3:1.213803,(((C6:0.190132,C4:0.190132):0.011461,C5:0.201593):0.745740,(C1:0.017299,C2:0.017299):0.930034):0.266470):0.786197);")
        recon = phylo.reconcile(tree, stree, gene2species)
        p = c_calcBirthDeathPrior(tree, stree, recon, l, u, maxdoom)
        p2 = calcBirthDeathPrior(tree, stree, recon, l, u, maxdoom)
        print "prior", p, p2
        fequal(p, p2)

        stree = treelib.parseNewick(
                "(((A:1,B:1):1,(C:1.5,D:1.5):0.5):.5,((E:.2,F:.2):.6):1.9);")
        tree = treelib.parseNewick("(((A1:1.000000,B1:1.000000):1.000000,(((C2:0.718949,C1:0.718949):0.168784,C3:0.887733):0.612267,D1:1.500000):0.500000):0.500000,((F8:0.122975,F5:0.122975):6.518970,(((E4:0.200000,F6:0.200000):5.257236,((E3:0.200000,F7:0.200000):4.029009,(E2:0.200000,F1:0.200000):4.029009):1.228227):0.306982,(((E5:0.200000,F3:0.200000):1.068443,(E6:0.200000,F2:0.200000):1.068443):1.094596,(E1:0.200000,F4:0.200000):2.163039):3.401179):0.877727):1.458055);")
        
        recon = phylo.reconcile(tree, stree, gene2species)
        p = c_calcBirthDeathPrior(tree, stree, recon, l, u, maxdoom)
        p2 = calcBirthDeathPrior(tree, stree, recon, l, u, maxdoom)
        print "prior", p, p2
        fequal(p, p2)


        # test for overflow
        stree = treelib.parseNewick("((A:1,B:1):1,C:2);")
        tree = treelib.parseNewick("((((C24:0.940136,C6:0.940136):0.140529,((((C37:0.374306,(C26:0.054540,C10:0.054540):0.319766):0.046428,(C15:0.009875,C29:0.009875):0.410860):0.112550,(C3:0.213709,C28:0.213709):0.319576):0.034152,C13:0.567437):0.513228):0.545124,((((C36:0.036428,C30:0.036428):1.402769,(((C33:0.038848,C19:0.038848):0.352795,(C9:0.282410,(C1:0.000411,C21:0.000411):0.281998):0.109233):0.452052,((C34:0.108366,C12:0.108366):0.332454,C35:0.440820):0.402875):0.595502):0.039525,((((((C40:0.082790,C23:0.082790):0.003327,(C11:0.021474,C14:0.021474):0.064643):0.031631,C31:0.117748):0.019433,C17:0.137181):0.619636,C39:0.756818):0.139581,(C4:0.160113,(C41:0.116482,C32:0.116482):0.043631):0.736286):0.582323):0.000255,(C5:0.389128,((C25:0.112569,C27:0.112569):0.127253,(C22:0.139232,C18:0.139232):0.100590):0.149306):1.089849):0.146811):0.299534,(C2:1.197153,(C7:0.690311,(C16:0.070431,((C20:0.000466,C8:0.000466):0.060700,C38:0.061165):0.009265):0.619881):0.506842):0.728170);")
        print "leaves", len(tree.leaves())
        recon = phylo.reconcile(tree, stree, gene2species)
        p = c_calcBirthDeathPrior(tree, stree, recon, l, u, maxdoom)
        p2 = calcBirthDeathPrior(tree, stree, recon, l, u, maxdoom)
        print "prior", p, p2
        fequal(p, p2)

        
        

    def test_birthDeathPriorFull(self):
        """test birth death prior with implied speciation nodes"""
        
        l = 2
        u = .5
        maxdoom = 10

        def gene2species(gene):
            return gene[:1].upper()

        stree = treelib.parseNewick("((A:1,B:1):1,((C:1,D:1):2,E:3):1);")
        tree = treelib.parseNewick("((((a1,a2),(a3,a4)),(b1,b2)),((c1,d1),(c2,c3)));")
        
        # test gene reconciling within species tree
        recon = phylo.reconcile(tree, stree, gene2species)
        p = c_calcBirthDeathPrior(tree, stree, recon, l, u, maxdoom)
        p2 = calcBirthDeathPrior(tree, stree, recon, l, u, maxdoom)
        print "prior", p, p2
        fequal(p, p2)



    def test_numRedundantTopologies(self):
        """test numRedundantTopologies"""

        def gene2species(gene):
            return gene[:1].upper()
        
        #stree = treelib.parseNewick("((A:1,B:1):1,((C:1,D:1):2,E:3):1);")
        
        tree = treelib.parseNewick("((((a1,b1),c1),((a2,b2),c2)),d1);")        
        self.assertEqual(numRedunantTopology(tree.root, gene2species), 4)

        tree = treelib.parseNewick("((a1,b1));")
        self.assertEqual(numRedunantTopology(tree.root, gene2species), 1)

        tree = treelib.parseNewick("((a1,a2),c1);")
        self.assertEqual(numRedunantTopology(tree.root, gene2species), 1)

        tree = treelib.parseNewick("((((a1,a3),b2),(a2,b1)),c1)")
        self.assertEqual(numRedunantTopology(tree.root, gene2species), 6)



    def test_birthDeathPrior_large(self):
        """test birth death prior for large trees"""
        
        l = 0.000732 
        u = 0.000859
        maxdoom = 20
        
        stree = treelib.read_tree("test/data/fungi.stree")
        gene2species = phylo.read_gene2species("test/data/fungi.smap")
        tree = treelib.read_tree("test/data/fungi/10169/10169.tree")
        recon = phylo.reconcile(tree, stree, gene2species)

        p = c_calcBirthDeathPrior(tree, stree, recon, l, u, maxdoom)
        print p
        self.assert_(p != -INF)

if __name__ == "__main__":
    unittest.main(testRunner=TestRunner())
