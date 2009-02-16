"""

    SPIDIR

    Test birth death functions

"""

import sys, unittest, ctypes


sys.path.append("python")
import spidir

from test import *

from rasmus.common import *
from rasmus.bio import phylo



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
                         for d in range(maxdoom+1))
            doomtable[i] = log(p)
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
    

def calcBirthDeathPrior(tree, stree, recon, birth, death, maxdoom):

    # TODO: must add implied speciation nodes for this to fully work

    events = phylo.labelEvents(tree, recon)
    phylo.addImpliedSpecNodes(tree, stree, recon, events)
    
    pstree, snodes, snodelookup = spidir.make_ptree(stree)

    treelib.drawTreeNames(tree, minlen=5)

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
                    
    

#=============================================================================

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
        doomtable = spidir.c_list(ctypes.c_float, doomtable)
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

        # TODO: simulate random trees and test their doom tables


    def test_birthDeathPrior(self):
        """test birth death prior (simple)"""
        
        l = 2
        u = .5
        maxdoom = 10

        stree = treelib.parseNewick("((A:1,B:1):1,((C:1,D:1):2,E:3):1);")
        #tree = treelib.parseNewick("((((a1,a2),(a3,a4)),(b1,b2)),((c1,d1),(c2,c3)));")
        tree = treelib.parseNewick("((((a1,a2),(a3,a4)),(b1,b2)),(((c1,d1),(c2,d2)),e1));")
        #tree = treelib.parseNewick("((a1,b1),(c1,d1));")
        def gene2species(gene):
            return gene[:1].upper()
        recon = phylo.reconcile(tree, stree, gene2species)
        events = phylo.labelEvents(tree, recon)

        ptree, nodes, nodelookup = spidir.make_ptree(tree)
        pstree, snodes, snodelookup = spidir.make_ptree(stree)

        ctree = spidir.tree2ctree(tree)
        cstree = spidir.tree2ctree(stree)
        recon2 = spidir.make_recon_array(tree, recon, nodes, snodelookup)
        events2 = spidir.make_events_array(nodes, events)

        print recon2
        print events2

        doomtable = spidir.c_list(ctypes.c_float, [0] * len(stree.nodes))

        spidir.calcDoomTable(cstree, l, u, maxdoom, doomtable)
        doomtable2 = calcDoomTable(stree, l, u, maxdoom)
        print list(doomtable)
        print list(doomtable2)
        
        p = spidir.birthDeathTreePrior(ctree, cstree,
                                       spidir.c_list(ctypes.c_int, recon2), 
                                       spidir.c_list(ctypes.c_int, events2),
                                       l, u, doomtable, maxdoom)
        p2 = calcBirthDeathPrior(tree, stree, recon, l, u, maxdoom)
        
        spidir.deleteTree(ctree)
        spidir.deleteTree(cstree)

        print "prior", p, p2
        fequal(p, p2)


    def test_birthDeathPriorFull(self):
        """test birth death prior with implied speciation nodes"""
        
        l = 2
        u = .5
        maxdoom = 10

        stree = treelib.parseNewick("((A:1,B:1):1,((C:1,D:1):2,E:3):1);")
        #tree = treelib.parseNewick("((((a1,a2),(a3,a4)),(b1,b2)),((c1,d1),(c2,c3)));")
        #tree = treelib.parseNewick("((((a1,a2),(a3,a4)),(b1,b2)),((c1,d1),(c2,d2)))")
        tree = treelib.parseNewick("((a1,b1),(c1,d1));")
        def gene2species(gene):
            return gene[:1].upper()
        recon = phylo.reconcile(tree, stree, gene2species)
        events = phylo.labelEvents(tree, recon)

        ptree, nodes, nodelookup = spidir.make_ptree(tree)
        pstree, snodes, snodelookup = spidir.make_ptree(stree)

        ctree = spidir.tree2ctree(tree)
        cstree = spidir.tree2ctree(stree)
        recon2 = spidir.make_recon_array(tree, recon, nodes, snodelookup)
        events2 = spidir.make_events_array(nodes, events)

        print recon2
        print events2

        doomtable = spidir.c_list(ctypes.c_float, [0] * len(stree.nodes))

        spidir.calcDoomTable(cstree, l, u, maxdoom, doomtable)
        doomtable2 = calcDoomTable(stree, l, u, maxdoom)
        print list(doomtable)
        print list(doomtable2)
        
        p = spidir.birthDeathTreePriorFull(ctree, cstree,
                                           spidir.c_list(ctypes.c_int, recon2), 
                                           spidir.c_list(ctypes.c_int, events2),
                                           l, u, doomtable, maxdoom)
        p2 = calcBirthDeathPrior(tree, stree, recon, l, u, maxdoom)
        
        spidir.deleteTree(ctree)
        spidir.deleteTree(cstree)

        print "prior", p, p2
        fequal(p, p2)


if __name__ == "__main__":
    unittest.main(testRunner=TestRunner())
