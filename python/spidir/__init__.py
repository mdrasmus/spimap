#
# Python module for SPIDIR
#


#
# requires rasmus library
#

import os
from math import *
from ctypes import *

from rasmus import treelib, util
from rasmus.bio import fasta

import pyspidir
libdir = os.path.join(os.path.dirname(__file__), "..", "..", "lib")
spidir = cdll.LoadLibrary(os.path.join(libdir, "libspidir.so"))


def export(lib, funcname, return_type, arg_types, scope=globals()):
    """Exports a C function"""

    scope[funcname] = lib.__getattr__(funcname)
    scope[funcname].restype = return_type
    scope[funcname].argtypes = arg_types


export(spidir, "inumHistories", c_int, [c_int])
export(spidir, "numHistories", c_double, [c_int])
export(spidir, "numTopologyHistories", c_double, [c_void_p])
export(spidir, "birthDeathCount", c_float, [c_int, c_float, c_float, c_float])
export(spidir, "makeTree", c_void_p, [c_int, POINTER(c_int)])
export(spidir, "setTreeDists", c_void_p, [c_void_p, POINTER(c_float)])
export(spidir, "deleteTree", c_int, [c_void_p])
export(spidir, "calcDoomTable", c_int, [c_void_p, c_float, c_float,
                                           c_int, POINTER(c_float)])
export(spidir, "birthDeathTreePrior", c_float,
       [c_void_p, c_void_p, POINTER(c_int), POINTER(c_int),
        c_float, c_float, POINTER(c_float), c_int])
export(spidir, "birthDeathTreePriorFull", c_float,
       [c_void_p, c_void_p, POINTER(c_int), POINTER(c_int),
        c_float, c_float, POINTER(c_float), c_int])
export(spidir, "sampleDupTimes", c_int, [c_void_p, c_void_p,
                                         POINTER(c_int), POINTER(c_int),
                                         c_float, c_float])


def c_list(c_type, lst):
    """Make a C array from a list"""
    list_type = c_type * len(lst)
    return list_type(* lst)


#=============================================================================


def make_ptree(tree):
    """Make parent tree array from tree"""

    nodes = []
    nodelookup = {}
    ptree = []
    
    def walk(node):
        for child in node.children:
            walk(child)
        nodes.append(node)
    walk(tree.root)
    
    def leafsort(a, b):
        if a.isLeaf():
            if b.isLeaf():
                return 0
            else:
                return -1
        else:
            if b.isLeaf():
                return 1
            else:
                return 0
    
    # bring leaves to front
    nodes.sort(cmp=leafsort)
    nodelookup = util.list2lookup(nodes)
    
    for node in nodes:
        if node == tree.root:
            ptree.append(-1)
        else:
            ptree.append(nodelookup[node.parent])
    
    assert nodes[-1] == tree.root
    
    return ptree, nodes, nodelookup


def ptree2ctree(ptree):
    """Makes a c++ Tree from a parent array"""
    
    pint = c_int * len(ptree)
    tree = makeTree(len(ptree), pint(* ptree))
    return tree


def tree2ctree(tree):
    """Make a c++ Tree from a treelib.Tree datastructure"""

    ptree, nodes, nodelookup = make_ptree(tree)
    dists = [x.dist for x in nodes]
    ctree = ptree2ctree(ptree)
    setTreeDists(ctree, c_list(c_float, dists))
    return ctree


def make_gene2species_array(stree, nodes, snodelookup, gene2species):
    gene2speciesarray = []
    for node in nodes:
        if node.isLeaf():
            gene2speciesarray.append(snodelookup[
                                     stree.nodes[gene2species(node.name)]])
        else:
            gene2speciesarray.append(-1)
    return gene2speciesarray


def make_recon_array(tree, recon, nodes, snodelookup):
    recon2 = []
    for node in nodes:
        recon2.append(snodelookup[recon[node]])
    return recon2


def make_events_array(nodes, events):
    mapping = {"gene": 0,
               "spec": 1,
               "dup": 2}
    return util.mget(mapping, util.mget(events, nodes))


def parsimony(aln, tree):    
    ptree, nodes, nodelookup = make_ptree(tree)
    leaves = [x.name for x in nodes if isinstance(x.name, str)]
    seqs = util.mget(aln, leaves)
    
    dists = pyspidir.parsimony(ptree, seqs)
    
    for i in xrange(len(dists)):
        nodes[i].dist = dists[i]


def sample_gene_rate(tree, stree, gene2species, params,
                     aln,
                     bgfreq, tsvratio,
                     nsamples=100):

    ptree, nodes, nodelookup = make_ptree(tree)
    pstree, snodes, snodelookup = make_ptree(stree)
    smap = make_gene2species_array(stree, nodes, snodelookup, gene2species)

    mu = [float(params[snode.name][0]) for snode in snodes]
    sigma = [float(params[snode.name][1]) for snode in snodes]
    alpha, beta = params['baserate']

    seqs = [aln[node.name] for node in nodes if isinstance(node.name, str)]

    generates = []
    def callback(generate):
        generates.append(generate)
    
    pyspidir.sample_gene_rate(nsamples,
                              ptree, 
                              pstree,
                              smap,
                              mu,
                              sigma,
                              alpha, beta,
                              bgfreq,
                              tsvratio,
                              seqs,
                              callback)

    return generates



def est_gene_rate(tree, stree, gene2species, params,
                  aln, bgfreq, tsvratio,
                  nsamples=1000):

    generates = sample_gene_rate(tree, stree, gene2species, params,
                                 aln,
                                 bgfreq, tsvratio,
                                 nsamples)
    generates.sort()
    low = generates[int(.05 * len(generates))]
    high = generates[int(.95 * len(generates))]
    
    return sum(generates) / float(nsamples), (low, high)

