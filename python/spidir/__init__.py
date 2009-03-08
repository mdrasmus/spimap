#
# Python module for SPIDIR
#


#
# requires rasmus library
#

import os
from math import *
from ctypes import *

# rasmus imports
from rasmus import treelib, util
from rasmus.bio import fasta, phylo

#import spidir C lib
libdir = os.path.join(os.path.dirname(__file__), "..", "..", "lib")
spidir = cdll.LoadLibrary(os.path.join(libdir, "libspidir.so"))


#=============================================================================
# ctypes help functions

def c_list(c_type, lst):
    """Make a C array from a list"""
    list_type = c_type * len(lst)
    return list_type(* lst)

#def c_matrix(c_type, mat):
#    """Make a C matrix from a list of lists (mat)"""
#
#    row_type = c_type * len(mat[0])
#    mat_type = row_type * len(mat)
#    mat = mat_type(* [row_type(* row) for row in mat])
#    return cast(mat, POINTER(POINTER(c_type)))


def export(lib, funcname, return_type, arg_types, scope=globals()):
    """Exports a C function with documentation"""

    cfunc = lib.__getattr__(funcname)
    cfunc.restype = return_type
    cfunc.argtypes = arg_types[0::2]

    def func(*args):
        return cfunc(*args)

    scope[funcname] = func

    # set documentation
    args_doc = arg_types[1::2]
    scope[funcname].__doc__ = "%s(%s)" % (funcname, ",".join(args_doc))


# additional C types
c_float_p = POINTER(c_float)
c_int_p = POINTER(c_int)
c_char_p_p = POINTER(c_char_p)

#=============================================================================
# wrap functions from c library

# basic tree functions
export(spidir, "deleteTree", c_int, [c_void_p, "tree"])
export(spidir, "makeTree", c_void_p, [c_int, "nnodes",
                                          c_int_p, "ptree"])
export(spidir, "setTreeDists", c_void_p, [c_void_p, "tree",
                                              c_float_p, "dists"])


# topology prior birthdeath functions
export(spidir, "inumHistories", c_int, [c_int, "nleaves"])
export(spidir, "numHistories", c_double, [c_int, "nleaves"])
export(spidir, "numTopologyHistories", c_double, [c_void_p, "tree"])
export(spidir, "birthDeathCount", c_float,
       [c_int, "ngenes", c_float, "time",
        c_float, "birth", c_float, "death"])
export(spidir, "calcDoomTable", c_int,
       [c_void_p, "tree", c_float, "birth", c_float, "death",
        c_int, "maxdoom", c_float_p, "doomtable"])
export(spidir, "birthDeathTreePrior", c_float,
       [c_void_p, "tree", c_void_p, "stree",
        c_int_p, "recon", c_int_p, "events",
        c_float, "birth", c_float, "death",
        c_float_p, "doomtable", c_int, "maxdoom"])
export(spidir, "birthDeathTreePriorFull", c_float,
       [c_void_p, "tree", c_void_p, "stree",
        c_int_p, "recon", c_int_p, "events",
        c_float, "birth", c_float, "death",
        c_float_p, "doomtable", c_int, "maxdoom"])
export(spidir, "sampleDupTimes", c_int,
       [c_void_p, "tree", c_void_p, "stree",
        c_int_p, "recon", c_int_p, "events",
        c_float, "birth", c_float, "death"])

# branch prior functions
export(spidir, "treelk", c_float,
       [c_int, "nnodes", c_int_p, "ptree", c_float_p, "dists",
        c_int, "nsnodes", c_int_p, "pstree", 
        c_int_p, "recon", c_int_p, "events",
        c_float_p, "mu", c_float_p, "sigma",
        c_float, "generate", 
        c_float, "predupprob", c_float, "dupprob", c_float, "lossprob",
        c_float, "alpha", c_float, "beta", c_int, "onlyduploss"])


# parsimony
export(spidir, "parsimony", c_void_p,
       [c_int, "nnodes", c_int, "ptree", c_int, "nseqs",
        c_char_p_p, "seqs", c_float, "dists",
        c_int, "buildAncestral", c_char_p_p, "ancetralSeqs"])


# sequence likelihood functions
export(spidir, "makeHkyMatrix", c_void_p,
       [c_float_p, "bgfreq", c_float, "kappa", c_float, "time",
        c_float_p, "matrix"])
export(spidir, "makeHkyDerivMatrix", c_void_p,
       [c_float_p, "bgfreq", c_float, "kappa", c_float, "time",
        c_float_p, "matrix"])
export(spidir, "makeHkyDeriv2Matrix", c_void_p,
       [c_float_p, "bgfreq", c_float, "kappa", c_float, "time",
        c_float_p, "matrix"])
export(spidir, "branchLikelihoodHky", c_float,
       [c_float_p, "probs1", c_float_p, "probs2",
        c_int, "seqlen", c_float_p, "bgfreq", c_float, "kappa",
        c_float, "time"])
export(spidir, "branchLikelihoodHkyDeriv", c_float,
       [c_float_p, "probs1", c_float_p, "probs2",
        c_int, "seqlen", c_float_p, "bgfreq",
        c_float, "kappa", c_float, "time"])       
export(spidir, "branchLikelihoodHkyDeriv2", c_float,
       [c_float_p, "probs1", c_float_p, "probs2",
        c_int, "seqlen", c_float_p, "bgfreq",
        c_float, "kappa", c_float, "time"])
export(spidir, "mleDistanceHky", c_float,
       [c_float_p, "probs1", c_float_p, "probs2",
        c_int, "seqlen", c_float_p, "bgfreq",
        c_float, "kappa", c_float, "t0", c_float, "t1"])
export(spidir, "calcSeqProbHky", c_float,
       [c_void_p, "tree", c_int, "nseqs", c_char_p_p, "seqs",
        c_float_p, "bgfreq", c_float, "kappa"])
export(spidir, "findMLBranchLengthsHky", c_float,
       [c_int, "nnodes", c_int_p, "ptree", c_int, "nseqs",
        c_char_p_p, "seqs", c_float_p, "dists",
        c_float_p, "bgfreq", c_float, "kappa",
        c_int, "maxiter", c_int, "parsinit"])

#=============================================================================
# additional python interface


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



def calc_birth_death_prior(tree, stree, recon, birth, death, maxdoom,
                              events=None):

    assert birth != death, "birth and death must be different rates"

    if events is None:
        events = phylo.labelEvents(tree, recon)

    ptree, nodes, nodelookup = make_ptree(tree)
    pstree, snodes, snodelookup = make_ptree(stree)

    ctree = tree2ctree(tree)
    cstree = tree2ctree(stree)
    recon2 = make_recon_array(tree, recon, nodes, snodelookup)
    events2 = make_events_array(nodes, events)

    doomtable = c_list(c_float, [0] * len(stree.nodes))
    calcDoomTable(cstree, birth, death, maxdoom, doomtable)
    
    p = birthDeathTreePriorFull(ctree, cstree,
                                c_list(c_int, recon2), 
                                c_list(c_int, events2),
                                birth, death, doomtable, maxdoom)
    deleteTree(ctree)
    deleteTree(cstree)

    return p

def make_hky_matrix(bgfreq, kappa, t):
    """
    Returns a HKY matrix

    bgfreq -- the background frequency A,C,G,T
    kappa  -- is transition/transversion ratio
    """
    
    matrix = [0.0] * 16
    matrix = c_list(c_float, matrix)
    makeHkyMatrix(c_list(c_float, bgfreq), kappa, t, matrix)
    return [matrix[0:4],
            matrix[4:8],
            matrix[8:12],
            matrix[12:16]]

def make_hky_deriv_matrix(bgfreq, kappa, t):
    """
    Returns a HKY Derivative matrix

    bgfreq -- the background frequency A,C,G,T
    kappa  -- is transition/transversion ratio
    """
    
    matrix = [0.0] * 16
    matrix = c_list(c_float, matrix)
    makeHkyDerivMatrix(c_list(c_float, bgfreq), kappa, t, matrix)
    return [matrix[0:4],
            matrix[4:8],
            matrix[8:12],
            matrix[12:16]]

def make_hky_deriv2_matrix(bgfreq, kappa, t):
    """
    Returns a HKY 2nd Derivative matrix

    bgfreq -- the background frequency A,C,G,T
    kappa  -- is transition/transversion ratio
    """
    
    matrix = [0.0] * 16
    matrix = c_list(c_float, matrix)
    makeHkyDeriv2Matrix(c_list(c_float, bgfreq), kappa, t, matrix)
    return [matrix[0:4],
            matrix[4:8],
            matrix[8:12],
            matrix[12:16]]


def branch_likelihood_hky(probs1, probs2, seqlen, bgfreq, kappa, time):
    return branchLikelihoodHky(c_list(c_float, probs1), c_list(c_float, probs2),
                               seqlen, c_list(c_float, bgfreq), kappa, time)

def branch_likelihood_hky_deriv(probs1, probs2, seqlen, bgfreq, kappa, time):
    return branchLikelihoodHkyDeriv(
        c_list(c_float, probs1), c_list(c_float, probs2), seqlen, 
        c_list(c_float, bgfreq), kappa, time)

def branch_likelihood_hky_deriv2(probs1, probs2, seqlen, bgfreq, kappa, time):
    return branchLikelihoodHkyDeriv2(
        c_list(c_float, probs1), c_list(c_float, probs2), seqlen, 
        c_list(c_float, bgfreq), kappa, time)

def mle_distance_hky(probs1, probs2, seqlen, bgfreq, kappa, t0, t1):    
    return mleDistanceHky(c_list(c_float, probs1), c_list(c_float, probs2),
                          seqlen, c_list(c_float, bgfreq), kappa,
                          t0, t1)


def calc_seq_likelihood_hky(tree, align, bgfreq, kappa):

    ctree = tree2ctree(tree)
    calign = (c_char_p * len(align))(* align)
    
    l = calcSeqProbHky(ctree, len(align), calign, c_list(c_float, bgfreq),
                       kappa)

    deleteTree(ctree)

    return l


def find_ml_branch_lengths_hky(tree, align, bgfreq, kappa, maxiter=20,
                               parsinit=True):

    ptree, nodes, nodelookup = make_ptree(tree)
    calign = (c_char_p * len(align))(* align)
    dists = c_list(c_float, [n.dist for n in nodes])

    l = findMLBranchLengthsHky(len(ptree), c_list(c_int, ptree),
                               len(align), calign,
                               dists, c_list(c_float, bgfreq), kappa,
                               maxiter, int(parsinit))
    
    for i, node in enumerate(nodes):
        node.dist = dists[i]
    
    return l



#=============================================================================
# pyspidir (older code)
# TODO: convert to ctypes


'''
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


'''
