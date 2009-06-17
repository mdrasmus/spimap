#!/usr/bin/env python
# Mon Mar 30 15:08:52 EDT 2009
# make training data to test my training code

from rasmus.common import *

while not os.path.exists("bin/gene-tree-sim"):
    os.chdir("..")

sys.path.append("python")

import spidir


#=============================================================================
# flies.gamma.uniform.param
assert os.system("PYTHONPATH=~/projects/spidir/python:$PYTHONPATH "
                 "bin/gene-tree-sim "
                 "-s test/data/flies.norm.stree "
                 "-p test/data/flies.gamma.uniform.param "
                 "-n 100 "                 
                 "-l 1000 "
                 "--bgfreq .258,.267,.266,.209 "
                 "--tsvratio 1.59 "
                 "-D 0.0002 "
                 "-L 0.0001 "
                 "-O test/data/flies-one2one.uniform "
                 "-I .info "                 
                 ) == 0


assert os.system("phylofiles test/data/flies-one2one.uniform/ .tree | "
                 "bin/make-branch-matrix "
                 "-s test/data/flies.norm.stree "
                 "-S test/data/flies.smap > test/data/flies-one2one.uniform.lens"
                 ) == 0

assert os.system("phylofiles test/data/flies-one2one.uniform/ .tree | "
                 "bin/spidir-prep "
                 "-A .align "
                 "-T .tree "
                 "-s test/data/flies.norm.stree "
                 "-S test/data/flies.smap > test/data/flies-one2one.uniform.treemat.tmp"
                 ) == 0


def noise_matrix(gene_sizes, lens):
    return [[poissonvariate(x * d) / float(d) for x in row]
            for d, row in izip(gene_sizes, lens)]        
        
dat = read_delim("test/data/flies-one2one.uniform.treemat.tmp")
species = dat[0][2:]
lens = map2(float, submatrix(dat, range(1, len(dat)),
                             range(2, len(dat[0]))))
gene_sizes = map(int, cget(dat, 1)[1:])
lens2 = noise_matrix(gene_sizes, lens)

for i in xrange(len(lens2)):
    for j in xrange(len(lens2[0])):
        dat[i+1][j+2] = lens2[i][j]

write_delim("test/data/flies-one2one.uniform.treemat", dat)



#=============================================================================
# flies.gamma.fastleaf.param
assert os.system("PYTHONPATH=~/projects/spidir/python:$PYTHONPATH "
                 "bin/gene-tree-sim "
                 "-s test/data/flies.norm.stree "
                 "-p test/data/flies.gamma.fastleaf.param "
                 "-n 100 "                 
                 "-l 1000 "
                 "--bgfreq .258,.267,.266,.209 "
                 "--tsvratio 1.59 "
                 "-D 0.0002 "
                 "-L 0.0001 "
                 "-O test/data/flies-one2one.fastleaf "
                 "-I .info "                 
                 ) == 0


assert os.system("phylofiles test/data/flies-one2one.fastleaf/ .tree | "
                 "xargs cat | "
                 "bin/make-branch-matrix "
                 "-s test/data/flies.norm.stree "
                 "-S test/data/flies.smap > test/data/flies-one2one.fastleaf.lens"
                 ) == 0



#=============================================================================
# lets look at what we generated

if 0:
    stree = readTree("test/data/flies.norm.stree")
    dat = read_delim("test/data/flies-one2one.fastleaf.lens")
    nodes = dat[0]
    lens = map2(float, dat[1:])
    gene_rates = map(sum, lens)
    rlens = [[v / g for v in row]
             for row, g in zip(lens, gene_rates)]

    sp = 1
    snodes = list(stree.postorder())
    l = [i/snodes[sp].dist for i in cget(rlens, sp)]
    p = plotdistrib(l)
    p.plotfunc(lambda x: gammaPdf(x, (4, 4)), 0, max(l), .001)
    
