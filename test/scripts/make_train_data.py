#!/usr/bin/env python
# Mon Mar 30 15:08:52 EDT 2009
# make training data to test my training code

from rasmus.common import *

while not os.path.exists("bin/gene-tree-sim"):
    os.chdir("..")

sys.path.append("python")

import spidir


#=============================================================================
# flies.param one2one
assert os.system("bin/gene-tree-sim "
                 "-s test/data/flies.stree "
                 "-p test/data/flies.param "
                 "-n 100 "                 
                 "-l 1000 "
                 "--bgfreq .258,.267,.266,.209 "
                 "--tsvratio 1.59 "
                 "-D 0.000001 "
                 "-L 0.000002 "
                 "-O test/data/flies-one2one "
                 "-I .info "                 
                 ) == 0

assert os.system("phylofiles test/data/flies-one2one/ .tree | "
                 "bin/spidir-prep "
                 "-A .align "
                 "-T .tree "
                 "-l test/data/flies.treemat "
                 "-s test/data/flies.stree "
                 "-S test/data/flies.smap > test/data/flies-one2one.treemat.tmp"
                 ) == 0

#=============================================================================
# flies.param 
assert os.system("bin/gene-tree-sim "
                 "-s test/data/flies.stree "
                 "-p test/data/flies.param "
                 "-n 100 "                 
                 "-l 1000  "
                 "--bgfreq .258,.267,.266,.209 "
                 "--tsvratio 1.59 "
                 "-D 0.0024 "
                 "-L 0.0024 "
                 "-O test/data/flies-duploss "
                 "-I .info "                 
                 ) == 0

