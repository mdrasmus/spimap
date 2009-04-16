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
assert os.system("bin/gene-tree-sim "
                 "-s test/data/flies.norm.stree "
                 "-p test/data/flies.nt.param "
                 "-n 1000 "
                 "-l 1000 "
                 "--bgfreq .258,.267,.266,.209 "
                 "--tsvratio 1.59 "
                 "-D 0.4 "
                 "-L 0.39 "
                 "-O test/data/flies.nt "
                 "-I .info "
                 ) == 0

