
import sys, os

import pygsl
import pygsl.sf

while "python" not in os.listdir("."):
    os.chdir("..")

sys.path.append("python")
import spidir

from rasmus.common import *
from rasmus.bio import phylo

from test import *

if os.system("which xpdf 2>/dev/null") != 0:
    rplot_set_viewer("display")


def exc_default(func, val, exc=Exception):
    """Specify a default value for when an exception occurs"""
    try:
        return func()
    except exc:
        return val


class TestProg (unittest.TestCase):
       
    def test_vert(self):
        """Test the main program"""

        outdir = "test/output/prog_large_trees"
        prep_dir(outdir)
        
        cmd = (#"valgrind "
                   "bin/spidir "
                   "-a test/data/verts/19520/19520.nt.mfa "
                   "-S test/data/verts/ensembl.smap "
                   "-s test/data/verts/ensembl.stree "
                   "-p test/data/verts/ensembl.param "
                   "-o %s/a "
                   "-k 1.59 "
                   "--duprate .4 "
                   "--lossrate .39 "
                   "--quicksamples 1 "
                   "-i 50 "
                   "--quickiter 1000 "
                   "-V 2 --log - " % outdir)

        print cmd
        self.assertEqual(os.system(cmd), 0)


    def test_fungi(self):
        """Test the main program"""

        outdir = "test/output/prog_large_trees"
        prep_dir(outdir)
        
        cmd = (#"valgrind "
                   "bin/spidir "
                   "-a test/data/fungi/10515/10515.align "
                   "-S test/data/fungi.smap "
                   "-s test/data/fungi.stree "
                   "-p test/data/fungi.param "
                   "-o %s/a "
                   "-k 3.483 "
                   "--duprate 0.000732 "
                   "--lossrate 0.000859 "
                   "--quicksamples 1 "
                   "-i 50 "
                   "--quickiter 1 "
                   "-V 2 --log - " % outdir)

        print cmd
        self.assertEqual(os.system(cmd), 0)
        
if __name__ == "__main__":
    unittest.main(testRunner=TestRunner())
