
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
       
    def test_prog(self):
        """Test the main program"""

        prep_dir("test/output/prog")

        treeids = os.listdir("test/data/flies.nt")[:1]

        for treeid in treeids:

            alignfile = "test/data/flies.nt/%s/%s.align" % (treeid, treeid)
            tree_correct = "test/data/flies.nt/%s/%s.tree" % (treeid, treeid)

            cmd = (#"valgrind "
                   "bin/spidir "
                   "-a %s "
                   "-S test/data/flies.smap "
                   "-s test/data/flies.norm.stree "
                   "-p test/data/flies.nt.param "
                   "-o test/output/prog/%s "
                   "-k 1.59 "
                   "--duprate .4 "
                   "--lossrate .39 "
                   "--no_spr_nbr --quicksamples 1 "
                   "-i 50 "
                   "--quickiter 1000 "
                   "--correct %s "
                   "-V 1 --log - ") \
                   % (alignfile, treeid, tree_correct)

            print cmd
            self.assertEqual(os.system(cmd), 0)

    def _test_prog_many(self):
        """Test the main program on many inputs"""

        prep_dir("test/output/prog_many")

        treeids = os.listdir("test/data/flies.nt")

        for treeid in treeids:

            alignfile = "test/data/flies.nt/%s/%s.align" % (treeid, treeid)
            tree_correct = "test/data/flies.nt/%s/%s.tree" % (treeid, treeid)

            cmd = ("bin/spidir "
                   "-a %s "
                   "-S test/data/flies.smap "
                   "-s test/data/flies.norm.stree "
                   "-p test/data/flies.nt.param "
                   "-o test/output/prog_many/%s "
                   "-k 1.59 "
                   "--bgfreq .258,.267,.266,.209 "
                   "--duprate .4 "
                   "--lossrate .39 "
                   "-i 1 "
                   "--quickiter 10 "
                   "--correct %s "
                   "-V 2 --log - ") \
                   % (alignfile, treeid, tree_correct)

            print "treeid:", treeid
            self.assertEqual(os.system(cmd), 0)
            
              
        
if __name__ == "__main__":
    unittest.main(testRunner=TestRunner())
