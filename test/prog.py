
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

        outdir = "test/output/prog"
        prep_dir(outdir)

        treeids = os.listdir("test/data/flies-duploss")[:1]

        for treeid in treeids:

            alignfile = "test/data/flies-duploss/%s/%s.align" % (treeid, treeid)
            tree_correct = "test/data/flies-duploss/%s/%s.tree" % (treeid, treeid)

            self.run_spidir(alignfile, tree_correct, outdir+"/"+treeid)

    def test_prog_boot(self):
        """Test the main program with bootstrapping"""

        outdir = "test/output/prog"
        prep_dir(outdir)

        treeids = os.listdir("test/data/flies-duploss")[:1]

        for treeid in treeids:

            alignfile = "test/data/flies-duploss/%s/%s.align" % (treeid, treeid)
            tree_correct = "test/data/flies-duploss/%s/%s.tree" % (treeid, treeid)

            self.run_spidir(alignfile, tree_correct, outdir+"/"+treeid,
                            boot=1, iter=1,
                            valgrind="valgrind --leak-check=full ")



    def test_prog_many(self):
        """Test the main program"""

        outdir = "test/output/prog"
        prep_dir(outdir)

        treeids = os.listdir("test/data/flies-duploss")

        for treeid in treeids:

            alignfile = "test/data/flies-duploss/%s/%s.align" % (treeid, treeid)
            tree_correct = "test/data/flies-duploss/%s/%s.tree" % (treeid, treeid)

            self.run_spidir(alignfile, tree_correct, outdir+"/"+treeid)


    def run_spidir(self, alignfile, tree_correct, outprefix, boot=1,
                   iter=50, valgrind=""):
        """Test the main program"""        
        
        cmd = (valgrind +
                   "bin/spidir "
                   "-a %s "
                   "-S test/data/flies.smap "
                   "-s test/data/flies.stree "
                   "-p test/data/flies.param "
                   "-o %s "
                   "--kappa 1.0 "
                   "--duprate .4 "
                   "--lossrate .39 "
                   "--quicksamples 1 "
                   "-b %d "
                   "-i %d "
                   "--quickiter 100 "
                   "--correct %s "
                   "-V 1 --log - ") \
                   % (alignfile, outprefix, boot, iter, tree_correct)

        print cmd
        self.assertEqual(os.system(cmd), 0)

              
        
if __name__ == "__main__":
    unittest.main(testRunner=TestRunner())