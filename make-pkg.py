#!/usr/bin/env python
# package spidir

files = ["bin",
         "include",
         "python",
         "src",
         "examples",

         "Makefile",
         "Makefile.dep",
         "setup.py",
         "README.txt"]

exclude = [".*\.o$", "src/old/.*", "bin/spimap"]

include = ["lib"]

#=============================================================================

import os, sys, re, shutil
from subprocess import call, Popen, PIPE

pkgdir = sys.argv[1]

exclude_expr = "|".join(exclude)

p = Popen(["find"] + files, stdout=PIPE)
all_files = [x.rstrip("\n") for x in p.stdout.readlines()]
all_files = [x for x in all_files
             if not re.match(exclude_expr, x)] + include

for f in all_files:
    dest = os.path.join(pkgdir, f)
    destdir = os.path.dirname(dest)
    print f, "-->", dest

    if os.path.isfile(f):
        # copy file
        if not os.path.exists(destdir):
            os.makedirs(destdir)
        shutil.copy(f, dest)
    else:
        # make dir
        if not os.path.exists(dest):
            os.makedirs(dest)

call(["tar", "zcvf", pkgdir + ".tar.gz", pkgdir])
