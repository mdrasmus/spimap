#!/bin/sh python
#
# This is an example of how to use SPIMAP to reconstruct gene families.
#
# In this example, we reconstruct simulated data.
#

# Make sure tools are compiled and installed before running the commands in 
# this tutorial.  This can be done with this command:
cd ..
make install
cd examples

# or you can install spimap into the prefix of your choice
cd ..
make install prefix=YOUR_INSTALL_PREFIX
cd examples

# or you can run from the source directory by setting these environment 
# variables:
cd ..
make
export PYTHONPATH=$PYTHONOATH:../python
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:../lib
cd examples


#=============================================================================
# 0. Simulate new gene families with parameters similar to real fungi species

# get help information
../bin/spimap-sim -h

#Usage: spimap-sim [options]
#
#Options:
#  -h, --help            show this help message and exit
#  -s <species tree newick file>, --stree=<species tree newick file>
#  -p <params>, --params=<params>
#  -g <gene rate>, --generate=<gene rate>
#  --generatefunc=<gene rate function>
#  -l <gene length in base pairs>, --genelen=<gene length in base pairs>
#  --genelenprobs==<gene length distribution>
#  -r <transition/transversion ratio>, --tsvratio=<transition/transversion ratio>
#  -b <A>,<C>,<G>,<T>, --bgfreq=<A>,<C>,<G>,<T>
#  -m <min # of genes per tree>, --minsize=<min # of genes per tree>
#  -M <max # of genes per tree>, --maxsize=<max # of genes per tree>
#  -D <duplication rate>, --duprate=<duplication rate>
#  -L <loss rate>, --lossrate==<loss rate>
#  -n <number of trees to produce>, --ntrees=<number of trees to produce>
#  --start=<starting number>
#  --nospecrates         do not use species rates
#  -O <output directory>, --outtree=<output directory>
#  -T <output tree extension>, --outtreeext=<output tree extension>
#  -A <output align extension>, --outalignext=<output align extension>
#  -F <output sequence FASTA extension>, --outseqext=<output sequence FASTA extension>
#  -I <output information extenstion>, --outinfoext=<output information extenstion>
#  --nodir               do not create sub-directories
#  --resume              

# simulate 100 fungi families with genes that are 1000bp in length
../bin/spimap-sim \
    -s config/fungi.stree \
    -p config/sim/fungi.params \
    -l 1000 \
    -m 4 \
    -D 0.000564 \
    -L 0.003056 \
    -n 100 \
    -O fungi-sim


# simulate 100 one-to-one fungi families with genes that are 1000bp in length
../bin/spimap-sim \
    -s config/fungi.stree \
    -p config/sim/fungi.params \
    -l 1000 \
    -D 0 \
    -L 0 \
    -n 100 \
    -O fungi-sim-one2ones




#=============================================================================
# 1. Learn substitution rates from training set of 1-to-1 orthologous gene 
# families (fungi-sim-one2ones/*/*.nt.tree)
# These trees should all contain exactly one gene per species and be 
# congruent to the species tree "config/fungi.stree".

# 1.a. Extract all branch lengths from the newick trees 
# "fungi-sim-one2ones/*/*.tree" into a branch length matrix file 
# "train/fungi.lens" using the script "spimap-prep-rates"

# get help information
../bin/spimap-prep-rates -h

#Usage: spimap-prep-rates [options]
#
#Options:
#  -h, --help            show this help message and exit
#  -T <tree file extension>, --tree=<tree file extension>
#  -A <alignment file extension>, --align=<alignment file extension>
#  -s <species tree newick file>, --stree=<species tree newick file>
#  -S <gene to species file>, --smap=<gene to species file>
#  -l <branch length matrix output file>, --lenmatrix=<branch length matrix output file>
#  -r, --reroot 
#

# extract branch lengths
mkdir -p train
find fungi-sim-one2ones -name '?.tree' -or -name '??.tree' | \
../bin/spimap-prep-rates \
    -T .tree \
    -s config/fungi.stree \
    -S config/fungi.smap \
    -l train/fungi-sim.lens


# 1.b. Use branch length matrix file "train/fungi-sim.lens" to learn 
# substitution rate parameters "train/fungi-sim.param"

# get help information
../bin/spimap-train-rates -h

#Usage: spimap-train-rates [options]
#
#Options:
#  -h, --help            show this help message and exit
#  -s <species tree newick file>, --stree=<species tree newick file>
#  -l <branch length matrix output file>, --lenmatrix=<branch length matrix output file>
#  -n <number of iterations>, --niter=<number of iterations>
#  -r <number of rate categories>, --nrates=<number of rate categories>
#  -p <output parameters file>, --params=<output parameters file>
#  -f <treelen percentile>, --filter=<treelen percentile>
#  -m <minimum subst. count on a branch>, --mincount=<minimum subst. count on a branch>


# train rates parameters
../bin/spimap-train-rates \
    -s config/fungi.stree \
    -l train/fungi-sim.lens \
    -p train/fungi-sim.params



#=============================================================================
# 2. Learn genome-wide duplication/loss rates from gene counts in gene family
# clusters.

# 2.a. Use the script "spimap-prep-duploss" to extract genes per species counts
# from FASTA alignments "fungi-sim/*/*.align"

# get help information
../bin/spimap-prep-duploss -h

#Usage: spimap-prep-duploss [options]
#
#Options:
#  -h, --help            show this help message and exit
#  -s <species tree newick file>, --stree=<species tree newick file>
#  -S <gene to species file>, --smap=<gene to species file>
#  -c <gene count matrix output file>, --countmatrix=<gene count matrix output file>
#  -p <gene partition file>, --part=<gene partition file>

# get gene counts
find fungi-sim/*/*.align | ../bin/spimap-prep-duploss \
    -s config/fungi.stree \
    -S config/fungi.smap \
    -c train/fungi-sim.counts
    

# 2.b. Use the script "spimap-train-duploss" to learn duplication and loss
# rates from gene counts "train/fungi-sim.counts".

# get help information
../bin/spimap-train-duploss -h

#Usage: spimap-train-duploss [options]
#
#Options:
#  -h, --help            show this help message and exit
#  -s <species tree newick file>, --stree=<species tree newick file>
#  -S <gene to species file>, --smap=<gene to species file>
#  -c <gene count matrix output file>, --countmatrix=<gene count matrix output file>
#  --maxgene=<maximum number of genes in ancestor>
#  --maxiter=<maximum number of ML iterations>
#  --birth=<initial birth rate>
#  --death=<initial death rate>
#  --step=<initial step size>
#  -r <start>,<step>,<stop>, --range=<start>,<step>,<stop>
#  -H <heatmap file>, --heatmap=<heatmap file>


../bin/spimap-train-duploss \
    -s config/fungi.stree \
    -S config/fungi.smap \
    -c train/fungi-sim.counts > train/fungi-sim.duploss


#=============================================================================
# 3. Reconstruct gene trees using trained parameters

# show help information
spimap -h

#Usage: spimap [OPTION]
#  -a,--align  <alignment fasta>
#    sequence alignment in fasta format
#  -S,--smap  <species map>
#    gene to species map
#  -s,--stree  <species tree>
#    species tree file in newick format
#  -p,--param  <params file>
#    substitution rate parameters file
#  -o,--output  <output filename prefix>
#    prefix for all output filenames
#
#Sequence evolution model
#  -k,--kappa  <transition/transversion ratio>
#    used for HKY model (default=estimate)
#  -f,--bgfreq  <A freq>,<C ferq>,<G freq>,<T freq>
#    background frequencies (default: estimate)
#
#Dup/loss evolution model
#  -D,--duprate  <duplication rate>
#    rate of a gene duplication (default=0.1)
#  -L,--lossrate  <loss rate>
#    probability of loss (default=0.1)
#
#Search
#  -i,--niter  <# iterations>
#    number of iterations
#  --quickiter  <quick iterations>
#    number of subproposals (default=50)
#
#Information
#  -V,--verbose  <verbosity level>
#    verbosity level 0=quiet, 1=low, 2=medium, 3=high


# run SPIMAP for family 0
spimap \
    -a fungi-sim/0/0.align \
    -s config/fungi.stree \
    -S config/fungi.smap \
    -p train/fungi-sim.params \
    -o fungi-sim/0/0.spimap \
    -D 0.000597 \
    -L 0.003022 \
    -i 100 \
    --quickiter 1000 \
    -V 1 --log -



# example of how to run SPIMAP for many families
/bin/ls fungi-sim | (while read FAMID; do
  NAME=fungi-sim/$FAMID/$FAMID
  spimap \
    -a $NAME.align \
    -s config/fungi.stree \
    -S config/fungi.smap \
    -p train/fungi-sim.params \
    -o $NAME.spimap \
    -D 0.000597 \
    -L 0.003022 \
    -i 1000 \
    --quickiter 1000 \
    -V 1
done)


# compute percent correct topologies
find fungi-sim -name '*.spimap.tree' | python -c '
import spidir
from rasmus import treelib
from compbio import phylo
import sys

ncorrect = 0
total = 0
for line in sys.stdin:
  treefile = line.rstrip()
  spimap_tree = treelib.read_tree(treefile)
  sim_tree = treelib.read_tree(treefile.replace(".spimap", ""))

  ncorrect += int(phylo.hash_tree(spimap_tree) == phylo.hash_tree(sim_tree))
  total += 1

print "total", total
print "correct", ncorrect
print "percent", ncorrect / float(total)
'


#=============================================================================
# Clean up
#
# Remove all outputs from example analysis
#

rm -f train/fungi*
find fungi-fams/ | egrep '\.tree$|\.log$' | xargs rm

