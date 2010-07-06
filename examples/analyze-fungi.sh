#!/bin/sh python
#
# This is an example of how to use SPIMAP to reconstruct gene families.
#
# In this example, we reconstruct the families of 16 fungi species
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
# 1. Learn substitution rates from training set of 1-to-1 orthologous gene 
# families (fungi-one2ones/*/*.nt.tree
# These trees should all contain exactly one gene per species and be 
# congruent to the species tree "config/fungi.stree".

# 1.a. Extract all branch lengths from the newick trees 
# "fungi-one2ones/*/*.mt.tree" into a branch length matrix file 
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
find fungi-one2ones -name '*.nt.tree' | \
../bin/spimap-prep-rates \
    -T .nt.tree \
    -A .nt.align \
    -s config/fungi.stree \
    -S config/fungi.smap \
    -l train/fungi.lens


# 1.b. Use branch length matrix file "train/fungi.lens" to learn substitution 
# rate parameters "train/fungi.param"

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
    -l train/fungi.lens \
    -p train/fungi.params



#=============================================================================
# 2. Learn genome-wide duplication/loss rates from gene counts in gene family
# clusters.

# 2.a. Use the script "spimap-prep-duploss" to extract genes per species counts
# from FASTA alignments "fungi-fams/*/*.nt.align"

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
find fungi-fams/*/*.nt.align | ../bin/spimap-prep-duploss \
    -s config/fungi.stree \
    -S config/fungi.smap \
    -c train/fungi.counts
    

# 2.b. Use the script "spimap-train-duploss" to learn duplication and loss
# rates from gene counts "train/fungi.counts".

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
    -c train/fungi.counts > train/fungi.duploss


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


# run SPIMAP for family 100
spimap \
    -a fungi-fams/100/100.nt.align \
    -s config/fungi.stree \
    -S config/fungi.smap \
    -p train/fungi.params \
    -o fungi-fams/100/100 \
    -D 0.000564 \
    -L 0.003056 \
    -i 100 \
    --quickiter 1000 \
    -V 1 --log -



# example of how to run SPIMAP for many families
ls fungi-fams | (while read FAMID; do
  NAME=fungi-fams/$FAMID/$FAMID
  spimap \
    -a $NAME.nt.align \
    -s config/fungi.stree \
    -S config/fungi.smap \
    -p train/fungi.params \
    -o $NAME \
    -D 0.000564 \
    -L 0.003056 \
    -i 100 \
    --quickiter 1000 \
    -V 1
done)



#=============================================================================
# Clean up
#
# Remove all outputs from example analysis
#

rm -f train/fungi*
find fungi-fams/ | egrep '\.tree$|\.log$' | xargs rm

