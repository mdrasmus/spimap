SPIMAP (Spieces informed max a poseriori)
http://compbio.mit.edu/spimap/
Matthew Rasmussen

=============================================================================
ABOUT

SPIMAP is phylogenetic program that uses species information to aide in
reconstructing gene trees.  It uses code from the SPIDIR phylogenetic library.

SPIMAP citation: Rasmussen, Kellis.  A Bayesian Approach for Fast and
Accurate Gene-tree Reconstruction. In prep 2010.

SPIDIR citation:
Rasmussen, Kellis. Accurate gene-tree reconstruction by learning
gene- and species-specific substitution rates across multiple complete genomes.
Genome Research. 2007

This package includes the C++ source of the SPIMAP program and SPIDIR library
as well as several library interfaces for C and python.


=============================================================================
INSTALL

NOTE: Makefile is installation will work best on UNIX or CYGWIN.


To compile the SPIMAP stand-alone program use the Makefile.

    make

To compile the SPIDIR C-library use:
    
    make libspidir    

To compile the SPIDIR python interface use:

    make pyspidir.so


Once compiled, to install the SPIMAP program (installs by default in /usr) use:

    make install

To specify your own installation path use:
    
    make install prefix=/usr/local


=============================================================================
USAGE

Running spimap with no arguments will print out its command-line usage:


Usage: spimap [OPTION]

  -a,--align  <alignment fasta>
    sequence alignment in fasta format

  -S,--smap  <species map>
    gene to species map

  -s,--stree  <species tree>
    species tree file in newick format

  -p,--param  <spidir params file>
    SPIDIR branch length parameters file

  -o,--output  <output filename prefix>
    prefix for all output filenames

Sequence model evolution
  -l,--lengths  (hky|parsimony)
    algorithm for determining branch lengths

  -r,--tsvratio  <transition/transversion ratio>
    used for HKY model (default=0.5)

  -f,--bgfreq  <A freq>,<C ferq>,<G freq>,<T freq>
    background frequencies (default=0.25,0.25,0.25,0.25

Miscellaneous
  -i,--niter  <# iterations>
    number of iterations

  -D,--dupprob  <duplication probability>
    probability of a node being a duplication (default=1.0)

  -P,--predupprob  <pre-duplication probability>
    probability of a node being a pre-duplication (default=0.01)

  -V,--verbose  <verbosity level>
    verbosity level 0=quiet, 1=low, 2=medium, 3=high

  ,--log  <log filename>
    log filename.  Use '-' to display on stdout.

  -v,--version  
    display version information

  -h,--help  
    display help information
