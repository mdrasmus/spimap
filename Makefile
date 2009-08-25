#
# SPIDIR (SPecies Informed DIstance-based Reconstruction) 
# Matt Rasmussen
# Copyright 2007-2009
#
# Makefile
#

# install prefix paths
prefix = /usr


# C++ compiler options
CXX = g++

CFLAGS := $(CFLAGS) \
    -Wall -fPIC \
    -Isrc


# matlab options
MEX = mex
MATLAB_DIR = /afs/csail/i386_linux24/matlab/2007a
MATLAB_CFLAGS = \
    -g \
    -I$(MATLAB_DIR)/extern/include/cpp 
MEX_EXT = mexglx

#=============================================================================
# optional CFLAGS

# profiling
ifdef PROFILE
	CFLAGS := $(CFLAGS) -pg
endif

# debugging
ifdef DEBUG
	CFLAGS := $(CFLAGS) -g
else
	CFLAGS := $(CFLAGS) -O3
endif


#=============================================================================
# SPIDIR program files

# program files
SPIDIR_PROG = bin/spidir
SPIDIR_DEBUG = bin/spidir-debug

SPIDIR_SRC = \
    src/spidir.cpp \
    src/seq_likelihood.cpp \
    src/hky.cpp \
    src/common.cpp \
    src/branch_prior.cpp \
    src/birthdeath.cpp \
    src/parsimony.cpp \
    src/phylogeny.cpp \
    src/search.cpp \
    src/Sequences.cpp \
    src/Tree.cpp \
    src/train.cpp
#    src/branch_sample.cpp 

SPIDIR_OBJS = $(SPIDIR_SRC:.cpp=.o)

PROG_SRC = src/spidir_main.cpp 
PROG_OBJS = src/spidir_main.o $(SPIDIR_OBJS)
PROG_LIBS = -lgsl -lgslcblas -lm
#`gsl-config --libs`

#=======================
# SPIDIR C-library files
LIBSPIDIR = lib/libspidir.a
LIBSPIDIR_SHARED = lib/libspidir.so
LIBSPIDIR_OBJS = $(SPIDIR_OBJS)


#====================
# SPIDIR matlab files 
MATLAB_OBJS = matlab/spidir_treelk.$(MEX_EXT) \
              matlab/spidir_display_tree.$(MEX_EXT) \
              matlab/spidir_genbranches.$(MEX_EXT) \
              matlab/spidir_mlhkydist.$(MEX_EXT) \
              matlab/spidir_neighborjoin.$(MEX_EXT) \
              matlab/spidir_reconcile.$(MEX_EXT) \
              matlab/spidir_readtree.$(MEX_EXT)

MATLAB_COMPILE = spidir_matlab_compile.m
MATLAB_COMPILE_RULES = \
              matlab/spidir_treelk.rule \
              matlab/spidir_display_tree.rule \
              matlab/spidir_genbranches.rule \
              matlab/spidir_mlhkydist.rule \
              matlab/spidir_neighborjoin.rule \
              matlab/spidir_reconcile.rule \
              matlab/spidir_readtree.rule

MATLAB_SRC = $(SPIDIR_SRC) src/matlab_interface.cpp



#=============================================================================
# targets

# default targets
all: $(SPIDIR_PROG) $(LIBSPIDIR) $(LIBSPIDIR_SHARED)

debug: $(SPIDIR_DEBUG)

# SPIDIR stand-alone program
$(SPIDIR_PROG): $(PROG_OBJS) 
	$(CXX) $(CFLAGS) $(PROG_OBJS) $(PROG_LIBS) -o $(SPIDIR_PROG)

$(SPIDIR_DEBUG): $(PROG_OBJS) 
	$(CXX) $(CFLAGS) $(PROG_OBJS) $(PROG_LIBS) -o $(SPIDIR_DEBUG)


#-----------------------------
# maximum likelihood program
maxml: maxml.o $(SPIDIR_OBJS)
	$(CXX) $(CFLAGS) maxml.o $(SPIDIR_OBJS) $(PROG_LIBS) -o maxml

#-----------------------------
# SPIDIR C-library
lib: $(LIBSPIDIR) $(LIBSPIDIR_SHARED)

$(LIBSPIDIR): $(LIBSPIDIR_OBJS)
	mkdir -p lib
	$(AR) -r $(LIBSPIDIR) $(LIBSPIDIR_OBJS)

$(LIBSPIDIR_SHARED): $(LIBSPIDIR_OBJS) 
	mkdir -p lib
	$(CXX) -o $(LIBSPIDIR_SHARED) -shared $(LIBSPIDIR_OBJS) $(PROG_LIBS)


#------------------------------
# SPIDIR matlab interface
matlab: $(MATLAB_OBJS) $(MATLAB_COMPILE)


$(MATLAB_OBJS): %.$(MEX_EXT): %.cpp
	$(MEX) $(MATLAB_CFLAGS) $(MATLAB_SRC) $< -o $@

# generate compile rules for windows
$(MATLAB_COMPILE): $(MATLAB_COMPILE_RULES)
$(MATLAB_COMPILE_RULES): %.rule: %.cpp
	echo "display('compiling $<...');" >> $(MATLAB_COMPILE)
	echo $(MEX) $(MATLAB_SRC) $< -o $(@:%.rule=%) >> $(MATLAB_COMPILE)
	touch $@



#=============================================================================
# basic rules

$(SPIDIR_OBJS): %.o: %.cpp
	$(CXX) -c $(CFLAGS) -o $@ $<


install: $(SPIDIR_PROG)
	cp $(SPIDIR_PROG) $(prefix)/bin


myinstall: $(SPIDIR_PROG) maxml
	cp $(SPIDIR_PROG) maxml ../bin


myinstall64: $(SPIDIR_PROG) maxml
	cp $(SPIDIR_PROG) maxml ../bin64


clean:
	rm -f $(PROG_OBJS) $(SPIDIR_PROG) $(LIBSPIDIR) \
              $(MATLAB_OBJS) maxml maxml.o \
              $(MATLAB_COMPILE) $(MATLAB_COMPILE_RULES)

clean-obj:
	rm -f $(PROG_OBJS)


#=============================================================================
# dependencies

dep:
	touch Makefile.dep
	makedepend -f Makefile.dep src/*.cpp src/*.h

Makefile.dep:
	touch Makefile.dep
	makedepend -f Makefile.dep src/*.cpp src/*.h

include Makefile.dep
