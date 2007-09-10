#
# SPIDIR (SPecies Informed DIstance-based Reconstruction) 
# Matt Rasmussen
# copyright 2007
#
# Makefile
#


CXX = g++
MEX = mex

CFLAGS = \
    -Wall \
    -I/usr/include/python2.4 \
    -I/util/include/python2.4


# matlab options
MATLAB_DIR = /afs/csail/i386_linux24/matlab/6.5r13sp1
MATLAB_CFLAGS = \
    -g \
    -I$(MATLAB_DIR)/extern/include/cpp \


#=============================================================================
# target names

# program files
SPIDIR_PROG = spidir

SPIDIR_SRC = \
    spidir.cpp \
    mldist.cpp \
    common.cpp \
    likelihood.cpp \
    parsimony.cpp \
    phylogeny.cpp \
    search.cpp \
    Tree.cpp \

SPIDIR_OBJS = \
    spidir.o \
    mldist.o \
    common.o \
    likelihood.o \
    parsimony.o \
    phylogeny.o \
    search.o \
    Tree.o \
    Sequences.o 

PROG_SRC = spidir_main.cpp 
PROG_OBJS = spidir_main.o $(SPIDIR_OBJS)
PROG_LIBS = 


# C-library files
LIBSPIDIR = lib/libspidir.a
LIBSPIDIR_OBJS = $(SPIDIR_OBJS)


# python files
PYTHON_MODULE = pyspidir.so
PYTHON_MODULE_OBJS = \
    pyspidir.o \
    $(SPIDIR_OBJS)


# matlab files
MATLAB_FUNCS = $(MATLAB_TREELK_FUNC) \
               $(MATLAB_DISPLAY_TREE_FUNC) \
               $(MATLAB_GENBRANCHES_FUNC)
               
MATLAB_OBJS = $(MATLAB_TREELK_OBJ) \
              $(MATLAB_DISPLAY_TREE_OBJ) \
              $(MATLAB_GENBRANCHES_OBJ)


MATLAB_TREELK_FUNC = spidir_treelk
MATLAB_TREELK_SRC = $(SPIDIR_SRC) matlab_interface.cpp matlab_treelk.cpp
MATLAB_TREELK_OBJ = spidir_treelk.mexglx

MATLAB_DISPLAY_TREE_FUNC = spidir_display_tree
MATLAB_DISPLAY_TREE_SRC = $(SPIDIR_SRC) matlab_interface.cpp matlab_display_tree.cpp
MATLAB_DISPLAY_TREE_OBJ = spidir_display_tree.mexglx

MATLAB_GENBRANCHES_FUNC = spidir_genbranches
MATLAB_GENBRANCHES_SRC = $(SPIDIR_SRC) matlab_interface.cpp matlab_genbranches.cpp
MATLAB_GENBRANCHES_OBJ = spidir_genbranches.mexglx



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
# targets

all: $(SPIDIR_PROG) $(LIBSPIDIR) $(PYTHON_MODULE) test_spidir

# stand-alone program
$(SPIDIR_PROG): $(PROG_OBJS)
	$(CXX) $(CFLAGS) $(PROG_OBJS) $(PROG_LIBS) -o $(SPIDIR_PROG)

maxml: maxml.o $(SPIDIR_OBJS)
	$(CXX) $(CFLAGS) maxml.o $(SPIDIR_OBJS) $(PROG_LIBS) -o maxml

# C-library
$(LIBSPIDIR): $(LIBSPIDIR_OBJS)
	mkdir -p lib
	$(AR) -r $(LIBSPIDIR) $(LIBSPIDIR_OBJS)

# testing program
test_spidir: $(SPIDIR_OBJS) test.o
	$(CXX) $(SPIDIR_OBJS) $(CFLAGS) test.o $(PROG_LIBS) -o test_spidir

# python module
$(PYTHON_MODULE): $(PYTHON_MODULE_OBJS)
	$(CXX) -shared $(PYTHON_MODULE_OBJS) -o $(PYTHON_MODULE)


# matlab interface
matlab_funcs: $(MATLAB_FUNCS)

# matlab treelk function
$(MATLAB_TREELK_FUNC): $(MATLAB_TREELK_SRC)
	$(MEX) $(MATLAB_CFLAGS) $(MATLAB_TREELK_SRC) -o $(MATLAB_TREELK_FUNC)

# matlab display tree function
$(MATLAB_DISPLAY_TREE_FUNC): $(MATLAB_DISPLAY_TREE_SRC)
	$(MEX) $(MATLAB_CFLAGS) $(MATLAB_DISPLAY_TREE_SRC) -o $(MATLAB_DISPLAY_TREE_FUNC)

$(MATLAB_GENBRANCHES_FUNC): $(MATLAB_GENBRANCHES_SRC)
	$(MEX) $(MATLAB_CFLAGS) $(MATLAB_GENBRANCHES_SRC) -o $(MATLAB_GENBRANCHES_FUNC)


#=============================================================================
# basic rules

$(PYTHON_MODULE_OBJS): %.o: %.cpp
	$(CXX) -c $(CFLAGS) -o $@ $<


install: $(SPIDIR_PROG) $(PYTHON_MODULE) test_spidir maxml
	cp $(SPIDIR_PROG) test_spidir maxml ../bin
	cp $(PYTHON_MODULE) ../python


clean:
	rm -rf $(PROG_OBJS) $(SPIDIR_PROG) $(LIBSPIDIR) \
                $(PYTHON_MODULE_OBJS) $(PYTHON_MODULE) \
                $(MATLAB_OBJS) \
	        test.o test_spidir
