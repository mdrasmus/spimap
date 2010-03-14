# 
# Make commands for matlab
#
# Sun Mar 14 14:50:08 EDT 2010
# The matlab interface is currently unmantained but is included here incase
# maintaince is resumed.
#



# matlab options
MEX = mex
MATLAB_DIR = /afs/csail/i386_linux24/matlab/2007a
MATLAB_CFLAGS = \
    -g \
    -I$(MATLAB_DIR)/extern/include/cpp 
MEX_EXT = mexglx


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


