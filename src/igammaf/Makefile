
CXX=g++
CFLAGS=-I libs

OBJS = gammaf.o \
       igammaf.o \
       normaldistr.o \
       libs/ap.o

igammaf.a: $(OBJS)
	$(AR) -r igammaf.a $(OBJS)

$(OBJS): %.o: %.cpp
	$(CXX) -c $(CFLAGS) -o $@ $<


