
CXX=g++
CFLAGS=-I libs -fPIC

OBJS = gammaf.o \
       igammaf.o \
       normaldistr.o \
       libs/ap.o

igammaf.a: $(OBJS)
	$(AR) -r igammaf.a $(OBJS)

$(OBJS): %.o: %.cpp
	$(CXX) -c $(CFLAGS) -o $@ $<

clean:
	rm -rf $(OBJS) igammaf.a
