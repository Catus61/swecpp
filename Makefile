CC=g++
CFLAGS=-g -l /home/jason/adolc_base/lib64
LDFLAGS=-llapack -lgsl -gslcblas -lblas -lm -ladolc -L /home/jason/adolc_base/lib64 -Wl,-rpath,/home/jason/adolc_base/lib64
BINARY=shallow_water

SOURCE_FILES=shallow_water.c\
             f_calc.c\
             init.c\
             io.c\
             jacobiandelay.c\
			 jacobiandrifter.c\
             drifter.c\
			 svdnudging.c

OBJECTS=$(SOURCE_FILES:.c=.o)

all: $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $(BINARY)

$(OBJECTS):
	$(CC) $(CFLAGS) -o $@ -c $(@:.o=.c)


clean:
	rm $(OBJECTS) $(BINARY)
