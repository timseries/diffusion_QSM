#default compiler is mpicxx
#use make bluegene for normal compilation
#use make bluegene profile to include profiling
CC		= mpicxx
CPPFLAGS	= -c
LDFLAGS		= 
LIBS		= 
SOURCES		= chimap.cc
OBJECTS		= $(SOURCES:.cc=.o)
EXECUTABLE	= chimap


ifdef bluegene
CC = mpixlcxx_r
endif

ifdef hpm
CPPFLAGS += -g
CPPFLAGS += -I/bgsys/ibmhpc/ppedev.hpct/include/
LDFLAGS += -L/bgsys/drivers/ppcfloor/bgpm/lib/
LDFLAGS += -L/bgsys/ibmhpc/ppedev.hpct/lib64 -lhpc_r -lbgpm -qsmp=omp
CPPFLAGS += -DHPM
endif

ifdef mpitrace
CPPFLAGS += -g
CPPFLAGS += -I/bgsys/ibmhpc/ppedev.hpct/include/
LDFLAGS += -L/bgsys/drivers/ppcfloor/bgpm/lib/
LDFLAGS += -L/bgsys/ibmhpc/ppedev.hpct/lib64 -lmpitrace
CPPFLAGS += -Dmpitrace
endif

ifdef profile
CPPFLAGS += -pg
LDFLAGS += -pg
endif


all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS) $(LIBS) 
	rm -rf *.o

.cc.o:
	$(CC) $(CPPFLAGS) $< -o $@

clean:
	rm -rf *.o $(EXECUTABLE)




