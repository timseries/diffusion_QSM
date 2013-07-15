#default compiler is mpicxx
#use make bluegene for normal compilation
#use make bluegene profile to include profiling
CC		= mpicxx
CPPFLAGS	= -c -o3
CPPFLAGS 	+= -Iinclude
LDFLAGS		= 
LIBS		= 

ifdef hybrid
SOURCEDIR     	= src/hybrid
SOURCES		= $(SOURCEDIR)/arghandler.cc $(SOURCEDIR)/modelmap.cc
SOURCES		+= $(SOURCEDIR)/output.cc $(SOURCEDIR)/util.cc
SOURCES		+= $(SOURCEDIR)/kernel.cc  $(SOURCEDIR)/process.cc
SOURCES		+= $(SOURCEDIR)/hybrid_main.cc
#SOURCES		= $(SOURCEDIR)/hybrid_main.cc
#SOURCES		+= $(SOURCEDIR)/process.cc $(SOURCEDIR)/kernel.cc
#SOURCES		+= $(SOURCEDIR)/modelmap.cc $(SOURCEDIR)/output.cc
#SOURCES		+= $(SOURCEDIR)/arghandler.cc $(SOURCEDIR)/util.cc

else
SOURCEDIR     	= src/legacy
SOURCES		= $(SOURCEDIR)/chimap.cc
endif

EXECUTABLE	= $(SOURCEDIR)/dqsm
OBJECTS		= $(SOURCES:.cc=.o)

ifdef bluegene
CC = mpixlcxx_r
endif

ifdef hpm
#CC = mpixlcxx
CPPFLAGS += -g
CPPFLAGS += -I/bgsys/ibmhpc/ppedev.hpct/include/
LDFLAGS += -L/bgsys/drivers/ppcfloor/bgpm/lib/
LDFLAGS += -L/bgsys/ibmhpc/ppedev.hpct/lib64 -lhpc_r -lbgpm -qsmp=omp
#LDFLAGS += -L/bgsys/ibmhpc/ppedev.hpct/lib64 -lhpc -lbgpm
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
	rm -rf $(SOURCEDIR)/*.o

.cc.o:
	$(CC) $(CPPFLAGS) $< -o $@

clean:
	rm -rf $(SOURCEDIR)/*.o $(SOURCEDIR)/$(EXECUTABLE)




