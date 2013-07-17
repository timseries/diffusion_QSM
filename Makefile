#default compiler is mpicxx
#use make bluegene for normal compilation
#use make bluegene profile to include profiling
CC		= mpicxx
CPPFLAGS	= -c -O3
CPPFLAGS 	+= -Iinclude
LDFLAGS		=
LIBS		= 

ifdef hybrid
SOURCEDIR     	= src/hybrid
SOURCES		= $(SOURCEDIR)/arghandler.cc $(SOURCEDIR)/modelmap.cc
SOURCES		+= $(SOURCEDIR)/output.cc $(SOURCEDIR)/util.cc
SOURCES		+= $(SOURCEDIR)/kernel.cc $(SOURCEDIR)/problem.cc
SOURCES		+= $(SOURCEDIR)/process.cc $(SOURCEDIR)/hybrid_main.cc
else
SOURCEDIR     	= src/legacy
SOURCES		= $(SOURCEDIR)/chimap.cc
endif

EXECUTABLE	= $(SOURCEDIR)/dqsm
OBJECTS		= $(SOURCES:.cc=.o)

ifdef openmp
CPPFLAGS += -DUSE_OPENMP
ifdef bluegene
CPPFLAGS += -qtm -qsmp=omp
LDFLAGS += -qtm -qsmp=omp
else
CPPFLAGS += -fopenmp
LDFLAGS += -fopenmp
endif
endif 

ifdef bluegene
CC = mpixlcxx_r

ifdef hpm_profile
CPPFLAGS += -g
CPPFLAGS += -I/bgsys/ibmhpc/ppedev.hpct/include/
LDFLAGS += -L/bgsys/drivers/ppcfloor/bgpm/lib/
LDFLAGS += -L/bgsys/ibmhpc/ppedev.hpct/lib64 -lhpc_r -lbgpm
#LDFLAGS += -L/bgsys/ibmhpc/ppedev.hpct/lib64 -lhpc -lbgpm
CPPFLAGS += -DHPM
endif

ifdef mpi_profile
CPPFLAGS += -g
CPPFLAGS += -I/bgsys/ibmhpc/ppedev.hpct/include/
#LDFLAGS += -L/bgsys/drivers/ppcfloor/bgpm/lib/
LDFLAGS += -L/bgsys/ibmhpc/ppedev.hpct/lib64 -lmpitrace
CPPFLAGS += -DMPI_PROFILE
endif

ifdef openmp_profile
CPPFLAGS += -g
CPPFLAGS += -I/bgsys/ibmhpc/ppedev.hpct/include/
LDFLAGS += -L/opt/ibmcmp/xlsmp/3.1/lib64 -lxlsmp_pomp
LDFLAGS += -L/bgsys/ibmhpc/ppedev.hpct/lib64 -lpomprof_probe
CPPFLAGS += -DOPENMP_PROFILE #never used
endif
endif

ifdef gprof_profile
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
	rm -rf $(SOURCEDIR)/*.o $(EXECUTABLE)




