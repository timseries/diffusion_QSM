#default compiler is mpicxx
#use make bluegene for normal compilation
#use make bluegene profile to include profiling
CC		= $(PREP) mpicxx
CPPFLAGS	= -c -DSINGLE_PRECISION
CPPFLAGS 	+= -Iinclude
LDFLAGS		=
LIBS		= 
SOURCES 	=

ifdef hybrid
SOURCEDIR     	= src/hybrid
ifdef USE_OPENCL
CPPFLAGS 	+= -DUSE_OPENCL -msse4a --fast-math
LDFLAGS 	+= -lm -lOpenCL
SOURCES		+= $(SOURCEDIR)/opencl_base.cc
endif
SOURCES		+= $(SOURCEDIR)/arghandler.cc $(SOURCEDIR)/modelmap.cc $(SOURCEDIR)/dataspec.cc
SOURCES		+= $(SOURCEDIR)/output.cc $(SOURCEDIR)/util.cc
SOURCES		+= $(SOURCEDIR)/kernel.cc $(SOURCEDIR)/problem.cc
SOURCES		+= $(SOURCEDIR)/process.cc $(SOURCEDIR)/hybrid_main.cc
else
SOURCEDIR     	= src/legacy
SOURCES		= $(SOURCEDIR)/chimap.cc
endif

EXECUTABLE	= $(SOURCEDIR)/dqsm
OBJECTS		= $(SOURCES:.cc=.o)


ifdef omp
CPPFLAGS += -DUSE_OPENMP
ifdef bluegene
CPPFLAGS += -qsmp=omp
LDFLAGS += -qsmp=omp
else
CPPFLAGS += -fopenmp
LDFLAGS += -fopenmp
endif
endif 

ifdef debug
#CPPFLAGS +=  -Wall -g
CPPFLAGS +=  -g
else
CPPFLAGS +=  -O3
endif 

ifdef bluegene
CC = mpixlcxx_r
CPPFLAGS += -qstrict
ifdef hpm_profile
CPPFLAGS += -g
CPPFLAGS += -I/bgsys/ibmhpc/ppedev.hpct/include/
LDFLAGS += -L/bgsys/drivers/ppcfloor/bgpm/lib/
ifdef omp
LDFLAGS += -L/bgsys/ibmhpc/ppedev.hpct/lib64 -lhpc_r -lbgpm
else
LDFLAGS += -L/bgsys/ibmhpc/ppedev.hpct/lib64 -lhpc -lbgpm -qsmp=omp
endif
ifdef fourier_spheres
CPPFLAGS += -I/usr/local/fftw/3.3.3-xl/include
endif
CPPFLAGS += -DHPM
endif

ifdef mpi_profile
CPPFLAGS += -g
CPPFLAGS += -I/bgsys/ibmhpc/ppedev.hpct/include/
#LDFLAGS += -L/bgsys/drivers/ppcfloor/bgpm/lib/
LDFLAGS += -L/bgsys/ibmhpc/ppedev.hpct/lib64 -lmpitrace
CPPFLAGS += -DMPI_PROFILE
endif

ifdef pomp_profile
CPPFLAGS += -g
CPPFLAGS += -I/bgsys/ibmhpc/ppedev.hpct/include/
LDFLAGS += -L/opt/ibmcmp/xlsmp/3.1/lib64 -lxlsmp_pomp
LDFLAGS += -L/bgsys/ibmhpc/ppedev.hpct/lib64 -lpomprof_probe
CPPFLAGS += -DOPENMP_PROFILE #never used in the code, just for consistency...
endif
endif

ifdef gmon_profile
CPPFLAGS += -pg
LDFLAGS += -pg
endif

ifdef fourier_spheres
CPPFLAGS += -DUSE_FOURIER_SPHERES
LDFLAGS += -lfftw3f -lm
endif


all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS) $(LIBS) 
	rm -rf $(SOURCEDIR)/*.o

.cc.o:
	$(CC) $(CPPFLAGS) $< -o $@

clean:
	rm -rf $(SOURCEDIR)/*.o $(EXECUTABLE)




