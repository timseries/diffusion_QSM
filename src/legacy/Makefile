#bluegene
#CC		= mpixlcxx_r 
#pc
CC		= mpicxx
#use -gp for both of these to enable profiling
CPPFLAGS	= -c
LDFLAGS		= 
LIBS		= 
SOURCES		= chimap.cc
OBJECTS		= $(SOURCES:.cc=.o)
EXECUTABLE	= chimap

all: $(SOURCES) $(EXECUTABLE)

profile: CPPFLAGS += -pg
profile: LDFLAGS += -pg
profile: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS) $(LIBS) 
	rm -rf *.o

.cc.o:
	$(CC) $(CPPFLAGS) $< -o $@

clean:
	rm -rf *.o $(EXECUTABLE)




