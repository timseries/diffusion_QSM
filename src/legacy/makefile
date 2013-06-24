
#CC		= mpixlcxx_r
CC		= mpic++
LIBS		= 
SOURCES		= chimap.cpp
OBJECTS		= $(SOURCES:.cpp=.o)
EXECUTABLE	= chimap

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS) $(LIBS) 
	rm -rf *.o

.cpp.o:
	$(CC) -c $(CPPFLAGS) $< -o $@

clean:
	rm -rf *.o $(EXECUTABLE)




