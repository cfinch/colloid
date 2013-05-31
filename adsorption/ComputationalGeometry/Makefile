# location of the Python header files
PYTHON_VERSION = 2.6
PYTHON_INCLUDE = /usr/include/python$(PYTHON_VERSION)

# location of the Boost Python include files and library

BOOST_INC = /usr/include
BOOST_LIB = /usr/lib

# compile mesh classes
HEADERS = surface_arrangement.h
SOURCES = surface_arrangement.C visualizer.C surface_arrangement_IO.C surface_arrangement_area.C
OBJECTS = surface_arrangement.o visualizer.o surface_arrangement_IO.o surface_arrangement_area.o
TARGET = surface.so

$(TARGET): $(OBJECTS)
	g++ -shared -Wl,--export-dynamic \
	$(OBJECTS) -L$(BOOST_LIB) -lboost_python -lCGAL \
	-L/usr/lib/python$(PYTHON_VERSION)/config -lpython$(PYTHON_VERSION) \
	-o $(TARGET)

$(OBJECTS): $(SOURCES) $(HEADERS)
	g++ -I$(PYTHON_INCLUDE) -I$(BOOST_INC) -fPIC -c $(SOURCES)

.PHONY : clean
clean : 
	rm $(OBJECTS) $(TARGET)
