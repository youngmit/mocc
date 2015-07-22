OBJS := pinmesh_cyl.o coremesh.o inputproc.o
CXX := g++
CXXFLAGS := -g
INCLUDES := -I. -I../lib/pugixml/src
CP := copy

mocc.exe: mocc.cpp Makefile pugixml.o $(OBJS)
	$(CXX) $(CXXFLAGS) mocc.cpp $(OBJS) pugixml.o -o mocc.exe
	$(CP) "./mocc.exe" "../bin/"
	
pugixml.o: ../lib/pugixml/src/pugixml.cpp
	$(CXX) $(INCLUDES) $(CXXFLAGS) -c $< -o $@

%.o: %.cpp
	$(CXX) $(INCLUDES) $(CXXFLAGS) -c $< -o $@
	
