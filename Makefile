OBJS = main.o SurfRemesh.o

CPPFLAGS = -O3 -fPIC -I$(TRIANGLE_DIR)/include -fopenmp

LIBS = $(TRIANGLE_DIR)/lib/triangle.o
LIBS += -lgomp

remesh:$(OBJS)
	g++ -o remesh $(OBJS) $(LIBS)

.o:.cpp
	g++ $(CPPFLAGS) $<

clean:
	\rm -rf *.o remesh

