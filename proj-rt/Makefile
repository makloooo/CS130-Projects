CXX = g++
CXXFLAGS = -Wall -g -O3 -std=c++11
LIBS = -lpng
SRCS = src/boolean.cpp src/camera.cpp src/dump_png.cpp src/flat_shader.cpp src/main.cpp src/parse.cpp src/phong_shader.cpp src/plane.cpp src/reflective_shader.cpp src/render_world.cpp src/sphere.cpp
OBJS := ${SRCS:.cpp=.o}

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	LIBS = -framework OpenGL -framework GLUT -lpng
	CXXFLAGS += -Wno-deprecated-declarations
endif

TARGET = ray_tracer

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $(TARGET) $(OBJS) $(LIBS)

clean:
	rm -f src/*.o src/*.d src/*~ src/*.bak $(TARGET)
	rm -rf grading

# Generate a list of dependencies for each cpp file
%.d: $(SRCS)
	@ $(CXX) $(CPPFLAGS) -MM $*.cpp | sed -e 's@^\(.*\)\.o:@\1.d \1.o:@' > $@

# At the end of the makefile
# Include the list of dependancies generated for each object file
# unless make was called with target clean
ifneq "$(MAKECMDGOALS)" "clean"
-include ${SRCS:.cpp=.d}
endif

