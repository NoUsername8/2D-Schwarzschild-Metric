# Files
S_FILES=Simulation.cpp

# Output
EXEC=simul

# Build settings
CC=g++
# SDL options
CC_SDL=-lSDL2
CC_TTF=-lSDL2_ttf
CC_GMP=-lgmpxx -lgmp

all:Build

Build:
	g++ -c math/Polar.cpp
	g++ -c math/Matrix.cpp
	g++ -c Metric.cpp
	g++ -c ScreenHelper.cpp
	$(CC) $(S_FILES) $(CC_GMP) -w $(CC_SDL) -w $(CC_TTF) -o $(EXEC) Matrix.o Polar.o ScreenHelper.o Metric.o

build_run:Build
	$(EXEC)

clean:
	rm -rf simul
