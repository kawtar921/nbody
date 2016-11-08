CC	=  gcc
CFLAGS	= -O0 -g -Wall -fopenmp
LDFLAGS = -pg -fopenmp -lm -lpthread -L/usr/X11/lib -lX11 
VERBOSE	=
TARGET	= nbody_brute_force 
OBJS	= ui.o xstuff.o nbody_tools.o nbody_alloc.o
USRLIBS   = -L/usr/X11/lib -lX11

#DISPLAY = -DDISPLAY
#DISPLAY =

all: program

program: nbody_brute_force.o $(OBJS)
	mpicc $(USRLIBS) nbody_brute_force.o $(OBJS) -o nbody_brute_force

nbody_brute_force: nbody_brute_force.c
	mpicc -c nbody_brute_force.o nbody_brute_force.c

ui.o: ui.c
	$(CC) $(CFLAGS) -c $< $(VERBOSE) $(DISPLAY)

xstuff.o: xstuff.c
	$(CC) $(CFLAGS) -c $< $(VERBOSE) $(DISPLAY)

body_tools.o: body_tools.c
	$(CC) $(CFLAGS) -c $< $(VERBOSE) $(DISPLAY)

nbody_alloc.o: nbody_alloc.c
	$(CC) $(CFLAGS) -c $< $(VERBOSE) $(DISPLAY)


clean:
	rm -f *.o $(TARGET)


