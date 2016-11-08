CC	= gcc-6
CFLAGS	= -O0 -g -Wall -fopenmp
LDFLAGS = -fopenmp -lm -lpthread -L/usr/X11/lib -lX11 
VERBOSE	=
TARGET	= nbody_brute_force nbody_barnes_hut
OBJS	= ui.o xstuff.o nbody_tools.o nbody_alloc.o

#DISPLAY = -DDISPLAY
#DISPLAY =

all: $(TARGET)

nbody_brute_force: nbody_brute_force.o $(OBJS)
	$(CC) $(VERBOSE) -o $@ $< $(OBJS) $(LDFLAGS)

nbody_barnes_hut: nbody_barnes_hut.o $(OBJS)
	$(CC) $(VERBOSE) -o $@ $< $(OBJS)  $(LDFLAGS)



%.o: %.c
	$(CC) $(CFLAGS) -c $< $(VERBOSE) $(DISPLAY)
clean:
	rm -f *.o $(TARGET)
