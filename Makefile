CC = gcc
CFLAGS = -Wall -O2
INCLUDES = -I/usr/include/htslib -I.
LIBS = -lhts -lz -lpthread

fastMPRA: fastMPRA.c
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $< $(LIBS)

clean:
	rm -f fastMPRA
