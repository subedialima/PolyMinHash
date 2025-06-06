
#CC = mpicxx
CC =   /usr/bin/mpicxx
CFLAGS = -std=c++20 -g -O0 -Wall
IFLAGS = -I/usr/include
LFLAGS = -L/usr/lib/x86_64-linux-gnu -lgeos_c -lpthread 
LIBS = -lgeos_c -lpthread

all: main.o parse_geodata.o geoutil.o util.o query.o brute_force.o mpi_gis.o
	$(CC) $(CFLAGS) -o spjoin main.o parse_geodata.o geoutil.o util.o query.o brute_force.o mpi_gis.o $(LFLAGS)

main.o: main.cpp
	$(CC) $(CFLAGS) $(IFLAGS) $(LFLAGS) -c main.cpp 

brute_force.o: brute_force.cpp query.h
	$(CC) $(CFLAGS) $(IFLAGS) -c brute_force.cpp
	
mpi_gis.o: mpi_gis.cpp mpi_gis.h
	$(CC) $(CFLAGS) $(IFLAGS) -c mpi_gis.cpp $(LFLAGS)	

geoutil.o: geoutil.cpp geoutil.h
	$(CC) $(CFLAGS) $(IFLAGS) -c geoutil.cpp $(LFLAGS)

parse_geodata.o: parse_geodata.cpp parse_geodata.h
	$(CC) $(CFLAGS) $(IFLAGS) $(LFLAGS) -c parse_geodata.cpp 

util.o: util.cpp util.h
	$(CC) $(CFLAGS) $(IFLAGS) -c util.cpp $(LFLAGS)

query.o: query.cpp query.h
	$(CC) $(CFLAGS) $(IFLAGS) -c query.cpp $(LFLAGS)

clean:
	rm -f *.o spjoin
