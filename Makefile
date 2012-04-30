FLAGS = -O3 -funroll-loops -I/usr/include/lapackpp/ -llapackpp
OBJS = ring.o matrix.o util.o gate.o circuit.o search.o

all: $(OBJS)
	g++ $(FLAGS) -o opt $(OBJS)

ring.o: ring.cpp
	g++ -c $(FLAGS) ring.cpp

matrix.o: matrix.cpp
	g++ -c $(FLAGS) matrix.cpp

util.o: util.cpp util.h
	g++ -c $(FLAGS) util.cpp

gate.o: gate.cpp
	g++ -c $(FLAGS) gate.cpp

circuit.o : circuit.cpp
	g++ -c $(FLAGS) circuit.cpp

search.o: search.cpp util.h
	g++ -c $(FLAGS) search.cpp

clean: 
	rm *.o
