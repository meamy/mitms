FLAGS = -O3 -I/usr/include/lapackpp/ -llapackpp
OBJS = ring.o matrix.o gate.o search.o

all: $(OBJS)
	g++ $(FLAGS) -o opt $(OBJS)

ring.o: ring.cpp
	g++ -c $(FLAGS) ring.cpp

matrix.o: matrix.cpp
	g++ -c $(FLAGS) matrix.cpp

gate.o: gate.cpp
	g++ -c $(FLAGS) gate.cpp

search.o: search.cpp
	g++ -c $(FLAGS) search.cpp

clean: 
	rm *.o
