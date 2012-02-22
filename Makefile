all: ring.o matrix.o gate.o search.o
	g++ -I/usr/include/lapackpp -llapackpp -o test gate.o ring.o matrix.o search.o

ring.o: ring.cpp
	g++ -c -I/usr/include/lapackpp -llapackpp ring.cpp

matrix.o: matrix.cpp
	g++ -c -I/usr/include/lapackpp -llapackpp matrix.cpp

gate.o: gate.cpp
	g++ -c -I/usr/include/lapackpp -llapackpp gate.cpp

search.o: search.cpp
	g++ -c -I/usr/include/lapackpp -llapackpp search.cpp

clean: 
	rm *.o
