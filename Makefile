all: gate.o ring.o main.o
	g++ -I/usr/include/lapackpp -llapackpp -o test gate.o ring.o main.o

ring.o: ring.cpp
	g++ -c -I/usr/include/lapackpp -llapackpp ring.cpp

gate.o: gate.cpp
	g++ -c -I/usr/include/lapackpp -llapackpp gate.cpp

main.o: main.cpp
	g++ -c -I/usr/include/lapackpp -llapackpp main.cpp

clean: 
	rm *.o
