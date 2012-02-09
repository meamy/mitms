all: gate.o main.o
	g++ -I/usr/include/lapackpp -llapackpp -o test gate.o main.o

gate.o: gate.cpp
	g++ -c -I/usr/include/lapackpp -llapackpp gate.cpp

main.o: main.cpp
	g++ -c -I/usr/include/lapackpp -llapackpp main.cpp

clean: 
	rm *.o
