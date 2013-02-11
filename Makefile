FLAGS = -std=c++0x -O3 -fpermissive -funroll-loops -L/usr/local/lib/ -L/usr/lib/ -I/usr/local/include/lapackpp/ -I/usr/include/lapackpp/ -llapackpp -lrt -pthread
OBJS = configs.o ring.o matrix.o gate.o circuit.o util.o database.o search.o main.o
MAINS = main.o

all: $(OBJS)
	g++ -o mitms $(OBJS) $(FLAGS)

configs.o: configs.cpp
	g++ -c $(FLAGS) configs.cpp

ring.o: ring.cpp
	g++ -c $(FLAGS) ring.cpp

matrix.o: matrix.cpp
	g++ -c $(FLAGS) matrix.cpp

gate.o: gate.cpp
	g++ -c $(FLAGS) gate.cpp

circuit.o: circuit.cpp
	g++ -c $(FLAGS) circuit.cpp

util.o: util.cpp
	g++ -c $(FLAGS) util.cpp

database.o: database.cpp
	g++ -c $(FLAGS) database.cpp

search.o: search.cpp
	g++ -c $(FLAGS) search.cpp

main.o: main.cpp
	g++ -c $(FLAGS) main.cpp 

clean: 
	rm *.o
