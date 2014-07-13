FLAGS = -std=c++0x -O3 -fpermissive -funroll-loops -I/usr/local/include/lapackpp/ -I/usr/include/lapackpp/ -llapackpp -lrt -pthread
OBJS = configs.o ring.o matrix.o gate.o circuit.o util.o database.o search.o main.o
MAINS = main.o

all: $(OBJS)
	g++ -o mitms $(OBJS) $(FLAGS)

configs.o: src/configs.cpp
	g++ -c $(FLAGS) src/configs.cpp

ring.o: src/ring.cpp
	g++ -c $(FLAGS) src/ring.cpp

matrix.o: src/matrix.cpp
	g++ -c $(FLAGS) src/matrix.cpp

gate.o: src/gate.cpp
	g++ -c $(FLAGS) src/gate.cpp

circuit.o: src/circuit.cpp
	g++ -c $(FLAGS) src/circuit.cpp

util.o: src/util.cpp
	g++ -c $(FLAGS) src/util.cpp

database.o: src/database.cpp
	g++ -c $(FLAGS) src/database.cpp

search.o: src/search.cpp
	g++ -c $(FLAGS) src/search.cpp

main.o: src/main.cpp
	g++ -c $(FLAGS) src/main.cpp 

clean: 
	rm *.o
