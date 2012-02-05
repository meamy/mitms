all:
	g++ -o test -I/usr/include/lapackpp -llapackpp gate.cpp
