all: main

main: geometry.o quadTables.o main.cpp
	g++ -I /usr/local/include/eigen geometry.cpp quadTables.cpp main.cpp -o main

geometry.o: geometry.cpp
	g++ -I /usr/local/include/eigen -c geometry.cpp

quadTables.o: quadTables.cpp
	g++ -I /usr/local/include/eigen -c quadTables.cpp

testQuadTables: quadTables.cpp scriptTestQuadrature.cpp
	g++ -I /usr/local/include/eigen quadTables.cpp scriptTestQuadrature.cpp -o scriptTestQuadrature

testLegendre: scriptTestLegendre.cpp
	g++ -L /usr/lib -std=c++17 scriptTestLegendre.cpp -o scriptTestLegendre

clean:
	rm *.o
