all: main

main: geometry.o quadTables.o main.cpp
	g++ -I /usr/local/include/eigen geometry.cpp quadTables.cpp main.cpp -o main

geometry.o: geometry.cpp
	g++ -I /usr/local/include/eigen -c geometry.cpp

quadTables.o: quadTables.cpp
	g++ -I /usr/local/include/eigen -c quadTables.cpp

legendrebasis.o: legendrebasis.cpp
	g++ -I /usr/local/include/eigen -L /usr/lib -c legendrebasis.cpp

matrixRoutines.o: matrixRoutines.cpp
	g++ -I /usr/local/include/eigen -L /usr/lib -c matrixRoutines.cpp

testsAndProjections.o: testsAndProjections.cpp
	g++ -I /usr/local/include/eigen -L /usr/lib -c testsAndProjections.cpp

# scripts down here

testQuadTables: quadTables.cpp scriptTestQuadrature.cpp
	g++ -I /usr/local/include/eigen quadTables.cpp scriptTestQuadrature.cpp -o scriptTestQuadrature

testLegendre: scriptTestLegendre.cpp
	g++ -L /usr/lib scriptTestLegendre.cpp -o scriptTestLegendre

scriptTestLegendreBasis: legendrebasis.o scriptTestLegendreBasis.cpp
	g++ -I /usr/local/include/eigen -L /usr/lib legendrebasis.cpp scriptTestLegendreBasis.cpp -o scriptTestLegendreBasis

scriptTestMatrixRoutines: matrixRoutines.o scriptTestMatrixRoutines.cpp
	g++ -I /usr/local/include/eigen -L /usr/lib matrixRoutines.cpp scriptTestMatrixRoutines.cpp -o scriptTestMatrixRoutines

scriptTestTestsAndProjecions: testsAndProjections.o quadTables.o matrixRoutines.o legendrebasis.o geometry.o
	 g++ -I /usr/local/include/eigen -L /usr/lib matrixRoutines.cpp quadTables.cpp legendrebasis.cpp geometry.cpp  -o scriptTestMatrixRoutines


clean:
	rm *.o
