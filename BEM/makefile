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

massMatrices.o: OperatorsAndPotentials.cpp legendrebasis.o
	g++ -I /usr/local/include/eigen -L /usr/lib legendrebasis.cpp -c OperatorsAndPotentials.cpp

errors.o: errors.cpp geometry.o quadTables.o legendrebasis.o
	g++ -I /usr/local/include/eigen -L /usr/lib geometry.cpp quadTables.cpp legendrebasis.cpp -c errors.cpp 

testPotentials.o: geometry.o legendrebasis.o quadTables.o matrixRoutines.o
	g++ -I /usr/local/include/eigen -L /usr/lib geometry.cpp quadTables.cpp legendrebasis.cpp matrixRoutines.cpp -c testPotentials.cpp 



# scripts down here

testGeometry: geometry.o
	g++ -I /usr/local/include/eigen geometry.cpp scriptTestGeometry.cpp -o scriptTestGeometry

testQuadTables: quadTables.cpp scriptTestQuadrature.cpp
	g++ -I /usr/local/include/eigen quadTables.cpp scriptTestQuadrature.cpp -o scriptTestQuadrature

testLegendre: scriptTestLegendre.cpp
	g++ -L /usr/lib scriptTestLegendre.cpp -o scriptTestLegendre

scriptTestLegendreBasis: legendrebasis.o scriptTestLegendreBasis.cpp
	g++ -I /usr/local/include/eigen -L /usr/lib legendrebasis.cpp scriptTestLegendreBasis.cpp -o scriptTestLegendreBasis

scriptTestMatrixRoutines: matrixRoutines.o scriptTestMatrixRoutines.cpp
	g++ -I /usr/local/include/eigen -L /usr/lib matrixRoutines.cpp scriptTestMatrixRoutines.cpp -o scriptTestMatrixRoutines

scriptTestTestsAndProjections: testsAndProjections.o quadTables.o matrixRoutines.o legendrebasis.o geometry.o massMatrices.o scriptTestTestsAndProjections.cpp
	 g++ -I /usr/local/include/eigen -L /usr/lib matrixRoutines.cpp quadTables.cpp legendrebasis.cpp geometry.cpp testsAndProjections.cpp OperatorsAndPotentials.cpp scriptTestTestsAndProjections.cpp  -o scriptTestTestsAndProjections

scriptTestMassMatrices: massMatrices.o geometry.o quadTables.o legendrebasis.o matrixRoutines.o
	g++ -I /usr/local/include/eigen -L /usr/lib  matrixRoutines.cpp quadTables.cpp legendrebasis.cpp geometry.cpp OperatorsAndPotentials.cpp scriptTestMassMatrices.cpp -o scriptTestMassMatrices

scriptTestErrors: errors.o geometry.o quadTables.o legendrebasis.o testsAndProjections.o matrixRoutines.o massMatrices.o scriptTestErrors.cpp
	g++ -I /usr/local/include/eigen -L /usr/lib errors.cpp geometry.cpp quadTables.cpp legendrebasis.cpp testsAndProjections.cpp matrixRoutines.cpp OperatorsAndPotentials.cpp scriptTestErrors.cpp -o scriptTestErrors

scriptTestDiffMatrix: geometry.o massMatrices.o matrixRoutines.o legendrebasis.o
	g++ -I /usr/local/include/eigen -L /usr/lib geometry.cpp OperatorsAndPotentials.cpp matrixRoutines.cpp legendrebasis.cpp scriptTestDiffMatrix.cpp -o scriptTestDiffMatrix 

scriptTestPotentials: geometry.o legendrebasis.o quadTables.o testPotentials.o matrixRoutines.o
	g++ -I /usr/local/include/eigen -L /usr/lib geometry.cpp legendrebasis.cpp quadTables.cpp testPotentials.cpp matrixRoutines.cpp scriptTesttestPotentials.cpp -o scriptTesttestPotentials

scriptReshape: scriptReshape.cpp
	g++ -I /usr/local/include/eigen -L /usr/lib scriptReshape.cpp -o scriptReshape

clean:
	rm *.o
