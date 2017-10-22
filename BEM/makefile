all: main

main: geometry.o quadTables.o main.cpp
	g++ -I /usr/local/include/eigen geometry.cpp quadTables.cpp main.cpp -o main

geometry.o: geometry.cpp
	g++ -I /usr/local/include/eigen -fopenmp -march=native -DNDEBUG -O3 -c geometry.cpp

quadTables.o: quadTables.cpp
	g++ -I /usr/local/include/eigen -fopenmp -march=native -DNDEBUG -O3 -c quadTables.cpp

legendrebasis.o: legendrebasis.cpp
	g++ -I /usr/local/include/eigen -L /usr/lib -fopenmp -march=native -DNDEBUG -O3 -c legendrebasis.cpp

matrixRoutines.o: matrixRoutines.cpp
	g++ -I /usr/local/include/eigen -L /usr/lib -fopenmp -march=native -DNDEBUG -O3 -c matrixRoutines.cpp

testsAndProjections.o: testsAndProjections.cpp
	g++ -I /usr/local/include/eigen -L /usr/lib -fopenmp -march=native -DNDEBUG -O3 -c testsAndProjections.cpp

massMatrices.o: OperatorsAndPotentials.cpp legendrebasis.o
	g++ -I /usr/local/include/eigen -L /usr/lib -fopenmp -march=native -DNDEBUG -O3 legendrebasis.cpp -c OperatorsAndPotentials.cpp

errors.o: errors.cpp geometry.o quadTables.o legendrebasis.o
	g++ -I /usr/local/include/eigen -L /usr/lib -fopenmp -march=native -DNDEBUG -O3 geometry.cpp quadTables.cpp legendrebasis.cpp -c errors.cpp 

testPotentials.o: geometry.o legendrebasis.o quadTables.o matrixRoutines.o
	g++ -I /usr/local/include/eigen -L /usr/lib -fopenmp -march=native -DNDEBUG -O3 geometry.cpp quadTables.cpp legendrebasis.cpp matrixRoutines.cpp -c testPotentials.cpp 

quadrature.o: matrixRoutines.o quadTables.o 
	g++ -I /usr/local/include/eigen -L /usr/lib -fopenmp -march=native -DNDEBUG -O3 matrixRoutines.cpp quadTables.cpp quadrature.cpp -c quadrature.cpp

Operators.o: matrixRoutines.o geometry.o legendrebasis.o
	g++ -I /usr/local/include/eigen -L /usr/lib -fopenmp -march=native -DNDEBUG -O3 matrixRoutines.cpp geometry.cpp legendrebasis.cpp Operators.cpp -c Operators.cpp 

# scripts down here

testGeometry: geometry.o
	g++ -I /usr/local/include/eigen -fopenmp -march=native -DNDEBUG -O3 geometry.cpp scriptTestGeometry.cpp -o scriptTestGeometry

testQuadTables: quadTables.cpp scriptTestQuadrature.cpp
	g++ -I /usr/local/include/eigen -fopenmp -march=native -DNDEBUG -O3 quadTables.cpp scriptTestQuadrature.cpp -o scriptTestQuadrature

testLegendre: scriptTestLegendre.cpp
	g++ -L /usr/lib -fopenmp -march=native -DNDEBUG -O3 scriptTestLegendre.cpp -o scriptTestLegendre

scriptTestLegendreBasis: legendrebasis.o scriptTestLegendreBasis.cpp
	g++ -I /usr/local/include/eigen -L /usr/lib -fopenmp -march=native -DNDEBUG -O3 legendrebasis.cpp scriptTestLegendreBasis.cpp -o scriptTestLegendreBasis

scriptTestMatrixRoutines: matrixRoutines.o scriptTestMatrixRoutines.cpp
	g++ -I /usr/local/include/eigen -L /usr/lib -fopenmp -march=native -DNDEBUG -O3 matrixRoutines.cpp scriptTestMatrixRoutines.cpp -o scriptTestMatrixRoutines

scriptTestTestsAndProjections: testsAndProjections.o quadTables.o matrixRoutines.o legendrebasis.o geometry.o massMatrices.o scriptTestTestsAndProjections.cpp
	 g++ -I /usr/local/include/eigen -L /usr/lib -fopenmp -march=native -DNDEBUG -O3 matrixRoutines.cpp quadTables.cpp legendrebasis.cpp geometry.cpp testsAndProjections.cpp OperatorsAndPotentials.cpp scriptTestTestsAndProjections.cpp  -o scriptTestTestsAndProjections

scriptTestMassMatrices: massMatrices.o geometry.o quadTables.o legendrebasis.o matrixRoutines.o
	g++ -I /usr/local/include/eigen -L /usr/lib -fopenmp -march=native -DNDEBUG -O3 matrixRoutines.cpp quadTables.cpp legendrebasis.cpp geometry.cpp OperatorsAndPotentials.cpp scriptTestMassMatrices.cpp -o scriptTestMassMatrices

scriptTestErrors: errors.o geometry.o quadTables.o legendrebasis.o testsAndProjections.o matrixRoutines.o massMatrices.o scriptTestErrors.cpp
	g++ -I /usr/local/include/eigen -L /usr/lib -fopenmp -march=native -DNDEBUG -O3 errors.cpp geometry.cpp quadTables.cpp legendrebasis.cpp testsAndProjections.cpp matrixRoutines.cpp OperatorsAndPotentials.cpp scriptTestErrors.cpp -o scriptTestErrors

scriptTestDiffMatrix: geometry.o massMatrices.o matrixRoutines.o legendrebasis.o
	g++ -I /usr/local/include/eigen -L /usr/lib -fopenmp -march=native -DNDEBUG -O3 geometry.cpp OperatorsAndPotentials.cpp matrixRoutines.cpp legendrebasis.cpp scriptTestDiffMatrix.cpp -o scriptTestDiffMatrix 

scriptTestPotentials: geometry.o legendrebasis.o quadTables.o testPotentials.o matrixRoutines.o
	g++ -I /usr/local/include/eigen -L /usr/lib -fopenmp -march=native -DNDEBUG -O3 geometry.cpp legendrebasis.cpp quadTables.cpp testPotentials.cpp matrixRoutines.cpp scriptTesttestPotentials.cpp -o scriptTesttestPotentials

scriptReshape: scriptReshape.cpp
	g++ -I /usr/local/include/eigen -L /usr/lib -fopenmp -march=native -DNDEBUG -O3 scriptReshape.cpp -o scriptReshape

scriptTestQuadrature: scriptTestQuadrature.cpp quadrature.cpp matrixRoutines.o quadTables.o
	g++ -I /usr/local/include/eigen -L /usr/lib -fopenmp -march=native -DNDEBUG -O3 quadrature.cpp matrixRoutines.cpp quadTables.cpp scriptTestQuadrature.cpp -o scriptTestQuadrature

scriptTestOperators: Operators.cpp matrixRoutines.o quadTables.o quadrature.o legendrebasis.o geometry.o OperatorsAndPotentials.cpp scriptTestOperators.cpp
	g++ -I /usr/local/include/eigen -L /usr/lib -fopenmp -march=native -DNDEBUG -O3 Operators.cpp matrixRoutines.cpp quadTables.cpp quadrature.cpp legendrebasis.cpp geometry.cpp OperatorsAndPotentials.cpp scriptTestOperators.cpp -o scriptTestOperators

scriptTestReshaping: matrixRoutines.o testsAndProjections.o geometry.o quadrature.o quadTables.o legendrebasis.o massMatrices.o scriptTestReshaping.cpp
	g++ -I /usr/local/include/eigen -L /usr/lib -fopenmp -march=native -DNDEBUG -O3 matrixRoutines.cpp testsAndProjections.cpp geometry.cpp quadrature.cpp quadTables.cpp legendrebasis.cpp OperatorsAndPotentials.cpp scriptTestReshaping.cpp -o scriptTestReshaping

scriptCalderonCalculusLaplace: geometry.o quadrature.o Operators.o massMatrices.o testsAndProjections.o testPotentials.o matrixRoutines.o quadTables.o legendrebasis.o
	g++ -I /usr/local/include/eigen -L /usr/lib -fopenmp -march=native -DNDEBUG -O3 -pg geometry.cpp quadrature.cpp Operators.cpp OperatorsAndPotentials.cpp testsAndProjections.cpp testPotentials.cpp matrixRoutines.cpp quadTables.cpp legendrebasis.cpp scriptCalderonCalculusLaplace.cpp -o scriptCalderonCalculusLaplace
	
scriptTestOpenMP: scriptTestOpenMP.cpp
	g++ -o scriptTestOpenMP -fopenmp -march=native -DNDEBUG -O3 scriptTestOpenMP.cpp  

clean:
	rm *.o
