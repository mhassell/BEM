all: main

main: geometry.o quadTables.o main.cpp
	g++ -I /usr/local/include/eigen geometry.cpp quadTables.cpp main.cpp -o main

geometry.o: geometry.cpp
	g++ -I /usr/local/include/eigen -c geometry.cpp

quadTables.o: quadTables.cpp
	g++ -I /usr/local/include/eigen -c quadTables.cpp

clean:
	rm *.o
