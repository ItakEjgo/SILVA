all: main

main:		main.cpp
	g++ -O3 -DNDEBUG -std=c++20 -Wall -I../../include/parlay -I../../include -I../../src -o main main.cpp -pthread

clean:
	rm -f test main
