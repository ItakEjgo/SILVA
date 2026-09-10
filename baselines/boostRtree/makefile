all: main

main:		main.cpp
	g++ -O3 -DNDEBUG -std=c++20 -Wall -I../../include/parlay -I../../include -I../../src -o main main.cpp -pthread

clean:
	rm -f test main


# g++ -O3 -DNDEBUG -std=c++17 -Wall -Iinclude -o main main.cpp --L/data/bhuan102/libspatialindex/lib64 lspatialindex
# g++ -g -O3 -DDEBUG -std=c++17 -Wall -Iinclude -o main main.cpp -L/data/bhuan102/libspatialindex/lib64 -lspatialindex