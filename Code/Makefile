objects = test.o

all: test

test:$(objects)
	g++ -O3 $(objects) -o test 
test.o:test.cpp
	g++ -std=c++2a -c -O3 test.cpp
.PHONY: clean
clean:
	-rm test *.o