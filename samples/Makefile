all: linalg-test

linalg-test: linalg-test.cpp ../linalg.h
	$(CXX) linalg-test.cpp -std=c++11 -o $@

clean:
	rm -f linalg-test