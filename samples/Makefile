all: linalg-test

linalg-test: linalg-test.cpp ../linalg.h
	$(CXX) linalg-test.cpp linalg-test2.cpp linalg-constant-expressions.cpp linalg-expression-validity.cpp -std=c++14 $(CXXLIBRARY) -o $@

clean:
	rm -f linalg-test