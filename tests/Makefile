sources = test-linalg.cpp test-expression-validity.cpp test-higher-order-functions.cpp test-constant-expressions.cpp test-quaternions.cpp 

all: test-linalg

test-linalg: $(sources) ../linalg.h
	$(CXX) $(sources) -std=c++11 $(CXXLIBRARY) -o $@

clean:
	rm -f test-linalg