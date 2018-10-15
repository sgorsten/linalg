sources = test-linalg.cpp test-constant-expressions.cpp test-expression-validity.cpp test-higher-order-functions.cpp test-quaternions.cpp test-scalar-accessors.cpp test-user-defined-conversions.cpp

all: test-linalg

test-linalg: $(sources) ../linalg.h
	$(CXX) $(sources) -std=c++11 $(CXXLIBRARY) -o $@

clean:
	rm -f test-linalg