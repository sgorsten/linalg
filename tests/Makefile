sources := test-linalg.cpp 
sources += test-constant-expressions.cpp 
sources += test-expression-validity.cpp
sources += test-higher-order-functions.cpp
sources += test-ostream-overloads.cpp
sources += test-quaternions.cpp 
sources += test-scalar-accessors.cpp 
sources += test-user-defined-conversions.cpp

all: test-linalg

test-linalg: $(sources) ../linalg.h
	$(CXX) $(sources) -std=c++11 $(CXXLIBRARY) -o $@

clean:
	rm -f test-linalg