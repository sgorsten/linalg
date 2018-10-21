sources:=test-linalg.cpp 
sources+=test-constant-expressions.cpp 
sources+=test-expression-validity.cpp
sources+=test-higher-order-functions.cpp
sources+=test-matrices.cpp 
sources+=test-ostream-overloads.cpp
sources+=test-quaternions.cpp 
sources+=test-scalar-accessors.cpp 
sources+=test-user-defined-conversions.cpp

all: time tests

time:
	time --format="  Including linalg.h in an otherwise empty file took %es" $(CXX) include-linalg.cpp -std=c++11 $(CXXLIBRARY) -c

tests: $(sources) ../linalg.h
	$(CXX) $(sources) -std=c++11 $(CXXLIBRARY) -o $@

generate-swizzles: generate-swizzles.cpp
	$(CXX) $< -std=c++11 $(CXXLIBRARY) -o $@
	./generate-swizzles
	
clean:
	rm -f include-linalg.o tests generate-swizzles