all: linalg-test

linalg-test:
	$(CXX) linalg-test.cpp linalg-test2.cpp test-user-defined-conversions.cpp test-matrices.cpp test-ostream-overloads.cpp test-constant-expressions.cpp test-higher-order-functions.cpp test-quaternions.cpp -std=c++11 -o $@

clean:
	rm -f linalg-test