#pragma once
#include "../linalg.h"
#include "thirdparty/doctest.h"
#include <random>
using namespace linalg::aliases;

// Tests which use random data will be repeated this many times
constexpr int reps = 5;

// Template aliases to make it easier to use linalg types in templated test cases
template<class T> using vec2 = linalg::vec<T,2>;
template<class T> using vec3 = linalg::vec<T,3>;
template<class T> using vec4 = linalg::vec<T,4>;

// Facility for retrieving random numbers
class random_number_generator
{
    std::mt19937 rng;
    std::normal_distribution<double> dist_double;
    std::normal_distribution<float> dist_float;
    std::uniform_int_distribution<int> dist_int;
    std::uniform_int_distribution<short> dist_short;
    std::uniform_int_distribution<unsigned> dist_uint;
    std::uniform_int_distribution<unsigned short> dist_ushort;
public:
    random_number_generator() : dist_int(-1000, 1000), dist_short(-100, 100), dist_uint(0, 1000), dist_ushort(0, 100) {}
    
    template<class T> T get() { return get((T*)nullptr); }

    double get(double *) { return dist_double(rng); }
    float get(float *) { return dist_float(rng); }
    int get(int *) { return dist_int(rng); }
    short get(short *) { return dist_short(rng); }
    unsigned int get(unsigned int *) { return dist_uint(rng); }
    unsigned short get(unsigned short *) { return dist_ushort(rng); }

    template<class T> linalg::vec<T,2> get(linalg::vec<T,2> *) { return linalg::vec<T,2>(get<T>(), get<T>()); }
    template<class T> linalg::vec<T,3> get(linalg::vec<T,3> *) { return linalg::vec<T,3>(get<T>(), get<T>(), get<T>()); }
    template<class T> linalg::vec<T,4> get(linalg::vec<T,4> *) { return linalg::vec<T,4>(get<T>(), get<T>(), get<T>(), get<T>()); }
    template<class T, int M> linalg::mat<T,M,2> get(linalg::mat<T,M,2> *) { return linalg::mat<T,M,2>(get<linalg::vec<T,M>>(), get<linalg::vec<T,M>>()); }
    template<class T, int M> linalg::mat<T,M,3> get(linalg::mat<T,M,3> *) { return linalg::mat<T,M,3>(get<linalg::vec<T,M>>(), get<linalg::vec<T,M>>(), get<linalg::vec<T,M>>()); }
    template<class T, int M> linalg::mat<T,M,4> get(linalg::mat<T,M,4> *) { return linalg::mat<T,M,4>(get<linalg::vec<T,M>>(), get<linalg::vec<T,M>>(), get<linalg::vec<T,M>>(), get<linalg::vec<T,M>>()); }

    template<class T> operator T () { return get<T>(); }
};

template<class T, int M> void check_approx_equal(const linalg::vec<T,M> & a, const linalg::vec<T,M> & b) { for(int j=0; j<M; ++j) CHECK( a[j] == doctest::Approx(b[j]) ); }
template<class T, int M, int N> void check_approx_equal(const linalg::mat<T,M,N> & a, const linalg::mat<T,M,N> & b) { for(int i=0; i<N; ++i) check_approx_equal(a[i], b[i]); }
