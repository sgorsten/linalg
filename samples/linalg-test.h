#pragma once
#include "../linalg.h"
#include "thirdparty/doctest.h"
#include <random>
using namespace linalg::aliases;

// Type lists to use in templated test cases
using floating_point_types = doctest::Types<double, float>;
using integral_types = doctest::Types<int, short, unsigned int, unsigned short>;
using signed_types = doctest::Types<double, float, int, short>;
using arithmetic_types = doctest::Types<double, float, int, short, unsigned int, unsigned short>;

// SFINAE based expression validity helpers
namespace detail
{
    // These function overloads exist if specific the expression T{} op U{} is valid, and return true_type
    template<class T, class U> inline decltype(std::declval<T>() + std::declval<U>(), std::true_type{}) has_op_add(int) { return {}; }
    template<class T, class U> inline decltype(std::declval<T>() - std::declval<U>(), std::true_type{}) has_op_sub(int) { return {}; }  
    template<class T, class U> inline decltype(std::declval<T>() * std::declval<U>(), std::true_type{}) has_op_mul(int) { return {}; }
    template<class T, class U> inline decltype(std::declval<T>() / std::declval<U>(), std::true_type{}) has_op_div(int) { return {}; }
    
    // These function overloads always exist, but have lowest selection priority, and return false_type
    template<class T, class U> inline std::false_type has_op_add(...) { return {}; }
    template<class T, class U> inline std::false_type has_op_sub(...) { return {}; } 
    template<class T, class U> inline std::false_type has_op_mul(...) { return {}; }
    template<class T, class U> inline std::false_type has_op_div(...) { return {}; }
}
template<class T, class U> constexpr bool has_op_add = decltype(detail::has_op_add<T,U>(0))::value;
template<class T, class U> constexpr bool has_op_sub = decltype(detail::has_op_sub<T,U>(0))::value;
template<class T, class U> constexpr bool has_op_mul = decltype(detail::has_op_mul<T,U>(0))::value;
template<class T, class U> constexpr bool has_op_div = decltype(detail::has_op_div<T,U>(0))::value;

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

    operator double () { return dist_double(rng); }
    operator float () { return dist_float(rng); }
    operator int () { return dist_int(rng); }
    operator short () { return dist_short(rng); }
    operator unsigned int () { return dist_uint(rng); }
    operator unsigned short () { return dist_ushort(rng); }
    template<class T> operator linalg::vec<T,2> () { return linalg::vec<T,2>((T)*this, (T)*this); }
    template<class T> operator linalg::vec<T,3> () { return linalg::vec<T,3>((T)*this, (T)*this, (T)*this); }
    template<class T> operator linalg::vec<T,4> () { return linalg::vec<T,4>((T)*this, (T)*this, (T)*this, *this); }
    template<class T, int M> operator linalg::mat<T,M,2> () { return linalg::mat<T,M,2>((linalg::vec<T,M>)*this, (linalg::vec<T,M>)*this); }
    template<class T, int M> operator linalg::mat<T,M,3> () { return linalg::mat<T,M,3>((linalg::vec<T,M>)*this, (linalg::vec<T,M>)*this, (linalg::vec<T,M>)*this); }
    template<class T, int M> operator linalg::mat<T,M,4> () { return linalg::mat<T,M,4>((linalg::vec<T,M>)*this, (linalg::vec<T,M>)*this, (linalg::vec<T,M>)*this, (linalg::vec<T,M>)*this); }
};