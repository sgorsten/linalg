#include "../linalg.h"
using namespace linalg::aliases;

// Test constexpr construction
namespace
{
    constexpr int2 a2 {1,2}, b2 {2,3};
    constexpr int3 a3 {1,2,4}, b3 {2,3,5};
    constexpr int4 a4 {1,2,4,8}, b4 {2,3,5,7};
    
    constexpr int2x2 a2x2 {{1,2},{2,4}}, b2x2 {{2,3},{3,5}};
    constexpr int3x3 a3x3 {{1,2,4},{2,4,8},{4,8,16}}, b3x3 {{2,3,5},{3,5,7},{5,7,11}};
    constexpr int4x4 a4x4 {{1,2,4,8},{2,4,8,16},{4,8,16,32},{8,16,32,64}}, b4x4 {{2,3,5,7},{3,5,7,11},{5,7,11,13},{7,11,13,17}};
}

// Check vec::operator[]
static_assert(a2[0] == 1, "linalg::vec<T,2>::operator[] should be constexpr");
static_assert(a2[1] == 2, "linalg::vec<T,2>::operator[] should be constexpr");
static_assert(a3[0] == 1, "linalg::vec<T,3>::operator[] should be constexpr");
static_assert(a3[1] == 2, "linalg::vec<T,3>::operator[] should be constexpr");
static_assert(a3[2] == 4, "linalg::vec<T,3>::operator[] should be constexpr");
static_assert(a4[0] == 1, "linalg::vec<T,4>::operator[] should be constexpr");
static_assert(a4[1] == 2, "linalg::vec<T,4>::operator[] should be constexpr");
static_assert(a4[2] == 4, "linalg::vec<T,4>::operator[] should be constexpr");
static_assert(a4[3] == 8, "linalg::vec<T,4>::operator[] should be constexpr");

// Check vec::operator==, !=, <, >, <=, >=
static_assert(a2 == a2, "linalg::vec<T,M>::operator== should be constexpr");
static_assert(a3 != b3, "linalg::vec<T,M>::operator!= should be constexpr");
static_assert(a4 <  b4, "linalg::vec<T,M>::operator< should be constexpr");
static_assert(b2 >  a2, "linalg::vec<T,M>::operator> should be constexpr");
static_assert(a3 <= b3, "linalg::vec<T,M>::operator<= should be constexpr");
static_assert(b4 >= a4, "linalg::vec<T,M>::operator>= should be constexpr");

// Check mat::operator[]
static_assert(a2x2[0][0] == 1,  "linalg::mat<T,M,2>::operator[] should be constexpr");
static_assert(a2x2[1][1] == 4,  "linalg::mat<T,M,2>::operator[] should be constexpr");
static_assert(a3x3[0][0] == 1,  "linalg::mat<T,M,3>::operator[] should be constexpr");
static_assert(a3x3[1][1] == 4,  "linalg::mat<T,M,3>::operator[] should be constexpr");
static_assert(a3x3[2][2] == 16, "linalg::mat<T,M,3>::operator[] should be constexpr");
static_assert(a4x4[0][0] == 1,  "linalg::mat<T,M,4>::operator[] should be constexpr");
static_assert(a4x4[1][1] == 4,  "linalg::mat<T,M,4>::operator[] should be constexpr");
static_assert(a4x4[2][2] == 16, "linalg::mat<T,M,4>::operator[] should be constexpr");
static_assert(a4x4[3][3] == 64, "linalg::mat<T,M,4>::operator[] should be constexpr");

// Check mat::operator==, !=, <, >, <=, >=
static_assert(a2x2 == a2x2, "linalg::mat<T,M,N>::operator== should be constexpr");
static_assert(a3x3 != b3x3, "linalg::mat<T,M,N>::operator!= should be constexpr");
static_assert(a4x4 <  b4x4, "linalg::mat<T,M,N>::operator< should be constexpr");
static_assert(b2x2 >  a2x2, "linalg::mat<T,M,N>::operator> should be constexpr");
static_assert(a3x3 <= b3x3, "linalg::mat<T,M,N>::operator<= should be constexpr");
static_assert(b4x4 >= a4x4, "linalg::mat<T,M,N>::operator>= should be constexpr");
