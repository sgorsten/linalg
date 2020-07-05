#include "test-linalg.h"

// Test constexpr construction
namespace
{
    constexpr int1 a1 {1}, b1 {2};
    constexpr int2 a2 {1,2}, b2 {2,3};
    constexpr int3 a3 {1,2,4}, b3 {2,3,5};
    constexpr int4 a4 {1,2,4,8}, b4 {2,3,5,7};
    
    constexpr int2x2 a2x2 {{1,2},{2,4}}, b2x2 {{2,3},{3,5}};
    constexpr int3x3 a3x3 {{1,2,4},{2,4,8},{4,8,16}}, b3x3 {{2,3,5},{3,5,7},{5,7,11}};
    constexpr int4x4 a4x4 {{1,2,4,8},{2,4,8,16},{4,8,16,32},{8,16,32,64}}, b4x4 {{2,3,5,7},{3,5,7,11},{5,7,11,13},{7,11,13,17}};

    constexpr int4x4 m {{1,2,3,4},{5,6,7,8},{9,10,11,12},{13,14,15,16}};
}


// Check vec::operator[]
static_assert(a1[0] == 1, "linalg::vec<T,1>::operator[] should be constexpr");
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

// Check swizzle<...>
static_assert(linalg::swizzle<3,1,2,0>(b4) == int4{7,3,5,2}, "swizzle<...> can reorder");
static_assert(linalg::swizzle<0,1,1,0>(a2) == int4{1,2,2,1}, "swizzle<...> can make larger vectors");
static_assert(linalg::swizzle<1,2>(a4) == int2{2,4}, "swizzle<...> can make smaller vectors");

// Check mat::operator==, !=, <, >, <=, >=
static_assert(a2x2 == a2x2, "linalg::mat<T,M,N>::operator== should be constexpr");
static_assert(a3x3 != b3x3, "linalg::mat<T,M,N>::operator!= should be constexpr");
static_assert(a4x4 <  b4x4, "linalg::mat<T,M,N>::operator< should be constexpr");
static_assert(b2x2 >  a2x2, "linalg::mat<T,M,N>::operator> should be constexpr");
static_assert(a3x3 <= b3x3, "linalg::mat<T,M,N>::operator<= should be constexpr");
static_assert(b4x4 >= a4x4, "linalg::mat<T,M,N>::operator>= should be constexpr");

// Check minelem
static_assert(minelem(int4(-1, 2, 3, 4)) == -1, "minelem should be constexpr");
static_assert(minelem(int4(1, -2, 3, 4)) == -2, "minelem should be constexpr");
static_assert(minelem(int4(1, 2, -3, 4)) == -3, "minelem should be constexpr");
static_assert(minelem(int4(1, 2, 3, -4)) == -4, "minelem should be constexpr");

// Check maxelem
static_assert(maxelem(int4(1, -2, -3, -4)) == 1, "maxelem should be constexpr");
static_assert(maxelem(int4(-1, 2, -3, -4)) == 2, "maxelem should be constexpr");
static_assert(maxelem(int4(-1, -2, 3, -4)) == 3, "maxelem should be constexpr");
static_assert(maxelem(int4(-1, -2, -3, 4)) == 4, "maxelem should be constexpr");

/*
// Check member .data()
static_assert(a2.data() == (const int *)a2._._, "linalg::vec<T,2>::data() should be constexpr");
static_assert(a3.data() == (const int *)a3._._, "linalg::vec<T,3>::data() should be constexpr");
static_assert(a4.data() == (const int *)a4._._, "linalg::vec<T,4>::data() should be constexpr");
static_assert(a2x2.data() == (const int *)a2x2._[0]._._, "linalg::mat<T,M,2>::data() should be constexpr");
static_assert(a3x3.data() == (const int *)a3x3._[0]._._, "linalg::mat<T,M,3>::data() should be constexpr");
static_assert(a4x4.data() == (const int *)a4x4._[0]._._, "linalg::mat<T,M,4>::data() should be constexpr");
*/

// Check submat/subvec
static_assert(linalg::submat<0,0, 4,4>(m) == int4x4{{1,2,3,4},{5,6,7,8},{9,10,11,12},{13,14,15,16}}, "...");
static_assert(linalg::submat<0,0, 3,3>(m) == int3x3{{1,2,3},{5,6,7},{9,10,11}}, "...");
static_assert(linalg::submat<0,0, 2,2>(m) == int2x2{{1,2},{5,6}}, "...");
    
static_assert(linalg::submat<0,0, 4,3>(m) == int4x3{{1,2,3,4},{5,6,7,8},{9,10,11,12}}, "...");
static_assert(linalg::submat<0,1, 4,4>(m) == int4x3{{5,6,7,8},{9,10,11,12},{13,14,15,16}}, "...");
    
static_assert(linalg::submat<0,0, 4,2>(m) == int4x2{{1,2,3,4},{5,6,7,8}}, "...");
static_assert(linalg::submat<0,1, 4,3>(m) == int4x2{{5,6,7,8},{9,10,11,12}}, "...");
static_assert(linalg::submat<0,2, 4,4>(m) == int4x2{{9,10,11,12},{13,14,15,16}}, "...");
    
static_assert(linalg::submat<0,0, 3,4>(m) == int3x4{{1,2,3},{5,6,7},{9,10,11},{13,14,15}}, "...");
static_assert(linalg::submat<1,0, 4,4>(m) == int3x4{{2,3,4},{6,7,8},{10,11,12},{14,15,16}}, "...");

static_assert(linalg::submat<0,0, 2,4>(m) == int2x4{{1,2},{5,6},{9,10},{13,14}}, "...");
static_assert(linalg::submat<1,0, 3,4>(m) == int2x4{{2,3},{6,7},{10,11},{14,15}}, "...");
static_assert(linalg::submat<2,0, 4,4>(m) == int2x4{{3,4},{7,8},{11,12},{15,16}}, "...");

/*
// Check constexpr iterators
#if !defined(_MSC_VER) || _MSC_VER > 1900
static_assert(linalg::end(a2) - linalg::begin(a2) == 2, "begin()/end() should be constexpr");
static_assert(linalg::end(a3) - linalg::begin(a3) == 3, "begin()/end() should be constexpr");
static_assert(linalg::end(a4) - linalg::begin(a4) == 4, "begin()/end() should be constexpr");
#endif
*/