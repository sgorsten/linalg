#include "../linalg.h"
using namespace linalg::aliases;

// For constexpr tests, we will use static_assert to check that the desired logic occurs at compile time
static constexpr int2 a2 {1,2}, b2 {2,3};
static constexpr int3 a3 {1,2,4}, b3 {2,3,5};
static constexpr int4 a4 {1,2,4,8}, b4 {2,3,5,7};

// Check operator[]
static_assert(a2[0] == 1, "linalg::vec<T,2>::operator[] should be constexpr");
static_assert(a2[1] == 2, "linalg::vec<T,2>::operator[] should be constexpr");
static_assert(a3[0] == 1, "linalg::vec<T,3>::operator[] should be constexpr");
static_assert(a3[1] == 2, "linalg::vec<T,3>::operator[] should be constexpr");
static_assert(a3[2] == 4, "linalg::vec<T,3>::operator[] should be constexpr");
static_assert(a4[0] == 1, "linalg::vec<T,4>::operator[] should be constexpr");
static_assert(a4[1] == 2, "linalg::vec<T,4>::operator[] should be constexpr");
static_assert(a4[2] == 4, "linalg::vec<T,4>::operator[] should be constexpr");
static_assert(a4[3] == 8, "linalg::vec<T,4>::operator[] should be constexpr");

// static_assert(a == a, "operator ==");

