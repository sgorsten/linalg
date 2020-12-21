#include "test-linalg.h"

// linalg::identity is a constant of type linalg::identity_t declared within a header file
// For C++11 compliance, we've declared it `static`, meaning that each translation unit technically has its own copy of linalg::identity.
// If we could assume C++17 support, I'd prefer to instead declare it `inline`. Either way, we need to doublecheck that no duplicate
// definition issues occur when we include linalg.h from multiple translation units.
TEST_CASE( "declared constants can be used from multiple translation units" )
{
    REQUIRE( sum(int2x2(linalg::identity)) == 2 );
    REQUIRE( sum(float3x3(linalg::identity)) == 3.0f );
    REQUIRE( sum(double4x4(linalg::identity)) == 4.0 );
}