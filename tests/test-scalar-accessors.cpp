#include "test-linalg.h"

TEST_CASE("Scalar accessors behave as intended")
{
    const int2 i2 {1,2};
    const int3 i3 {3,4,5};
    const int4 i4 {6,7,8,9};

    // Check xyzw accessors
    CHECK(i2.x == 1);
    CHECK(i2.y == 2);
    CHECK(i3.x == 3);
    CHECK(i3.y == 4);    
    CHECK(i3.z == 5);    
    CHECK(i4.x == 6);
    CHECK(i4.y == 7);    
    CHECK(i4.z == 8);
    CHECK(i4.w == 9);

    // Check rgba accessors
    CHECK(i2.r == 1);
    CHECK(i2.g == 2);
    CHECK(i3.r == 3);
    CHECK(i3.g == 4);
    CHECK(i3.b == 5);
    CHECK(i4.r == 6);
    CHECK(i4.g == 7);
    CHECK(i4.b == 8);
    CHECK(i4.a == 9);

    // Check stpq accessors
    CHECK(i2.s == 1);
    CHECK(i2.t == 2);
    CHECK(i3.s == 3);
    CHECK(i3.t == 4);
    CHECK(i3.p == 5);
    CHECK(i4.s == 6);
    CHECK(i4.t == 7);
    CHECK(i4.p == 8);
    CHECK(i4.q == 9);

    // Accessors should behave like scalars in linalg ops
    CHECK(i2.x + i3 == int3{4,5,6});
    CHECK(i4 * i3.y == int4{24,28,32,36});
}