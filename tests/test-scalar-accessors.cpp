#include "test-linalg.h"

TEST_CASE("Scalar accessors behave as intended")
{
    int2 i2 {1,2};
    int3 i3 {3,4,5};
    int4 i4 {6,7,8,9};

    const int2 ci2 {1,2};
    const int3 ci3 {3,4,5};
    const int4 ci4 {6,7,8,9};

    // Check address-of operator for non-const
    CHECK(&i2.x == &i2[0]);
    CHECK(&i2.y == &i2[1]);
    CHECK(&i3.x == &i3[0]);
    CHECK(&i3.y == &i3[1]);
    CHECK(&i3.z == &i3[2]);
    CHECK(&i4.x == &i4[0]);
    CHECK(&i4.y == &i4[1]);
    CHECK(&i4.z == &i4[2]);
    CHECK(&i4.w == &i4[3]);

    // Check address-of operator for const
    CHECK(&ci2.x == &ci2[0]);
    CHECK(&ci2.y == &ci2[1]);
    CHECK(&ci3.x == &ci3[0]);
    CHECK(&ci3.y == &ci3[1]);
    CHECK(&ci3.z == &ci3[2]);
    CHECK(&ci4.x == &ci4[0]);
    CHECK(&ci4.y == &ci4[1]);
    CHECK(&ci4.z == &ci4[2]);
    CHECK(&ci4.w == &ci4[3]);

    // Check xyzw accessors for non-const
    CHECK(i2.x == 1);
    CHECK(i2.y == 2);
    CHECK(i3.x == 3);
    CHECK(i3.y == 4);    
    CHECK(i3.z == 5);    
    CHECK(i4.x == 6);
    CHECK(i4.y == 7);    
    CHECK(i4.z == 8);
    CHECK(i4.w == 9);

    // Check xyzw accessors for const
    CHECK(ci2.x == 1);
    CHECK(ci2.y == 2);
    CHECK(ci3.x == 3);
    CHECK(ci3.y == 4);    
    CHECK(ci3.z == 5);    
    CHECK(ci4.x == 6);
    CHECK(ci4.y == 7);    
    CHECK(ci4.z == 8);
    CHECK(ci4.w == 9);

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

    // Accessors should allow assignment
    i4.x = 10;
    i4.y = i4.z;
    CHECK(i4[0] == 10);
    CHECK(i4[1] == i4[2]);
}

TEST_CASE("Scalar accessors assign as scalars, not as vectors")
{
    // Assigning via braced initializer list should assign one element, should not reassign entire vector
    float4 a {1,2,3,4}, b {5,6,7,8};

    SUBCASE("Assignment of number assigns as a scalar")
    {
        a.x = 9;
        CHECK(a[0] == 9);
        CHECK(a[1] == 2);
        CHECK(a[2] == 3);
        CHECK(a[3] == 4);    
    }

    SUBCASE("Assignment of braced number assigns as a scalar")
    {
        a.y = {9};
        CHECK(a[0] == 1);
        CHECK(a[1] == 9);
        CHECK(a[2] == 3);
        CHECK(a[3] == 4);    
    }

    SUBCASE("Assignment from same accessor")
    {
        a.z = b.z;
        CHECK(a[0] == 1);
        CHECK(a[1] == 2);
        CHECK(a[2] == 7);
        CHECK(a[3] == 4);    
    }

    SUBCASE("Assignment from different accessor")
    {
        a.w = b.y;
        CHECK(a[0] == 1);
        CHECK(a[1] == 2);
        CHECK(a[2] == 3);
        CHECK(a[3] == 6);    
    }
}