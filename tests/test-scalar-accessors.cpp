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

TEST_CASE("Test swizzle reads")
{
    const float2 a {1,2};
    CHECK(a.xy == float2{1,2});
    CHECK(a.yx == float2{2,1});

    const float3 b {3,4,5};
    CHECK(b.xy == float2{3,4}); 
    CHECK(b.xz == float2{3,5});
    CHECK(b.yx == float2{4,3});
    CHECK(b.yz == float2{4,5}); 
    CHECK(b.zx == float2{5,3}); 
    CHECK(b.zy == float2{5,4});
    CHECK(b.xyz == float3{3,4,5});
    CHECK(b.xzy == float3{3,5,4});
    CHECK(b.yxz == float3{4,3,5});
    CHECK(b.yzx == float3{4,5,3});
    CHECK(b.zxy == float3{5,3,4});
    CHECK(b.zyx == float3{5,4,3});

    const float4 c {6,7,8,9};
    CHECK(c.xy == float2{6,7});
    CHECK(c.xz == float2{6,8});
    CHECK(c.xw == float2{6,9});
    CHECK(c.yx == float2{7,6});
    CHECK(c.yz == float2{7,8});
    CHECK(c.yw == float2{7,9});
    CHECK(c.zx == float2{8,6});
    CHECK(c.zy == float2{8,7});
    CHECK(c.zw == float2{8,9});
    CHECK(c.wx == float2{9,6});
    CHECK(c.wy == float2{9,7});
    CHECK(c.wz == float2{9,8});
    CHECK(c.xyz == float3{6,7,8});
    CHECK(c.xyw == float3{6,7,9});
    CHECK(c.xzy == float3{6,8,7});
    CHECK(c.xzw == float3{6,8,9});
    CHECK(c.xwy == float3{6,9,7});
    CHECK(c.xwz == float3{6,9,8});
    CHECK(c.yxz == float3{7,6,8});
    CHECK(c.yxw == float3{7,6,9});
    CHECK(c.yzx == float3{7,8,6});
    CHECK(c.yzw == float3{7,8,9});
    CHECK(c.ywx == float3{7,9,6});
    CHECK(c.ywz == float3{7,9,8});
    CHECK(c.zxy == float3{8,6,7});
    CHECK(c.zxw == float3{8,6,9});
    CHECK(c.zyx == float3{8,7,6});
    CHECK(c.zyw == float3{8,7,9});
    CHECK(c.zwx == float3{8,9,6});
    CHECK(c.zwy == float3{8,9,7});
    CHECK(c.wxy == float3{9,6,7});
    CHECK(c.wxz == float3{9,6,8});
    CHECK(c.wyx == float3{9,7,6});
    CHECK(c.wyz == float3{9,7,8});
    CHECK(c.wzx == float3{9,8,6});
    CHECK(c.wzy == float3{9,8,7});
    CHECK(c.xyzw == float4{6,7,8,9});
    CHECK(c.xywz == float4{6,7,9,8});
    CHECK(c.xzyw == float4{6,8,7,9});
    CHECK(c.xzwy == float4{6,8,9,7});
    CHECK(c.xwyz == float4{6,9,7,8});
    CHECK(c.xwzy == float4{6,9,8,7});
    CHECK(c.yxzw == float4{7,6,8,9});
    CHECK(c.yxwz == float4{7,6,9,8});
    CHECK(c.yzxw == float4{7,8,6,9});
    CHECK(c.yzwx == float4{7,8,9,6});
    CHECK(c.ywxz == float4{7,9,6,8});
    CHECK(c.ywzx == float4{7,9,8,6});
    CHECK(c.zxyw == float4{8,6,7,9});
    CHECK(c.zxwy == float4{8,6,9,7});
    CHECK(c.zyxw == float4{8,7,6,9});
    CHECK(c.zywx == float4{8,7,9,6});
    CHECK(c.zwxy == float4{8,9,6,7});
    CHECK(c.zwyx == float4{8,9,7,6});
    CHECK(c.wxyz == float4{9,6,7,8});
    CHECK(c.wxzy == float4{9,6,8,7});
    CHECK(c.wyxz == float4{9,7,6,8});
    CHECK(c.wyzx == float4{9,7,8,6});
    CHECK(c.wzxy == float4{9,8,6,7});
    CHECK(c.wzyx == float4{9,8,7,6});
}

TEST_CASE("Test swizzle writes")
{
    float4 a {1,2,3,4}, b {5,6,7,8};

    SUBCASE("Can assign through swizzle with explicit type")
    {
        a.zy = float2{9,10};
        CHECK(a[0] == 1);
        CHECK(a[1] == 10);
        CHECK(a[2] == 9);
        CHECK(a[3] == 4);
    }

    SUBCASE("Can assign through swizzle with braced list type")
    {
        a.zy = {9,10};
        CHECK(a[0] == 1);
        CHECK(a[1] == 10);
        CHECK(a[2] == 9);
        CHECK(a[3] == 4);
    }

    SUBCASE("Can assign from same swizzle")
    {
        a.wz = b.wz;
        CHECK(a[0] == 1);
        CHECK(a[1] == 2);
        CHECK(a[2] == 7);
        CHECK(a[3] == 8); 
    }

    SUBCASE("Can assign from different swizzle")
    {
        a.wz = b.xy;
        CHECK(a[0] == 1);
        CHECK(a[1] == 2);
        CHECK(a[2] == 6);
        CHECK(a[3] == 5);    
    }
}