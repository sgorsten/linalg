#include "test-linalg.h"
#include <array>

namespace linalg
{
    // Emulate construction of linalg types from std::array
    template<class T> struct converter<vec<T,2>, std::array<T,2>> { vec<T,2> operator() (const std::array<T,2> & a) const { return {a[0], a[1]}; } };
    template<class T> struct converter<vec<T,3>, std::array<T,3>> { vec<T,3> operator() (const std::array<T,3> & a) const { return {a[0], a[1], a[2]}; } };
    template<class T> struct converter<vec<T,4>, std::array<T,4>> { vec<T,4> operator() (const std::array<T,4> & a) const { return {a[0], a[1], a[2], a[3]}; } };
}

TEST_CASE("Test user-defined conversion from std::array<T,M> to linalg::vec<T,M>")
{
    float2 a = std::array<float,2>{1,2};
    float3 b = std::array<float,3>{3,4,5};
    float4 c = std::array<float,4>{6,7,8,9};

    CHECK(a == float2{1,2});
    CHECK(b == float3{3,4,5});
    CHECK(c == float4{6,7,8,9});
}

namespace linalg
{
    // Emulate construction of vectors from pointers to const elements
    template<class T> struct converter<vec<T,2>, const T *> { vec<T,2> operator() (const T * p) const { return {p[0], p[1]}; } };
    template<class T> struct converter<vec<T,3>, const T *> { vec<T,3> operator() (const T * p) const { return {p[0], p[1], p[2]}; } };
    template<class T> struct converter<vec<T,4>, const T *> { vec<T,4> operator() (const T * p) const { return {p[0], p[1], p[2], p[3]}; } };

    // Emulate construction of matrices from pointers to const elements
    template<class T, int M> struct converter<mat<T,M,2>, const T *> { converter<vec<T,M>, const T *> c; mat<T,M,2> operator() (const T * p) const { return {c(p+M*0), c(p+M*1)}; } };
    template<class T, int M> struct converter<mat<T,M,3>, const T *> { converter<vec<T,M>, const T *> c; mat<T,M,3> operator() (const T * p) const { return {c(p+M*0), c(p+M*1), c(p+M*2)}; } };
    template<class T, int M> struct converter<mat<T,M,4>, const T *> { converter<vec<T,M>, const T *> c; mat<T,M,4> operator() (const T * p) const { return {c(p+M*0), c(p+M*1), c(p+M*2), c(p+M*3)}; } };

    // Enable construction of linalg types from pointers to non-const elements
    template<class T, int M> struct converter<vec<T,M>, T *> : converter<vec<T,M>, const T *> {};
    template<class T, int M, int N> struct converter<mat<T,M,N>, T *> : converter<mat<T,M,N>, const T *> {};
}

TEST_CASE("Test construction of linalg types from pointers to elements, using user-defined conversions")
{
    float elements[] {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
    float * ptr_m = elements;
    const float * ptr_c = elements;

    CHECK(float2(ptr_m) == float2{1,2});
    CHECK(float3(ptr_m) == float3{1,2,3});
    CHECK(float4(ptr_m) == float4{1,2,3,4});
    CHECK(float2x2(ptr_m) == float2x2{{1,2},{3,4}});
    CHECK(float3x3(ptr_m) == float3x3{{1,2,3},{4,5,6},{7,8,9}});
    CHECK(float4x4(ptr_m) == float4x4{{1,2,3,4},{5,6,7,8},{9,10,11,12},{13,14,15,16}});

    CHECK(float2(ptr_c) == float2{1,2});
    CHECK(float3(ptr_c) == float3{1,2,3});
    CHECK(float4(ptr_c) == float4{1,2,3,4});
    CHECK(float2x2(ptr_c) == float2x2{{1,2},{3,4}});
    CHECK(float3x3(ptr_c) == float3x3{{1,2,3},{4,5,6},{7,8,9}});
    CHECK(float4x4(ptr_c) == float4x4{{1,2,3,4},{5,6,7,8},{9,10,11,12},{13,14,15,16}});
}

// Define some "external" types
struct Vec2 { float x,y; };
struct Vec3 { float x,y,z; };
struct Vec4 { float x,y,z,w; };
struct Mat2 { float m[4]; };
struct Mat3 { float m[9]; };
struct Mat4 { float m[16]; };
struct Quat { float w,x,y,z; };

namespace linalg
{
    // Define conversion from linalg types to "external" types
    template<> struct converter<Vec2, float2> { constexpr Vec2 operator() (const float2 & v) const { return {v[0], v[1]}; } };
    template<> struct converter<Vec3, float3> { constexpr Vec3 operator() (const float3 & v) const { return {v[0], v[1], v[2]}; } };
    template<> struct converter<Vec4, float4> { constexpr Vec4 operator() (const float4 & v) const { return {v[0], v[1], v[2], v[3]}; } };
    template<> struct converter<Mat2, float2x2> { constexpr Mat2 operator() (const float2x2 & m) const { return {m[0][0], m[0][1], m[1][0], m[1][1]}; } };
    template<> struct converter<Mat3, float3x3> { constexpr Mat3 operator() (const float3x3 & m) const { return {m[0][0], m[0][1], m[0][2], m[1][0], m[1][1], m[1][2], m[2][0], m[2][1], m[2][2]}; } };
    template<> struct converter<Mat4, float4x4> { constexpr Mat4 operator() (const float4x4 & m) const { return {m[0][0], m[0][1], m[0][2], m[0][3], m[1][0], m[1][1], m[1][2], m[1][3], m[2][0], m[2][1], m[2][2], m[2][3], m[3][0], m[3][1], m[3][2], m[3][3]}; } };
    template<> struct converter<Quat, quatf> { constexpr Quat operator() (const quatf & q) const { return {q.w, q.x, q.y, q.z}; } };

    // Define conversion from "external" types to linalg types
    template<> struct converter<float2, Vec2> { constexpr float2 operator() (const Vec2 & v) const { return {v.x, v.y}; } };
    template<> struct converter<float3, Vec3> { constexpr float3 operator() (const Vec3 & v) const { return {v.x, v.y, v.z}; } };
    template<> struct converter<float4, Vec4> { constexpr float4 operator() (const Vec4 & v) const { return {v.x, v.y, v.z, v.w}; } };
    template<> struct converter<float2x2, Mat2> { constexpr float2x2 operator() (const Mat2 & m) const { return {{m.m[0], m.m[1]}, {m.m[2], m.m[3]}}; } };
    template<> struct converter<float3x3, Mat3> { constexpr float3x3 operator() (const Mat3 & m) const { return {{m.m[0], m.m[1], m.m[2]}, {m.m[3], m.m[4], m.m[5]}, {m.m[6], m.m[7], m.m[8]}}; } };
    template<> struct converter<float4x4, Mat4> { constexpr float4x4 operator() (const Mat4 & m) const { return {{m.m[0], m.m[1], m.m[2], m.m[3]}, {m.m[4], m.m[5], m.m[6], m.m[7]}, {m.m[8], m.m[9], m.m[10], m.m[11]}, {m.m[12], m.m[13], m.m[14], m.m[15]}}; } };
    template<> struct converter<quatf, Quat> { constexpr quatf operator() (const Quat & q) const { return {q.x, q.y, q.z, q.w}; } };
}

TEST_CASE("Test user-defined conversions for linalg::vec<T,2>")
{
    float2 a {1,2};
    Vec2 b = a;
    CHECK(b.x == 1);
    CHECK(b.y == 2);

    float2 c = b;
    CHECK(c == a);
}

TEST_CASE("Test user-defined conversions for linalg::vec<T,3>")
{
    float3 a {1,2,3};
    Vec3 b = a;
    CHECK(b.x == 1);
    CHECK(b.y == 2);
    CHECK(b.z == 3);

    float3 c = b;
    CHECK(c == a);
}

TEST_CASE("Test user-defined conversions for linalg::vec<T,4>")
{
    float4 a {1,2,3,4};
    Vec4 b = a;
    CHECK(b.x == 1);
    CHECK(b.y == 2);
    CHECK(b.z == 3);
    CHECK(b.w == 4);

    float4 c = b;
    CHECK(c == a);
}

TEST_CASE("Test user-defined conversions for linalg::mat<T,M,2>")
{
    float2x2 a {{1,2},{3,4}};
    Mat2 b = a;
    CHECK(b.m[0] == 1);
    CHECK(b.m[1] == 2);
    CHECK(b.m[2] == 3);
    CHECK(b.m[3] == 4);

    float2x2 c = b;
    CHECK(c == a);
}

TEST_CASE("Test user-defined conversions for linalg::mat<T,M,3>")
{
    float3x3 a {{1,2,3},{4,5,6},{7,8,9}};
    Mat3 b = a;
    CHECK(b.m[0] == 1);
    CHECK(b.m[1] == 2);
    CHECK(b.m[2] == 3);
    CHECK(b.m[3] == 4);
    CHECK(b.m[4] == 5);
    CHECK(b.m[5] == 6);
    CHECK(b.m[6] == 7);
    CHECK(b.m[7] == 8);
    CHECK(b.m[8] == 9);

    float3x3 c = b;
    CHECK(c == a);
}

TEST_CASE("Test user-defined conversions for linalg::mat<T,M,4>")
{
    float4x4 a {{1,2,3,4},{5,6,7,8},{9,10,11,12},{13,14,15,16}};
    Mat4 b = a;
    CHECK(b.m[0] == 1);
    CHECK(b.m[1] == 2);
    CHECK(b.m[2] == 3);
    CHECK(b.m[3] == 4);
    CHECK(b.m[4] == 5);
    CHECK(b.m[5] == 6);
    CHECK(b.m[6] == 7);
    CHECK(b.m[7] == 8);
    CHECK(b.m[8] == 9);
    CHECK(b.m[9] == 10);
    CHECK(b.m[10] == 11);
    CHECK(b.m[11] == 12);
    CHECK(b.m[12] == 13);
    CHECK(b.m[13] == 14);
    CHECK(b.m[14] == 15);
    CHECK(b.m[15] == 16);

    float4x4 c = b;
    CHECK(c == a);
}

TEST_CASE("Test user-defined conversions for linalg::quat<T>")
{
    quatf a {5,6,7,8};
    Quat b = a;
    CHECK(b.w == 8);
    CHECK(b.x == 5);
    CHECK(b.y == 6);
    CHECK(b.z == 7);    

    quatf c = b;
    CHECK(c == a);
}