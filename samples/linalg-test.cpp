#include "../linalg.h"
using namespace linalg::aliases;

#define CATCH_CONFIG_MAIN
#include "thirdparty/catch.hpp"

#include <random>
#include <algorithm>

template<class T, int M> void require_zero(const linalg::vec<T,M> & v) { for(int i=0; i<M; ++i) REQUIRE( v[i] == 0 ); }
template<class T, int M, int N> void require_zero(const linalg::mat<T,M,N> & m) { for(int j=0; j<N; ++j) require_zero(m[j]); }

TEST_CASE( "linalg types default construct to zero" ) 
{
    require_zero( uint2() );
    require_zero( int3x4() );
    require_zero( ushort4() );
    require_zero( float4x2() );
    require_zero( double2x3() );
}

TEST_CASE( "equality and inequality operators" )
{
    // Two vectors compare equal if all elements are equal
    REQUIRE(   float4(1,2,3,4) == float4(1,2,3,4)  );
    REQUIRE( !(float4(1,2,3,4) == float4(5,2,3,4)) );
    REQUIRE( !(float4(1,2,3,4) == float4(1,5,3,4)) );
    REQUIRE( !(float4(1,2,3,4) == float4(1,2,5,4)) );
    REQUIRE( !(float4(1,2,3,4) == float4(1,2,3,5)) );

    // Two vectors compare inequal if at least one element is inequal
    REQUIRE( !(float4(1,2,3,4) != float4(1,2,3,4)) );
    REQUIRE(   float4(1,2,3,4) != float4(5,2,3,4)  );
    REQUIRE(   float4(1,2,3,4) != float4(1,5,3,4)  );
    REQUIRE(   float4(1,2,3,4) != float4(1,2,5,4)  );
    REQUIRE(   float4(1,2,3,4) != float4(1,2,3,5)  );
}

TEST_CASE( "vector and matrix can be constructed from pointer to const elements" )
{
    // Vectors load elements in order
    const float f2[] = {1,2}, f3[] = {2,3,4}, f4[] = {5,6,7,8};
    REQUIRE( float2(f2) == float2(1,2) );
    REQUIRE( float3(f3) == float3(2,3,4) );
    REQUIRE( float4(f4) == float4(5,6,7,8) );

    // Matrices load elements in column major order
    const float f2x2[] = {1,2,3,4}, f2x3[] = {2,3,4,5,6,7}, f2x4[] = {3,4,5,6,7,8,9,10};
    REQUIRE( float2x2(f2x2) == float2x2({1,2},{3,4}) );
    REQUIRE( float2x3(f2x3) == float2x3({2,3},{4,5},{6,7}) );
    REQUIRE( float2x4(f2x4) == float2x4({3,4},{5,6},{7,8},{9,10}) );

    // Should be possible to load the same data in multiple ways
    REQUIRE( float3x2(f2x3) == float3x2({2,3,4},{5,6,7}) );
    REQUIRE( float4x2(f2x4) == float4x2({3,4,5,6},{7,8,9,10}) );
}

TEST_CASE( "operator overloads produce correct results" )
{
    // All operators can be applied to vector types, and results are computed elementwise
    REQUIRE( (float2(2,3) +  float2(4,12)) == float2(  6,  15) );
    REQUIRE( (float2(2,3) -  float2(4,12)) == float2( -2,  -9) );
    REQUIRE( (float2(2,3) *  float2(4,12)) == float2(  8,  36) );
    REQUIRE( (float2(2,3) /  float2(4,12)) == float2(0.5,0.25) );
    REQUIRE( (int2(27,31) %  int2(5,8))    == int2( 2,  7) );
    REQUIRE( (int2(27,31) |  int2(5,8))    == int2(31, 31) );
    REQUIRE( (int2(27,31) ^  int2(5,8))    == int2(30, 23) );
    REQUIRE( (int2(27,31) &  int2(5,8))    == int2( 1,  8) );
    REQUIRE( (int2(14,35) << int2(2,3))    == int2(56,280) );
    REQUIRE( (int2(14,35) >> int2(2,3))    == int2( 3,  4) );
}

TEST_CASE( "integer promotion rules apply" )
{
    // uint8_t promotes to int
    REQUIRE(+byte2()              == int2());
    REQUIRE(-byte3()              == int3());
    REQUIRE(~byte4()              == int4(~0));
    REQUIRE(!byte2()              == bool2(true)); // operator! always returns a vector of bools
    REQUIRE((byte2() +  byte2())  == int2());
    REQUIRE((byte3() -  byte3())  == int3());
    REQUIRE((byte4() *  byte4())  == int4());
    REQUIRE((byte2() /  byte2(1)) == int2());
    REQUIRE((byte3() %  byte3(1)) == int3());
    REQUIRE((byte4() |  byte4())  == int4());
    REQUIRE((byte2() ^  byte2())  == int2());
    REQUIRE((byte3() &  byte3())  == int3());
    REQUIRE((byte4() << byte4())  == int4());
    REQUIRE((byte2() >> byte2())  == int2());

    // int16_t promotes to int
    REQUIRE(+short3()               == int3());
    REQUIRE(-short4()               == int4());
    REQUIRE(~short2()               == int2(~0));
    REQUIRE(!short3()               == bool3(true)); // operator! always returns a vector of bools
    REQUIRE((short3() +  short3())  == int3());
    REQUIRE((short4() -  short4())  == int4());
    REQUIRE((short2() *  short2())  == int2());
    REQUIRE((short3() /  short3(1)) == int3());
    REQUIRE((short4() %  short4(1)) == int4());
    REQUIRE((short2() |  short2())  == int2());
    REQUIRE((short3() ^  short3())  == int3());
    REQUIRE((short4() &  short4())  == int4());
    REQUIRE((short2() << short2())  == int2());
    REQUIRE((short3() >> short3())  == int3());

    // uint16_t promotes to int
    REQUIRE(+ushort4()                == int4());
    REQUIRE(-ushort2()                == int2());
    REQUIRE(~ushort3()                == int3(~0));
    REQUIRE(!ushort4()                == bool4(true)); // operator! always returns a vector of bools
    REQUIRE((ushort4() +  ushort4())  == int4());
    REQUIRE((ushort2() -  ushort2())  == int2());
    REQUIRE((ushort3() *  ushort3())  == int3());
    REQUIRE((ushort4() /  ushort4(1)) == int4());
    REQUIRE((ushort2() %  ushort2(1)) == int2());
    REQUIRE((ushort3() |  ushort3())  == int3());
    REQUIRE((ushort4() ^  ushort4())  == int4());
    REQUIRE((ushort2() &  ushort2())  == int2());
    REQUIRE((ushort3() << ushort3())  == int3());
    REQUIRE((ushort4() >> ushort4())  == int4());

    // int is not promoted
    REQUIRE(+int2()             == int2());
    REQUIRE(-int3()             == int3());
    REQUIRE(~int4()             == int4(~0));
    REQUIRE(!int2()             == bool2(true)); // operator! always returns a vector of bools
    REQUIRE((int2() +  int2())  == int2());
    REQUIRE((int3() -  int3())  == int3());
    REQUIRE((int4() *  int4())  == int4());
    REQUIRE((int2() /  int2(1)) == int2());
    REQUIRE((int3() %  int3(1)) == int3());
    REQUIRE((int4() |  int4())  == int4());
    REQUIRE((int2() ^  int2())  == int2());
    REQUIRE((int3() &  int3())  == int3());
    REQUIRE((int4() << int4())  == int4());
    REQUIRE((int2() >> int2())  == int2());

    // unsigned is not promoted
    REQUIRE(+uint3()              == uint3());
    REQUIRE(-uint4()              == uint4()); // NOTE: Will produce a warning about unary minus applied to unsigned type, this is probably desired behavior
    REQUIRE(~uint2()              == uint2(~0));
    REQUIRE(!uint3()              == bool3(true)); // operator! always returns a vector of bools
    REQUIRE((uint3() +  uint3())  == uint3());
    REQUIRE((uint4() -  uint4())  == uint4());
    REQUIRE((uint2() *  uint2())  == uint2());
    REQUIRE((uint3() /  uint3(1)) == uint3());
    REQUIRE((uint4() %  uint4(1)) == uint4());
    REQUIRE((uint2() |  uint2())  == uint2());
    REQUIRE((uint3() ^  uint3())  == uint3());
    REQUIRE((uint4() &  uint4())  == uint4());
    REQUIRE((uint2() << uint2())  == uint2());
    REQUIRE((uint3() >> uint3())  == uint3());

    // float is not promoted
    REQUIRE(+float4()              == float4());
    REQUIRE(-float2()              == float2());
    REQUIRE(!float4()              == bool4(true)); // operator! always returns a vector of bools
    REQUIRE((float4() + float4()) == float4());
    REQUIRE((float2() - float2()) == float2());
    REQUIRE((float3() * float3(1)) == float3());
    REQUIRE((float4() / float4(1)) == float4());
}

TEST_CASE( "elementwise comparison functions produce correct results" )
{
    REQUIRE( equal  (float3(1,2,3), float3(4,-2,3)) == bool3(false, false, true ) );
    REQUIRE( nequal (float3(1,2,3), float3(4,-2,3)) == bool3(true,  true,  false) );
    REQUIRE( less   (float3(1,2,3), float3(4,-2,3)) == bool3(true,  false, false) );
    REQUIRE( greater(float3(1,2,3), float3(4,-2,3)) == bool3(false, true,  false) );
    REQUIRE( lequal (float3(1,2,3), float3(4,-2,3)) == bool3(true,  false, true ) );
    REQUIRE( gequal (float3(1,2,3), float3(4,-2,3)) == bool3(false, true,  true ) );
}

TEST_CASE( "no unintended ADL on operator +=" )
{
    std::vector<int3> tris_a = {{0,1,2}, {0,2,3}, {0,3,4}}, tris_b;
    tris_b = std::move(tris_a); // This line is known to cause problems if linalg::operator+= is allowed to match too broadly.
    REQUIRE( tris_b.size() == 3 );
    REQUIRE( tris_a.size() == 0 );
}

TEST_CASE( "unary functions behave as intended" )
{
    // Unary functions should apply elementwise to their arguments
    REQUIRE( abs  (float3(1.1f, -2.3f, 3.5f)) == float3(std::abs  (1.1f), std::abs  (-2.3f), std::abs  (3.5f)) );
    REQUIRE( floor(float3(1.1f, -2.3f, 3.5f)) == float3(std::floor(1.1f), std::floor(-2.3f), std::floor(3.5f)) );
    REQUIRE( ceil (float3(1.1f, -2.3f, 3.5f)) == float3(std::ceil (1.1f), std::ceil (-2.3f), std::ceil (3.5f)) );
    REQUIRE( exp  (float3(1.1f, -2.3f, 3.5f)) == float3(std::exp  (1.1f), std::exp  (-2.3f), std::exp  (3.5f)) );
    REQUIRE( log  (float3(1.1f, +2.3f, 3.5f)) == float3(std::log  (1.1f), std::log  (+2.3f), std::log  (3.5f)) );
    REQUIRE( log10(float3(1.1f, +2.3f, 3.5f)) == float3(std::log10(1.1f), std::log10(+2.3f), std::log10(3.5f)) );
    REQUIRE( sqrt (float3(1.1f, +2.3f, 3.5f)) == float3(std::sqrt (1.1f), std::sqrt (+2.3f), std::sqrt (3.5f)) );
    REQUIRE( sin  (float3(1.1f, -2.3f, 3.5f)) == float3(std::sin  (1.1f), std::sin  (-2.3f), std::sin  (3.5f)) );
    REQUIRE( cos  (float3(1.1f, -2.3f, 3.5f)) == float3(std::cos  (1.1f), std::cos  (-2.3f), std::cos  (3.5f)) );
    REQUIRE( tan  (float3(1.1f, -2.3f, 3.5f)) == float3(std::tan  (1.1f), std::tan  (-2.3f), std::tan  (3.5f)) );
    REQUIRE( asin (float3(0.3f, -0.6f, 0.9f)) == float3(std::asin (0.3f), std::asin (-0.6f), std::asin (0.9f)) );
    REQUIRE( acos (float3(0.3f, -0.6f, 0.9f)) == float3(std::acos (0.3f), std::acos (-0.6f), std::acos (0.9f)) );
    REQUIRE( atan (float3(1.1f, -2.3f, 3.5f)) == float3(std::atan (1.1f), std::atan (-2.3f), std::atan (3.5f)) );
    REQUIRE( sinh (float3(1.1f, -2.3f, 3.5f)) == float3(std::sinh (1.1f), std::sinh (-2.3f), std::sinh (3.5f)) );
    REQUIRE( cosh (float3(1.1f, -2.3f, 3.5f)) == float3(std::cosh (1.1f), std::cosh (-2.3f), std::cosh (3.5f)) );
    REQUIRE( tanh (float3(1.1f, -2.3f, 3.5f)) == float3(std::tanh (1.1f), std::tanh (-2.3f), std::tanh (3.5f)) );
    REQUIRE( round(float3(1.1f, -2.3f, 3.5f)) == float3(std::round(1.1f), std::round(-2.3f), std::round(3.5f)) );

    // Unary functions should retain element type and vector/matrix size
    REQUIRE( abs(int4(-5)) == int4(5) );
    REQUIRE( floor(float2(7.7f)) == float2(std::floor(7.7f)) );
    REQUIRE( ceil (float3(7.7f)) == float3(std::ceil (7.7f)) );
    REQUIRE( exp  (float4(7.7f)) == float4(std::exp  (7.7f)) );
    REQUIRE( log  (double2(7.7)) == double2(std::log  (7.7)) );
    REQUIRE( log10(double3(7.7)) == double3(std::log10(7.7)) );
    REQUIRE( sqrt (double4(7.7)) == double4(std::sqrt (7.7)) );
    REQUIRE( sin  (float2x2(7.7f)) == float2x2(std::sin (7.7f)) );
    REQUIRE( cos  (float2x3(7.7f)) == float2x3(std::cos (7.7f)) );
    REQUIRE( tan  (float2x4(7.7f)) == float2x4(std::tan (7.7f)) );
    REQUIRE( asin (float3x2(0.5f)) == float3x2(std::asin(0.5f)) );
    REQUIRE( acos (float3x3(0.5f)) == float3x3(std::acos(0.5f)) );
    REQUIRE( atan (float3x4(7.7f)) == float3x4(std::atan(7.7f)) );
    REQUIRE( sinh (float4x2(7.7f)) == float4x2(std::sinh(7.7f)) );
    REQUIRE( cosh (float4x3(7.7f)) == float4x3(std::cosh(7.7f)) );
    REQUIRE( tanh (float4x4(7.7f)) == float4x4(std::tanh(7.7f)) );
}

TEST_CASE( "relational operators model LessThanComparable" )
{
    // Should not compare equal to begin with
    std::vector<int3> points = {{3,5,2}, {1,2,6}, {3,2,2}, {8,2,5}, {4,5,8}, {3,5,8}, {1,2,2}, {4,2,9}};
    std::vector<int3> ordered_points = {{1,2,2}, {1,2,6}, {3,2,2}, {3,5,2}, {3,5,8}, {4,2,9}, {4,5,8}, {8,2,5}};
    REQUIRE( points != ordered_points );

    // After sorting, points should be in ascending order, and match ordered_points
    std::sort(begin(points), end(points));
    REQUIRE( points == ordered_points );

    // Sorting in descending order should produce the reverse ordering
    std::sort(begin(points), end(points), [](const int3 & a, const int3 & b) { return a > b; });
    std::reverse(begin(ordered_points), end(ordered_points));
    REQUIRE( points == ordered_points );
}

TEST_CASE( "matrix multiplication produces correct result dimensions" )
{
    REQUIRE( mul(float2x2(), float2()) == float2() );
    REQUIRE( mul(float2x3(), float3()) == float2() );
    REQUIRE( mul(float2x4(), float4()) == float2() );
    REQUIRE( mul(float3x2(), float2()) == float3() );
    REQUIRE( mul(float3x3(), float3()) == float3() );
    REQUIRE( mul(float3x4(), float4()) == float3() );
    REQUIRE( mul(float4x2(), float2()) == float4() );
    REQUIRE( mul(float4x3(), float3()) == float4() );
    REQUIRE( mul(float4x4(), float4()) == float4() );

    REQUIRE( mul(float2x2(), float2x2()) == float2x2() );
    REQUIRE( mul(float2x3(), float3x2()) == float2x2() );
    REQUIRE( mul(float2x4(), float4x2()) == float2x2() );
    REQUIRE( mul(float2x2(), float2x3()) == float2x3() );
    REQUIRE( mul(float2x3(), float3x3()) == float2x3() );
    REQUIRE( mul(float2x4(), float4x3()) == float2x3() );
    REQUIRE( mul(float2x2(), float2x4()) == float2x4() );
    REQUIRE( mul(float2x3(), float3x4()) == float2x4() );
    REQUIRE( mul(float2x4(), float4x4()) == float2x4() );
    REQUIRE( mul(float3x2(), float2x2()) == float3x2() );
    REQUIRE( mul(float3x3(), float3x2()) == float3x2() );
    REQUIRE( mul(float3x4(), float4x2()) == float3x2() );
    REQUIRE( mul(float3x2(), float2x3()) == float3x3() );
    REQUIRE( mul(float3x3(), float3x3()) == float3x3() );
    REQUIRE( mul(float3x4(), float4x3()) == float3x3() );
    REQUIRE( mul(float3x2(), float2x4()) == float3x4() );
    REQUIRE( mul(float3x3(), float3x4()) == float3x4() );
    REQUIRE( mul(float3x4(), float4x4()) == float3x4() );
    REQUIRE( mul(float4x2(), float2x2()) == float4x2() );
    REQUIRE( mul(float4x3(), float3x2()) == float4x2() );
    REQUIRE( mul(float4x4(), float4x2()) == float4x2() );
    REQUIRE( mul(float4x2(), float2x3()) == float4x3() );
    REQUIRE( mul(float4x3(), float3x3()) == float4x3() );
    REQUIRE( mul(float4x4(), float4x3()) == float4x3() );
    REQUIRE( mul(float4x2(), float2x4()) == float4x4() );
    REQUIRE( mul(float4x3(), float3x4()) == float4x4() );
    REQUIRE( mul(float4x4(), float4x4()) == float4x4() );

    // Outer product of vec<T,M> and vec<T,N> is equivalent to product of Mx1 and 1xN matrices
    REQUIRE( outerprod(float2(), float2()) == float2x2() );
    REQUIRE( outerprod(float2(), float3()) == float2x3() );
    REQUIRE( outerprod(float2(), float4()) == float2x4() );
    REQUIRE( outerprod(float3(), float2()) == float3x2() );
    REQUIRE( outerprod(float3(), float3()) == float3x3() );
    REQUIRE( outerprod(float3(), float4()) == float3x4() );
    REQUIRE( outerprod(float4(), float2()) == float4x2() );
    REQUIRE( outerprod(float4(), float3()) == float4x3() );
    REQUIRE( outerprod(float4(), float4()) == float4x4() );
}

TEST_CASE( "matrix inverse is correct for trivial cases" )
{
    const float2x2 id2 {{1,0},{0,1}};
    const float3x3 id3 {{1,0,0},{0,1,0},{0,0,1}};
    const float4x4 id4 {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
    REQUIRE( diagonal(id2) == float2(1,1) );
    REQUIRE( diagonal(id3) == float3(1,1,1) );
    REQUIRE( diagonal(id4) == float4(1,1,1,1) );
    REQUIRE( transpose(id2) == id2 );
    REQUIRE( transpose(id3) == id3 );
    REQUIRE( transpose(id4) == id4 );
    REQUIRE( inverse(id2) == id2 );
    REQUIRE( inverse(id3) == id3 );
    REQUIRE( inverse(id4) == id4 );
    REQUIRE( adjugate(id2) == id2 );
    REQUIRE( adjugate(id3) == id3 );
    REQUIRE( adjugate(id4) == id4 );
    REQUIRE( determinant(id2) == 1.0f );
    REQUIRE( determinant(id3) == 1.0f );
    REQUIRE( determinant(id4) == 1.0f );
}

TEST_CASE( "matrix inverse is correct for general case" )
{
    const float4x4 mat {{1,2,3,4}, {5,-6,7,8}, {9,10,-11,12}, {13,14,15,-16}};
    const float4x4 inv = inverse(mat);
    const float4x4 id = mul(mat, inv);
    for(int j=0; j<4; ++j)
    {
        for(int i=0; i<4; ++i)
        {
            if(i == j) REQUIRE( id[j][i] == Approx(1.0f) );
            else REQUIRE( id[j][i] == Approx(0.0f) );
        }
    }
}

TEST_CASE( "rotation quaternions roundtrip with rotation matrices" )
{
    std::mt19937 engine;
    std::normal_distribution<float> f;
    std::normal_distribution<double> d;

    for(int i=0; i<1000; ++i)
    {
        float4 q = normalize(float4(f(engine), f(engine), f(engine), f(engine)));
        if(q.w < 0) q = -q; // rotation_quat(float3x3) always returns a quat with q.w >= 0
        float4 q2 = rotation_quat(qmat(q));

        REQUIRE( abs(q.x - q2.x) < 0.0001f );
        REQUIRE( abs(q.y - q2.y) < 0.0001f );
        REQUIRE( abs(q.z - q2.z) < 0.0001f );
        REQUIRE( abs(q.w - q2.w) < 0.0001f );
    }

    for(int i=0; i<1000; ++i)
    {
        double4 q = normalize(double4(d(engine), d(engine), d(engine), d(engine)));
        if(q.w < 0) q = -q; // rotation_quat(double3x3) always returns a quat with q.w >= 0
        double4 q2 = rotation_quat(qmat(q));

        REQUIRE( abs(q.x - q2.x) < 0.00000001 );
        REQUIRE( abs(q.y - q2.y) < 0.00000001 );
        REQUIRE( abs(q.z - q2.z) < 0.00000001 );
        REQUIRE( abs(q.w - q2.w) < 0.00000001 );
    }
}

TEST_CASE( "hashing works as expected" )
{
    // std::hash specializations should take their specified type and return size_t
    REQUIRE( typeid(std::hash<int2     >()(int2     {})) == typeid(size_t) );
    REQUIRE( typeid(std::hash<float3   >()(float3   {})) == typeid(size_t) );
    REQUIRE( typeid(std::hash<double4  >()(double4  {})) == typeid(size_t) );
    REQUIRE( typeid(std::hash<int2x4   >()(int2x4   {})) == typeid(size_t) );
    REQUIRE( typeid(std::hash<float3x2 >()(float3x2 {})) == typeid(size_t) );
    REQUIRE( typeid(std::hash<double4x3>()(double4x3{})) == typeid(size_t) );

    // Small list of items which are known to have no duplicate hashes
    const float3 points[] = {
        {1,2,3}, {1,3,2}, {2,1,3}, {2,3,1}, {3,1,2}, {3,2,1},
        {1,1,1}, {2,2,2}, {3,3,3}, {4,4,4}, {5,5,5}, {6,6,6},
        {0,0,0}, {0,0,1}, {0,1,0}, {1,0,0}, {0,1,1}, {1,1,0}
    };
    const std::hash<float3> h = {};
    for(auto & a : points)
    {
        for(auto & b : points)
        {
            if(a == b) REQUIRE( h(a) == h(b) );
            else REQUIRE( h(a) != h(b) ); // Note: Not required to be different in the general case
        }
    }
}

template<class T> void take(const T &) {}
#define MATCH(TYPE, ...) static_assert(std::is_same<TYPE, decltype(__VA_ARGS__)>::value, #TYPE " != " #__VA_ARGS__); take(__VA_ARGS__)

TEST_CASE( "templates instantiate correctly", "" ) 
{
    // Declare some variables to test functions requiring an lvalue
    const float2 cf2; const float3 cf3; const float4 cf4;
    float2 f2; float3 f3; float4 f4; int2 i2; int3 i3; int4 i4;
    float fs[] = {1,2,3,4};

    // Exercise vec<T,2>
    MATCH(float2, float2());
    MATCH(float2, float2(1,2) );
    MATCH(float2, float2(5) );
    MATCH(float2, float2(fs) );
    MATCH(float2, float2(int2(3,4)) );
    MATCH(const float&, cf2[1] );
    MATCH(float&, f2[1] );

    // Exercise vec<T,3>
    MATCH(float3, float3() );
    MATCH(float3, float3(1,2,3) );
    MATCH(float3, float3(float2(),4) );
    MATCH(float3, float3(5) );
    MATCH(float3, float3(fs) );
    MATCH(float3, float3(int3(3,4,5)) );
    MATCH(const float&, cf3[1] );
    MATCH(float&, f3[1] );
    MATCH(float2 &, f3.xy() );

    // Exercise vec<T,4>
    MATCH(float4, float4() );
    MATCH(float4, float4(1,2,3,4) );
    MATCH(float4, float4(float3(),4) );
    MATCH(float4, float4(5) );
    MATCH(float4, float4(fs) );
    MATCH(float4, float4(int4(3,4,5,6)) );
    MATCH(const float&, cf4[1] );
    MATCH(float&, f4[1] );
    MATCH(float3 &, f4.xyz() );

    // TODO: Exercise mat<T,M,N> for N=2,3,4

    // Exercise sequence functions
    for(float & f : f4) take(f);
    for(float f : float4()) take(f);
    for(float4 & f : float4x4()) take(f);

    // Exercise relational operators
    MATCH(bool, int2() == int2() );
    MATCH(bool, float3() == float3() );
    MATCH(bool, double4() == double4() );
    MATCH(bool, short2() != short2() );
    MATCH(bool, int2() < int2() );
    MATCH(bool, float3() < float3() );
    MATCH(bool, double4() < double4() );
    MATCH(bool, int2() > int2() );
    MATCH(bool, float3() <= float3() );
    MATCH(bool, double4() >= double4() );

    // Exercise binary operators
    MATCH(float2 , float2 () +  float2 () );
    MATCH(float3 , float3 () -  float3 () );
    MATCH(double4, double4() *  double4() );
    MATCH(double2, double2() /  double2() );
    MATCH(int3   , int3   () %  int3   (1) );
    MATCH(int4   , int4   () |  int4   () );
    MATCH(int2   , short2 () ^  short2 () );
    MATCH(int3   , short3 () &  short3 () );
    MATCH(int3   , int3   () << int3   () );
    MATCH(int4   , int4   () >> int4   () );

    MATCH(float2 , float2 () +  float   () );
    MATCH(float3 , float3 () -  float   () );
    MATCH(double4, double4() *  double  () );
    MATCH(double2, double2() /  double  () );
    MATCH(int3   , int3   () %  int     (1) );
    MATCH(int4   , int4   () |  int     () );
    MATCH(int2   , short2 () ^  short   () );
    MATCH(int3   , short3 () &  short   () );
    MATCH(int3   , int3   () << int     () );
    MATCH(uint4  , uint4  () >> unsigned() );

    MATCH(float2 , float   () +  float2 () );
    MATCH(float3 , float   () -  float3 () );
    MATCH(double4, double  () *  double4() );
    MATCH(double2, double  () /  double2() );
    MATCH(int3   , int     () %  int3   (1) );
    MATCH(int4   , int     () |  int4   () );
    MATCH(int2   , short   () ^  short2 () );
    MATCH(int3   , short   () &  short3 () );
    MATCH(int3   , int     () << int3   () );
    MATCH(uint4  , unsigned() >> uint4  () );

    MATCH(float2&, f2 += float2() );
    MATCH(float3&, f3 -= float3() );
    MATCH(float4&, f4 *= float4() );
    MATCH(float2&, f2 /= float2() );
    MATCH(int2&  , i2 %= int2(1) );
    MATCH(int3&  , i3 |= int3() );
    MATCH(int4&  , i4 ^= int4() );
    MATCH(int2&  , i2 &= int2() );
    MATCH(int3&  , i3 <<= int3() );
    MATCH(int4&  , i4 >>= int4() );

    MATCH(float2&, f2 += float() );
    MATCH(float3&, f3 -= float() );
    MATCH(float4&, f4 *= float() );
    MATCH(float2&, f2 /= float() );
    MATCH(int2&  , i2 %= int(1) );
    MATCH(int3&  , i3 |= int() );
    MATCH(int4&  , i4 ^= int() );
    MATCH(int2&  , i2 &= int() );
    MATCH(int3&  , i3 <<= int() );
    MATCH(int4&  , i4 >>= int() );

    MATCH(float2 , min  (float2(), float2()) );
    MATCH(float3 , max  (float3(), float3()) );
    MATCH(float2 , fmod (float2(), float2()) );
    MATCH(float3 , pow  (float3(), float3()) );
    MATCH(float4 , atan2(float4(), float4()) );
    MATCH(float3 , clamp(float3(), float3(), float3()) );

    // Exercise componentwise relational ops
    MATCH(bool2, equal  (float2(), float2()) );
    MATCH(bool3, nequal (float3(), float3()) );
    MATCH(bool4, less   (double4(), double4()) );
    MATCH(bool2, greater(double2(), double2()) );
    MATCH(bool3, lequal (int3(), int3()) );
    MATCH(bool4, gequal (int4(), int4()) );

    // Exercise reduction functions
    MATCH(bool, any(bool3()) );
    MATCH(bool, all(bool4()) );
    MATCH(int, sum(int2()) );
    MATCH(float, product(float4()) );

    // Exercise selection functions
    MATCH(int, argmin(float2()) );
    MATCH(int, argmax(float3()) );
    MATCH(float, minelem(float4()) );
    MATCH(float, maxelem(float2()) );

    // Exercise vector algebra functions
    MATCH(float, cross(float2(), float2()) );
    MATCH(float3, cross(float3(), float3()) );
    MATCH(float, dot(float4(), float4()) );
    MATCH(float, length2(float2()) );
    MATCH(float, length(float3()) );
    MATCH(float, distance2(float4(), float4()) );
    MATCH(float, distance(float2(), float2()) );
    MATCH(float3, normalize(float3()) );
    MATCH(float4, lerp(float4(), float4(), float()) );
    MATCH(float4, nlerp(float4(), float4(), float()) );
    MATCH(float4, slerp(float4(), float4(), float()) );

    // Exercise quaternion algebra functions
    MATCH(float4, qconj(float4()) );
    MATCH(float4, qinv(float4()) );
    MATCH(float4, qmul(float4(), float4()) );
    MATCH(float3, qxdir(float4()) );
    MATCH(float3, qydir(float4()) );
    MATCH(float3, qzdir(float4()) );
    MATCH(float3, qrot(float4(), float3()) );
    MATCH(float , qangle(float4()) );
    MATCH(float3, qaxis(float4()) );
    MATCH(float4, qnlerp(float4(), float4(), float()) );
    MATCH(float4, qslerp(float4(), float4(), float()) );

    // Exercise factory functions
    MATCH(float4, rotation_quat(float3(), float()) );
    MATCH(float4x4, translation_matrix(float3()) );
    MATCH(float4x4, rotation_matrix(float4()) );
    MATCH(float4x4, scaling_matrix(float3()) );
    MATCH(float4x4, pose_matrix(float4(), float3()) );
    MATCH(float4x4, linalg::frustum_matrix(float(), float(), float(), float(), float(), float()) );
    MATCH(float4x4, linalg::perspective_matrix(float(), float(), float(), float()) );
}