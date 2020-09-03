#include "test-linalg.h"

TEST_CASE( "matrix transpose produces correct result dimensions" )
{
    REQUIRE( transpose(float1x1()) == float1x1() );
    REQUIRE( transpose(float1x2()) == float2x1() );
    REQUIRE( transpose(float1x3()) == float3x1() );
    REQUIRE( transpose(float1x4()) == float4x1() );
    REQUIRE( transpose(float2x1()) == float1x2() );
    REQUIRE( transpose(float2x2()) == float2x2() );
    REQUIRE( transpose(float2x3()) == float3x2() );
    REQUIRE( transpose(float2x4()) == float4x2() );
    REQUIRE( transpose(float3x1()) == float1x3() );
    REQUIRE( transpose(float3x2()) == float2x3() );
    REQUIRE( transpose(float3x3()) == float3x3() );
    REQUIRE( transpose(float3x4()) == float4x3() );
    REQUIRE( transpose(float4x1()) == float1x4() );
    REQUIRE( transpose(float4x2()) == float2x4() );
    REQUIRE( transpose(float4x3()) == float3x4() );
    REQUIRE( transpose(float4x4()) == float4x4() );

    REQUIRE( transpose(float1()) == float1x1() );
    REQUIRE( transpose(float2()) == float1x2() );
    REQUIRE( transpose(float3()) == float1x3() );
    REQUIRE( transpose(float4()) == float1x4() );
}

TEST_CASE( "matrix multiplication produces correct result dimensions" )
{
    REQUIRE( mul(float1x1(), float1()) == float1() );
    REQUIRE( mul(float1x2(), float2()) == float1() );
    REQUIRE( mul(float1x3(), float3()) == float1() );
    REQUIRE( mul(float1x4(), float4()) == float1() );
    REQUIRE( mul(float2x1(), float1()) == float2() );
    REQUIRE( mul(float2x2(), float2()) == float2() );
    REQUIRE( mul(float2x3(), float3()) == float2() );
    REQUIRE( mul(float2x4(), float4()) == float2() );
    REQUIRE( mul(float3x1(), float1()) == float3() );
    REQUIRE( mul(float3x2(), float2()) == float3() );
    REQUIRE( mul(float3x3(), float3()) == float3() );
    REQUIRE( mul(float3x4(), float4()) == float3() );
    REQUIRE( mul(float4x1(), float1()) == float4() );
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
    REQUIRE( outerprod(float1(), float1()) == float1x1() );
    REQUIRE( outerprod(float1(), float2()) == float1x2() );
    REQUIRE( outerprod(float1(), float3()) == float1x3() );
    REQUIRE( outerprod(float1(), float4()) == float1x4() );
    REQUIRE( outerprod(float2(), float1()) == float2x1() );
    REQUIRE( outerprod(float2(), float2()) == float2x2() );
    REQUIRE( outerprod(float2(), float3()) == float2x3() );
    REQUIRE( outerprod(float2(), float4()) == float2x4() );
    REQUIRE( outerprod(float3(), float1()) == float3x1() );
    REQUIRE( outerprod(float3(), float2()) == float3x2() );
    REQUIRE( outerprod(float3(), float3()) == float3x3() );
    REQUIRE( outerprod(float3(), float4()) == float3x4() );
    REQUIRE( outerprod(float4(), float1()) == float4x1() );
    REQUIRE( outerprod(float4(), float2()) == float4x2() );
    REQUIRE( outerprod(float4(), float3()) == float4x3() );
    REQUIRE( outerprod(float4(), float4()) == float4x4() );

    // Row vector x matrix products can be emulated using the transpose function
    REQUIRE( mul(transpose(float1()), float1x1()) == transpose(float1()) );
    REQUIRE( mul(transpose(float2()), float2x1()) == transpose(float1()) );
    REQUIRE( mul(transpose(float3()), float3x1()) == transpose(float1()) );
    REQUIRE( mul(transpose(float4()), float4x1()) == transpose(float1()) );
    REQUIRE( mul(transpose(float1()), float1x2()) == transpose(float2()) );
    REQUIRE( mul(transpose(float2()), float2x2()) == transpose(float2()) );
    REQUIRE( mul(transpose(float3()), float3x2()) == transpose(float2()) );
    REQUIRE( mul(transpose(float4()), float4x2()) == transpose(float2()) );
    REQUIRE( mul(transpose(float1()), float1x3()) == transpose(float3()) );
    REQUIRE( mul(transpose(float2()), float2x3()) == transpose(float3()) );
    REQUIRE( mul(transpose(float3()), float3x3()) == transpose(float3()) );
    REQUIRE( mul(transpose(float4()), float4x3()) == transpose(float3()) );
    REQUIRE( mul(transpose(float1()), float1x4()) == transpose(float4()) );
    REQUIRE( mul(transpose(float2()), float2x4()) == transpose(float4()) );
    REQUIRE( mul(transpose(float3()), float3x4()) == transpose(float4()) );
    REQUIRE( mul(transpose(float4()), float4x4()) == transpose(float4()) );
}

TEST_CASE_TEMPLATE( "matrix diagonal and trace are correct", T, double, float, int, short )
{
    random_number_generator rng;

    SUBCASE("1x1 matrices")
    {
        for(int i=0; i<reps; ++i)
        {
            const auto mat = rng.get<linalg::mat<T,1,1>>();
            CHECK(diagonal(mat) == linalg::vec<T,1>{mat[0][0]});
            CHECK(trace(mat) == mat[0][0]);
        }
    }

    SUBCASE("2x2 matrices")
    {
        for(int i=0; i<reps; ++i)
        {
            const auto mat = rng.get<linalg::mat<T,2,2>>();
            CHECK(diagonal(mat) == linalg::vec<T,2>{mat[0][0], mat[1][1]});
            CHECK(trace(mat) == mat[0][0] + mat[1][1]);
        }
    }

    SUBCASE("3x3 matrices")
    {
        for(int i=0; i<reps; ++i)
        {
            const auto mat = rng.get<linalg::mat<T,3,3>>();
            CHECK(diagonal(mat) == linalg::vec<T,3>{mat[0][0], mat[1][1], mat[2][2]});
            CHECK(trace(mat) == mat[0][0] + mat[1][1] + mat[2][2]);
        }
    }

    SUBCASE("4x4 matrices")
    {
        for(int i=0; i<reps; ++i)
        {
            const auto mat = rng.get<linalg::mat<T,4,4>>();
            CHECK(diagonal(mat) == linalg::vec<T,4>{mat[0][0], mat[1][1], mat[2][2], mat[3][3]});
            CHECK(trace(mat) == mat[0][0] + mat[1][1] + mat[2][2] + mat[3][3]);
        }
    }
}

TEST_CASE( "matrix adjugate and determinant are correct" )
{
    std::mt19937 rng;
    std::uniform_int_distribution<int> dist(-100,100);

    SUBCASE("1x1 matrices")
    {
        for(int i=0; i<reps; ++i)
        {
            const int1x1 id = linalg::identity;
            const int1x1 mat {
                {dist(rng)}
            };
            const int1x1 adj = adjugate(mat);
            const int det = determinant(mat);
            CHECK(mul(mat,adj) == id*det);
            CHECK(mul(adj,mat) == id*det);
            CHECK(mul(mat,adj) == det*id);
            CHECK(mul(adj,mat) == det*id);
            CHECK(determinant(adj) == 1);
        }
    }

    SUBCASE("2x2 matrices")
    {
        for(int i=0; i<reps; ++i)
        {
            const int2x2 id = linalg::identity;
            const int2x2 mat {
                {dist(rng), dist(rng)},
                {dist(rng), dist(rng)}
            };
            const int2x2 adj = adjugate(mat);
            const int det = determinant(mat);
            CHECK(mul(mat,adj) == id*det);
            CHECK(mul(adj,mat) == id*det);
            CHECK(mul(mat,adj) == det*id);
            CHECK(mul(adj,mat) == det*id);
            CHECK(determinant(adj) == det);
        }
    }

    SUBCASE("3x3 matrices")
    {
        for(int i=0; i<reps; ++i)
        {
            const int3x3 id = linalg::identity;
            const int3x3 mat {
                {dist(rng), dist(rng), dist(rng)},
                {dist(rng), dist(rng), dist(rng)},
                {dist(rng), dist(rng), dist(rng)}
            };
            const int3x3 adj = adjugate(mat);
            const int det = determinant(mat);
            CHECK(mul(mat,adj) == id*det);
            CHECK(mul(adj,mat) == id*det);
            CHECK(mul(mat,adj) == det*id);
            CHECK(mul(adj,mat) == det*id);
            CHECK(determinant(adj) == det*det);
        }
    }

    SUBCASE("4x4 matrices")
    {
        for(int i=0; i<reps; ++i)
        {
            const int4x4 id = linalg::identity;
            const int4x4 mat {
                {dist(rng), dist(rng), dist(rng), dist(rng)},
                {dist(rng), dist(rng), dist(rng), dist(rng)},
                {dist(rng), dist(rng), dist(rng), dist(rng)},
                {dist(rng), dist(rng), dist(rng), dist(rng)}        
            };
            const int4x4 adj = adjugate(mat);
            const int det = determinant(mat);
            CHECK(mul(mat,adj) == id*det);
            CHECK(mul(adj,mat) == id*det);
            CHECK(mul(mat,adj) == det*id);
            CHECK(mul(adj,mat) == det*id);
            CHECK(determinant(adj) == det*det*det);
        }
    }
}

TEST_CASE_TEMPLATE( "comatrix can be used to transform bivectors", T, float, double )
{
    random_number_generator rng;
    for(int i=0; i<reps; ++i)
    {
        const linalg::mat<T,3,3> m = rng;
        for(int j=0; j<reps; ++j)
        {
            const vec3<T> a = rng, b = rng;
            check_approx_equal( cross(mul(m,a), mul(m,b)), mul(comatrix(m),cross(a,b)) );
        }
    }
}

TEST_CASE( "matrix inverse is correct for trivial cases" )
{
    const float1x1 id1 {{1}};
    const float2x2 id2 {{1,0},{0,1}};
    const float3x3 id3 {{1,0,0},{0,1,0},{0,0,1}};
    const float4x4 id4 {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
    REQUIRE( diagonal(id1) == float1(1) );
    REQUIRE( diagonal(id2) == float2(1,1) );
    REQUIRE( diagonal(id3) == float3(1,1,1) );
    REQUIRE( diagonal(id4) == float4(1,1,1,1) );
    REQUIRE( transpose(id1) == id1 );
    REQUIRE( transpose(id2) == id2 );
    REQUIRE( transpose(id3) == id3 );
    REQUIRE( transpose(id4) == id4 );
    REQUIRE( inverse(id1) == id1 );
    REQUIRE( inverse(id2) == id2 );
    REQUIRE( inverse(id3) == id3 );
    REQUIRE( inverse(id4) == id4 );
    REQUIRE( adjugate(id1) == id1 );
    REQUIRE( adjugate(id2) == id2 );
    REQUIRE( adjugate(id3) == id3 );
    REQUIRE( adjugate(id4) == id4 );
    REQUIRE( determinant(id1) == 1.0f );
    REQUIRE( determinant(id2) == 1.0f );
    REQUIRE( determinant(id3) == 1.0f );
    REQUIRE( determinant(id4) == 1.0f );
}

TEST_CASE_TEMPLATE( "matrix inverse is correct for general case", T, float, double )
{
    const linalg::mat<T,4,4> mat {{1,2,3,4}, {5,-6,7,8}, {9,10,-11,12}, {13,14,15,-16}};
    const linalg::mat<T,4,4> inv = inverse(mat);
    const linalg::mat<T,4,4> id = mul(mat, inv);
    for(int j=0; j<4; ++j)
    {
        for(int i=0; i<4; ++i)
        {
            if(i == j) REQUIRE( id[j][i] == doctest::Approx(1.0f) );
            else REQUIRE( id[j][i] == doctest::Approx(0.0f) );
        }
    }
}

TEST_CASE_TEMPLATE( "linalg::identity functions correctly", T, double, float, int, short, unsigned int, unsigned short )
{
    const linalg::mat<T,1,1> a1 {linalg::identity}, b1 {{1}}, c1 {};
    const linalg::mat<T,2,2> a2 {linalg::identity}, b2 {{1,0},{0,1}}, c2 {};
    const linalg::mat<T,3,3> a3 {linalg::identity}, b3 {{1,0,0},{0,1,0},{0,0,1}}, c3 {};
    const linalg::mat<T,4,4> a4 {linalg::identity}, b4 {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}}, c4 {};

    REQUIRE(a1 == b1);
    REQUIRE(a2 == b2);
    REQUIRE(a3 == b3);
    REQUIRE(a4 == b4);
    REQUIRE(a1 != c1);
    REQUIRE(a2 != c2);
    REQUIRE(a3 != c3);
    REQUIRE(a4 != c4);
}
