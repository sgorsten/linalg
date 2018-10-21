#include "test-linalg.h"

TEST_CASE( "matrix multiplication produces correct result dimensions" )
{
    REQUIRE( float2x2() * float2() == float2() );
    REQUIRE( float2x3() * float3() == float2() );
    REQUIRE( float2x4() * float4() == float2() );
    REQUIRE( float3x2() * float2() == float3() );
    REQUIRE( float3x3() * float3() == float3() );
    REQUIRE( float3x4() * float4() == float3() );
    REQUIRE( float4x2() * float2() == float4() );
    REQUIRE( float4x3() * float3() == float4() );
    REQUIRE( float4x4() * float4() == float4() );

    REQUIRE( float2x2() * float2x2() == float2x2() );
    REQUIRE( float2x3() * float3x2() == float2x2() );
    REQUIRE( float2x4() * float4x2() == float2x2() );
    REQUIRE( float2x2() * float2x3() == float2x3() );
    REQUIRE( float2x3() * float3x3() == float2x3() );
    REQUIRE( float2x4() * float4x3() == float2x3() );
    REQUIRE( float2x2() * float2x4() == float2x4() );
    REQUIRE( float2x3() * float3x4() == float2x4() );
    REQUIRE( float2x4() * float4x4() == float2x4() );
    REQUIRE( float3x2() * float2x2() == float3x2() );
    REQUIRE( float3x3() * float3x2() == float3x2() );
    REQUIRE( float3x4() * float4x2() == float3x2() );
    REQUIRE( float3x2() * float2x3() == float3x3() );
    REQUIRE( float3x3() * float3x3() == float3x3() );
    REQUIRE( float3x4() * float4x3() == float3x3() );
    REQUIRE( float3x2() * float2x4() == float3x4() );
    REQUIRE( float3x3() * float3x4() == float3x4() );
    REQUIRE( float3x4() * float4x4() == float3x4() );
    REQUIRE( float4x2() * float2x2() == float4x2() );
    REQUIRE( float4x3() * float3x2() == float4x2() );
    REQUIRE( float4x4() * float4x2() == float4x2() );
    REQUIRE( float4x2() * float2x3() == float4x3() );
    REQUIRE( float4x3() * float3x3() == float4x3() );
    REQUIRE( float4x4() * float4x3() == float4x3() );
    REQUIRE( float4x2() * float2x4() == float4x4() );
    REQUIRE( float4x3() * float3x4() == float4x4() );
    REQUIRE( float4x4() * float4x4() == float4x4() );

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
            const int1x1 mat {{dist(rng)}};
            const int1x1 adj = adjugate(mat);
            const int det = determinant(mat);
            CHECK(mat * adj == id * det);
            CHECK(adj * mat == id * det);
            CHECK(mat * adj == det * id);
            CHECK(adj * mat == det * id);
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
            CHECK(mat * adj == id * det);
            CHECK(adj * mat == id * det);
            CHECK(mat * adj == det * id);
            CHECK(adj * mat == det * id);
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
            CHECK(mat * adj == id * det);
            CHECK(adj * mat == id * det);
            CHECK(mat * adj == det * id);
            CHECK(adj * mat == det * id);
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
            CHECK(mat * adj == id * det);
            CHECK(adj * mat == id * det);
            CHECK(mat * adj == det * id);
            CHECK(adj * mat == det * id);
        }
    }
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

TEST_CASE_TEMPLATE( "matrix inverse is correct for general case", T, float, double )
{
    const linalg::mat<T,4,4> mat {{1,2,3,4}, {5,-6,7,8}, {9,10,-11,12}, {13,14,15,-16}};
    const linalg::mat<T,4,4> inv = inverse(mat);
    const linalg::mat<T,4,4> id = mat * inv;
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
    const linalg::mat<T,2,2> a2 {linalg::identity}, b2 {{1,0},{0,1}}, c2 {};
    const linalg::mat<T,3,3> a3 {linalg::identity}, b3 {{1,0,0},{0,1,0},{0,0,1}}, c3 {};
    const linalg::mat<T,4,4> a4 {linalg::identity}, b4 {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}}, c4 {};

    REQUIRE(a2 == b2);
    REQUIRE(a3 == b3);
    REQUIRE(a4 == b4);
    REQUIRE(a2 != c2);
    REQUIRE(a3 != c3);
    REQUIRE(a4 != c4);
}


float3 transform_point(const float4x4 & m, const float3 & p) { const auto r = m*float4(p,1); return r.xyz/r.w; }

TEST_CASE( "Projection matrices behave as intended" )
{
    const float n = 0.1f, f = 10.0f;
    const float nx0 = -0.9f*n, ny0 = -0.6f*n, nx1 = 0.8f*n, ny1 = 0.7f*n, ncx = (nx0+nx1)/2, ncy = (ny0+ny1)/2;
    const float fx0 = -0.9f*f, fy0 = -0.6f*f, fx1 = 0.8f*f, fy1 = 0.7f*f, fcx = (fx0+fx1)/2, fcy = (fy0+fy1)/2;

    // Right handed OpenGL convention, x-right, y-up, z-back
    const float4x4 gl_rh = frustum_matrix(nx0, nx1, ny0, ny1, n, f, linalg::neg_z, linalg::neg_one_to_one); 
    check_approx_equal( transform_point(gl_rh, float3(ncx, ncy, -n)), float3( 0,  0, -1) );
    check_approx_equal( transform_point(gl_rh, float3(ncx, ny0, -n)), float3( 0, -1, -1) );
    check_approx_equal( transform_point(gl_rh, float3(ncx, ny1, -n)), float3( 0, +1, -1) );
    check_approx_equal( transform_point(gl_rh, float3(nx0, ncy, -n)), float3(-1,  0, -1) );
    check_approx_equal( transform_point(gl_rh, float3(nx1, ncy, -n)), float3(+1,  0, -1) );
    check_approx_equal( transform_point(gl_rh, float3(fcx, fcy, -f)), float3( 0,  0, +1) );
    check_approx_equal( transform_point(gl_rh, float3(fcx, fy0, -f)), float3( 0, -1, +1) );
    check_approx_equal( transform_point(gl_rh, float3(fcx, fy1, -f)), float3( 0, +1, +1) );
    check_approx_equal( transform_point(gl_rh, float3(fx0, fcy, -f)), float3(-1,  0, +1) );
    check_approx_equal( transform_point(gl_rh, float3(fx1, fcy, -f)), float3(+1,  0, +1) );

    // Left handed OpenGL convention, x-right, y-up, z-forward
    const float4x4 gl_lh = frustum_matrix(nx0, nx1, ny0, ny1, n, f, linalg::pos_z, linalg::neg_one_to_one);
    check_approx_equal( transform_point(gl_lh, float3(ncx, ncy, +n)), float3( 0,  0, -1) );
    check_approx_equal( transform_point(gl_lh, float3(ncx, ny0, +n)), float3( 0, -1, -1) );
    check_approx_equal( transform_point(gl_lh, float3(ncx, ny1, +n)), float3( 0, +1, -1) );
    check_approx_equal( transform_point(gl_lh, float3(nx0, ncy, +n)), float3(-1,  0, -1) );
    check_approx_equal( transform_point(gl_lh, float3(nx1, ncy, +n)), float3(+1,  0, -1) );
    check_approx_equal( transform_point(gl_lh, float3(fcx, fcy, +f)), float3( 0,  0, +1) );
    check_approx_equal( transform_point(gl_lh, float3(fcx, fy0, +f)), float3( 0, -1, +1) );
    check_approx_equal( transform_point(gl_lh, float3(fcx, fy1, +f)), float3( 0, +1, +1) );
    check_approx_equal( transform_point(gl_lh, float3(fx0, fcy, +f)), float3(-1,  0, +1) );
    check_approx_equal( transform_point(gl_lh, float3(fx1, fcy, +f)), float3(+1,  0, +1) );

    // Right handed Vulkan convention, x-right, y-down, z-forward
    const float4x4 vk_rh = frustum_matrix(nx0, nx1, ny0, ny1, n, f, linalg::pos_z, linalg::zero_to_one);
    check_approx_equal( transform_point(vk_rh, float3(ncx, ncy, +n)), float3( 0,  0, 0) );
    check_approx_equal( transform_point(vk_rh, float3(ncx, ny0, +n)), float3( 0, -1, 0) );
    check_approx_equal( transform_point(vk_rh, float3(ncx, ny1, +n)), float3( 0, +1, 0) );
    check_approx_equal( transform_point(vk_rh, float3(nx0, ncy, +n)), float3(-1,  0, 0) );
    check_approx_equal( transform_point(vk_rh, float3(nx1, ncy, +n)), float3(+1,  0, 0) );
    check_approx_equal( transform_point(vk_rh, float3(fcx, fcy, +f)), float3( 0,  0, 1) );
    check_approx_equal( transform_point(vk_rh, float3(fcx, fy0, +f)), float3( 0, -1, 1) );
    check_approx_equal( transform_point(vk_rh, float3(fcx, fy1, +f)), float3( 0, +1, 1) );
    check_approx_equal( transform_point(vk_rh, float3(fx0, fcy, +f)), float3(-1,  0, 1) );
    check_approx_equal( transform_point(vk_rh, float3(fx1, fcy, +f)), float3(+1,  0, 1) );

    // Left handed Vulkan convention, x-right, y-down, z-back
    const float4x4 vk_lh = frustum_matrix(nx0, nx1, ny0, ny1, n, f, linalg::neg_z, linalg::zero_to_one); 
    check_approx_equal( transform_point(vk_lh, float3(ncx, ncy, -n)), float3( 0,  0, 0) );
    check_approx_equal( transform_point(vk_lh, float3(ncx, ny0, -n)), float3( 0, -1, 0) );
    check_approx_equal( transform_point(vk_lh, float3(ncx, ny1, -n)), float3( 0, +1, 0) );
    check_approx_equal( transform_point(vk_lh, float3(nx0, ncy, -n)), float3(-1,  0, 0) );
    check_approx_equal( transform_point(vk_lh, float3(nx1, ncy, -n)), float3(+1,  0, 0) );
    check_approx_equal( transform_point(vk_lh, float3(fcx, fcy, -f)), float3( 0,  0, 1) );
    check_approx_equal( transform_point(vk_lh, float3(fcx, fy0, -f)), float3( 0, -1, 1) );
    check_approx_equal( transform_point(vk_lh, float3(fcx, fy1, -f)), float3( 0, +1, 1) );
    check_approx_equal( transform_point(vk_lh, float3(fx0, fcy, -f)), float3(-1,  0, 1) );
    check_approx_equal( transform_point(vk_lh, float3(fx1, fcy, -f)), float3(+1,  0, 1) );
}
