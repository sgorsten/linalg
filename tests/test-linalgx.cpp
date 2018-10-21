#include "test-linalg.h"
#include "../linalgx.h"

TEST_CASE_TEMPLATE( "rotation quaternions roundtrip with rotation matrices", T, float, double )
{
    random_number_generator rng;
    for(int i=0; i<reps; ++i)
    {
        const quat<T> qq = rng, q = normalize(qq);        
        const quat<T> q2 = rotation_quat(qmat(q));
        check_approx_equal(q, dot(q, q2) > 0 ? q2 : -q2);
    }
}

TEST_CASE_TEMPLATE( "90 degree rotation matrices round trip with rotation quaternions", T, float, double )
{
    const linalg::mat<T,3,3> matrices[] {
        {{+1,0,0},{0,+1,0},{0,0,+1}},
        {{+1,0,0},{0,-1,0},{0,0,-1}},
        {{+1,0,0},{0,0,+1},{0,-1,0}},
        {{+1,0,0},{0,0,-1},{0,+1,0}},
        {{-1,0,0},{0,+1,0},{0,0,-1}},
        {{-1,0,0},{0,-1,0},{0,0,+1}},
        {{-1,0,0},{0,0,+1},{0,+1,0}},
        {{-1,0,0},{0,0,-1},{0,-1,0}},
        {{0,+1,0},{+1,0,0},{0,0,-1}},
        {{0,+1,0},{-1,0,0},{0,0,+1}},
        {{0,+1,0},{0,0,+1},{+1,0,0}},
        {{0,+1,0},{0,0,-1},{-1,0,0}},
        {{0,-1,0},{+1,0,0},{0,0,+1}},
        {{0,-1,0},{-1,0,0},{0,0,-1}},
        {{0,-1,0},{0,0,+1},{-1,0,0}},
        {{0,-1,0},{0,0,-1},{+1,0,0}},
        {{0,0,+1},{+1,0,0},{0,+1,0}},
        {{0,0,+1},{-1,0,0},{0,-1,0}},
        {{0,0,+1},{0,+1,0},{-1,0,0}},
        {{0,0,+1},{0,-1,0},{+1,0,0}},
        {{0,0,-1},{+1,0,0},{0,-1,0}},
        {{0,0,-1},{-1,0,0},{0,+1,0}},
        {{0,0,-1},{0,+1,0},{+1,0,0}},
        {{0,0,-1},{0,-1,0},{-1,0,0}},
    };
    for(auto & m : matrices)
    {
        check_approx_equal(m, qmat(rotation_quat(m)));
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
