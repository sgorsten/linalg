#include "test-linalg.h"

float3 transform_point(const float4x4 & m, const float3 & p) { const auto r = mul(m,float4(p,1)); return r.xyz()/r.w; }

void test_lookat_matrix(const float3 & eye, const float3 & center, const float3 & view_y, linalg::fwd_axis fwd,
    const float3 & expected_origin, const float3 & expected_xdir, const float3 & expected_ydir, const float3 & expected_zdir)
{
    DOCTEST_CAPTURE(eye);
    DOCTEST_CAPTURE(center);
    DOCTEST_CAPTURE(view_y);
    DOCTEST_CAPTURE(fwd);

    auto view_matrix = lookat_matrix(eye, center, view_y, linalg::neg_z);
    CHECK( transform_point(view_matrix, center) == approx(float3{0,0,0}) );
    CHECK( transform_point(view_matrix, center + normalize(expected_xdir)) == approx(float3{1,0,0}) );
    CHECK( transform_point(view_matrix, center + normalize(expected_ydir)) == approx(float3{0,1,0}) );
    CHECK( transform_point(view_matrix, center + normalize(expected_zdir)) == approx(float3{0,0,1}) );
}

TEST_CASE( "Lookat matrices behave as intended" )
{
    const struct lookat_matrix_test_case
    {
        float3 eye, center, view_y;
        linalg::fwd_axis fwd;
        float3 expected_xdir, expected_ydir, expected_zdir;
    } cases[]
    {
        // Tests for a y-up world with an x-right, y-up, z-back NDC
        {{ 0, 0,+5}, { 0, 0,-5}, {0,1,0}, linalg::neg_z,    {+1, 0, 0}, { 0,+1, 0}, { 0, 0,+1}}, // Looking forward through origin along -z
        {{+5, 0, 0}, {-5, 0, 0}, {0,1,0}, linalg::neg_z,    { 0, 0,-1}, { 0,+1, 0}, {+1, 0, 0}}, // Looking forward through origin along -x
        {{ 0, 0,-5}, { 0, 0,+5}, {0,1,0}, linalg::neg_z,    {-1, 0, 0}, { 0,+1, 0}, { 0, 0,-1}}, // Looking forward through origin along +z
        {{-5, 0, 0}, {+5, 0, 0}, {0,1,0}, linalg::neg_z,    { 0, 0,+1}, { 0,+1, 0}, {-1, 0, 0}}, // Looking forward through origin along +x

        {{ 0,+3,+4}, { 0, 0, 0}, {0,1,0}, linalg::neg_z,    {+1, 0, 0}, { 0,+8,-6}, { 0,+6,+8}}, // Looking slightly down through origin along -z
        {{+4,+3, 0}, { 0, 0, 0}, {0,1,0}, linalg::neg_z,    { 0, 0,-1}, {-6,+8, 0}, {+8,+6, 0}}, // Looking slightly down through origin along -x
        {{ 0,+3,-4}, { 0, 0, 0}, {0,1,0}, linalg::neg_z,    {-1, 0, 0}, { 0,+8,+6}, { 0,+6,-8}}, // Looking slightly down through origin along +z
        {{-4,+3, 0}, { 0, 0, 0}, {0,1,0}, linalg::neg_z,    { 0, 0,+1}, {+6,+8, 0}, {-8,+6, 0}}, // Looking slightly down through origin along +x

        {{ 0,-3,+4}, { 0, 0, 0}, {0,1,0}, linalg::neg_z,    {+1, 0, 0}, { 0,+8,+6}, { 0,-6,+8}}, // Looking slightly up through origin along -z
        {{+4,-3, 0}, { 0, 0, 0}, {0,1,0}, linalg::neg_z,    { 0, 0,-1}, {+6,+8, 0}, {+8,-6, 0}}, // Looking slightly up through origin along -x
        {{ 0,-3,-4}, { 0, 0, 0}, {0,1,0}, linalg::neg_z,    {-1, 0, 0}, { 0,+8,-6}, { 0,-6,-8}}, // Looking slightly up through origin along +z
        {{-4,-3, 0}, { 0, 0, 0}, {0,1,0}, linalg::neg_z,    { 0, 0,+1}, {-6,+8, 0}, {-8,-6, 0}}, // Looking slightly up through origin along +x

        // Tests for a z-up world with an x-right, y-up, z-back NDC
        {{ 0,+5, 0}, { 0,-5, 0}, {0,0,1}, linalg::neg_z,    {-1, 0, 0}, { 0, 0,+1}, { 0,+1, 0}}, // Looking forward through origin along -y
        {{+5, 0, 0}, {-5, 0, 0}, {0,0,1}, linalg::neg_z,    { 0,+1, 0}, { 0, 0,+1}, {+1, 0, 0}}, // Looking forward through origin along -x
        {{ 0,-5, 0}, { 0,+5, 0}, {0,0,1}, linalg::neg_z,    {+1, 0, 0}, { 0, 0,+1}, { 0,-1, 0}}, // Looking forward through origin along +y
        {{-5, 0, 0}, {+5, 0, 0}, {0,0,1}, linalg::neg_z,    { 0,-1, 0}, { 0, 0,+1}, {-1, 0, 0}}, // Looking forward through origin along +x

        {{ 0,+4,+3}, { 0, 0, 0}, {0,0,1}, linalg::neg_z,    {-1, 0, 0}, { 0,-6,+8}, { 0,+8,+6}}, // Looking slightly down through origin along -y
        {{+4, 0,+3}, { 0, 0, 0}, {0,0,1}, linalg::neg_z,    { 0,+1, 0}, {-6, 0,+8}, {+8, 0,+6}}, // Looking slightly down through origin along -x
        {{ 0,-4,+3}, { 0, 0, 0}, {0,0,1}, linalg::neg_z,    {+1, 0, 0}, { 0,+6,+8}, { 0,-8,+6}}, // Looking slightly down through origin along +y
        {{-4, 0,+3}, { 0, 0, 0}, {0,0,1}, linalg::neg_z,    { 0,-1, 0}, {+6, 0,+8}, {-8, 0,+6}}, // Looking slightly down through origin along +x

        {{ 0,+4,-3}, { 0, 0, 0}, {0,0,1}, linalg::neg_z,    {-1, 0, 0}, { 0,+6,+8}, { 0,+8,-6}}, // Looking slightly up through origin along -y
        {{+4, 0,-3}, { 0, 0, 0}, {0,0,1}, linalg::neg_z,    { 0,+1, 0}, {+6, 0,+8}, {+8, 0,-6}}, // Looking slightly up through origin along -x
        {{ 0,-4,-3}, { 0, 0, 0}, {0,0,1}, linalg::neg_z,    {+1, 0, 0}, { 0,-6,+8}, { 0,-8,-6}}, // Looking slightly up through origin along +y
        {{-4, 0,-3}, { 0, 0, 0}, {0,0,1}, linalg::neg_z,    { 0,-1, 0}, {-6, 0,+8}, {-8, 0,-6}}, // Looking slightly up through origin along +x

        // Tests for a y-down world with an x-right, y-down, z-forward NDC
        {{ 0, 0,+5}, { 0, 0,-5}, {0,1,0}, linalg::pos_z,    {-1, 0, 0}, { 0,+1, 0}, { 0, 0,-1}}, // Looking forward through origin along -z
        {{+5, 0, 0}, {-5, 0, 0}, {0,1,0}, linalg::pos_z,    { 0, 0,+1}, { 0,+1, 0}, {-1, 0, 0}}, // Looking forward through origin along -x
        {{ 0, 0,-5}, { 0, 0,+5}, {0,1,0}, linalg::pos_z,    {+1, 0, 0}, { 0,+1, 0}, { 0, 0,+1}}, // Looking forward through origin along +z
        {{-5, 0, 0}, {+5, 0, 0}, {0,1,0}, linalg::pos_z,    { 0, 0,-1}, { 0,+1, 0}, {+1, 0, 0}}, // Looking forward through origin along +x

        {{ 0,+3,+4}, { 0, 0, 0}, {0,1,0}, linalg::pos_z,    {-1, 0, 0}, { 0,+8,-6}, { 0,-6,-8}}, // Looking slightly up through origin along -z
        {{+4,+3, 0}, { 0, 0, 0}, {0,1,0}, linalg::pos_z,    { 0, 0,+1}, {-6,+8, 0}, {-8,-6, 0}}, // Looking slightly up through origin along -x
        {{ 0,+3,-4}, { 0, 0, 0}, {0,1,0}, linalg::pos_z,    {+1, 0, 0}, { 0,+8,+6}, { 0,-6,+8}}, // Looking slightly up through origin along +z
        {{-4,+3, 0}, { 0, 0, 0}, {0,1,0}, linalg::pos_z,    { 0, 0,-1}, {+6,+8, 0}, {+8,-6, 0}}, // Looking slightly up through origin along +x

        {{ 0,-3,+4}, { 0, 0, 0}, {0,1,0}, linalg::pos_z,    {-1, 0, 0}, { 0,+8,+6}, { 0,+6,-8}}, // Looking slightly down through origin along -z
        {{+4,-3, 0}, { 0, 0, 0}, {0,1,0}, linalg::pos_z,    { 0, 0,+1}, {+6,+8, 0}, {-8,+6, 0}}, // Looking slightly down through origin along -x
        {{ 0,-3,-4}, { 0, 0, 0}, {0,1,0}, linalg::pos_z,    {+1, 0, 0}, { 0,+8,-6}, { 0,+6,+8}}, // Looking slightly down through origin along +z
        {{-4,-3, 0}, { 0, 0, 0}, {0,1,0}, linalg::pos_z,    { 0, 0,-1}, {-6,+8, 0}, {+8,+6, 0}}, // Looking slightly down through origin along +x

        // Tests for a z-up world with an x-right, y-down, z-forward NDC
        {{ 0,+5, 0}, { 0,-5, 0}, {0,0,-1}, linalg::pos_z,    {-1, 0, 0}, { 0, 0,-1}, { 0,-1, 0}}, // Looking forward through origin along -y
        {{+5, 0, 0}, {-5, 0, 0}, {0,0,-1}, linalg::pos_z,    { 0,+1, 0}, { 0, 0,-1}, {-1, 0, 0}}, // Looking forward through origin along -x
        {{ 0,-5, 0}, { 0,+5, 0}, {0,0,-1}, linalg::pos_z,    {+1, 0, 0}, { 0, 0,-1}, { 0,+1, 0}}, // Looking forward through origin along +y
        {{-5, 0, 0}, {+5, 0, 0}, {0,0,-1}, linalg::pos_z,    { 0,-1, 0}, { 0, 0,-1}, {+1, 0, 0}}, // Looking forward through origin along +x

        {{ 0,+4,+3}, { 0, 0, 0}, {0,0,-1}, linalg::pos_z,    {-1, 0, 0}, { 0,+6,-8}, { 0,-8,-6}}, // Looking slightly down through origin along -y
        {{+4, 0,+3}, { 0, 0, 0}, {0,0,-1}, linalg::pos_z,    { 0,+1, 0}, {+6, 0,-8}, {-8, 0,-6}}, // Looking slightly down through origin along -x
        {{ 0,-4,+3}, { 0, 0, 0}, {0,0,-1}, linalg::pos_z,    {+1, 0, 0}, { 0,-6,-8}, { 0,+8,-6}}, // Looking slightly down through origin along +z
        {{-4, 0,+3}, { 0, 0, 0}, {0,0,-1}, linalg::pos_z,    { 0,-1, 0}, {-6, 0,-8}, {+8, 0,-6}}, // Looking slightly down through origin along +x

        {{ 0,+4,-3}, { 0, 0, 0}, {0,0,-1}, linalg::pos_z,    {-1, 0, 0}, { 0,-6,-8}, { 0,-8,+6}}, // Looking slightly up through origin along -y
        {{+4, 0,-3}, { 0, 0, 0}, {0,0,-1}, linalg::pos_z,    { 0,+1, 0}, {-6, 0,-8}, {-8, 0,+6}}, // Looking slightly up through origin along -x
        {{ 0,-4,-3}, { 0, 0, 0}, {0,0,-1}, linalg::pos_z,    {+1, 0, 0}, { 0,+6,-8}, { 0,+8,+6}}, // Looking slightly up through origin along +y
        {{-4, 0,-3}, { 0, 0, 0}, {0,0,-1}, linalg::pos_z,    { 0,-1, 0}, {+6, 0,-8}, {+8, 0,+6}}, // Looking slightly up through origin along +x
    };
    for(auto & c : cases)
    {
        const float3 eye = c.eye, center = c.center, view_y = c.view_y;
        const auto fwd = c.fwd;
        const float3 expected_xdir = normalize(c.expected_xdir), expected_ydir = normalize(c.expected_ydir), expected_zdir = normalize(c.expected_zdir);

        CAPTURE(eye);
        CAPTURE(center);
        CAPTURE(view_y);
        CAPTURE(fwd);

        CHECK(transform_point(lookat_matrix(eye, center, view_y, fwd), eye) == approx(float3{0,0,0}));

        {
            CAPTURE(expected_xdir);
            CHECK(transform_point(lookat_matrix(eye, center, view_y, fwd), eye - expected_xdir) == approx(float3{-1,0,0}));
            CHECK(transform_point(lookat_matrix(eye, center, view_y, fwd), eye + expected_xdir) == approx(float3{+1,0,0}));
        }

        {
            CAPTURE(expected_ydir);
            CHECK(transform_point(lookat_matrix(eye, center, view_y, fwd), eye - expected_ydir) == approx(float3{0,-1,0}));
            CHECK(transform_point(lookat_matrix(eye, center, view_y, fwd), eye + expected_ydir) == approx(float3{0,+1,0}));
        }

        {
            CAPTURE(expected_zdir);
            CHECK(transform_point(lookat_matrix(eye, center, view_y, fwd), eye - expected_zdir) == approx(float3{0,0,-1}));
            CHECK(transform_point(lookat_matrix(eye, center, view_y, fwd), eye + expected_zdir) == approx(float3{0,0,+1}));
        }
    }
}

TEST_CASE( "Projection matrices behave as intended" )
{
    const float n = 0.1f, f = 10.0f;
    const float nx0 = -0.9f*n, ny0 = -6*n, nx1 = 8*n, ny1 = 0.7f*n, ncx = (nx0+nx1)/2, ncy = (ny0+ny1)/2;
    const float fx0 = -0.9f*f, fy0 = -6*f, fx1 = 8*f, fy1 = 0.7f*f, fcx = (fx0+fx1)/2, fcy = (fy0+fy1)/2;

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