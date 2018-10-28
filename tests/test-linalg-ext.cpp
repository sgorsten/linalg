#include "test-linalg.h"
#include "../linalg-ext.h"

namespace linalg
{
    static_assert(coord_system{coord_axis::forward, coord_axis::left,    coord_axis::up     }.is_right_handed(), "forward-left-up coordinate system is right handed");
    static_assert(coord_system{coord_axis::forward, coord_axis::left,    coord_axis::down   }.is_left_handed(),  "forward-left-down coordinate system is left handed");
    static_assert(coord_system{coord_axis::forward, coord_axis::right,   coord_axis::up     }.is_left_handed(),  "forward-right-up coordinate system is left handed");
    static_assert(coord_system{coord_axis::forward, coord_axis::right,   coord_axis::down   }.is_right_handed(), "forward-right-down coordinate system is right handed");
    static_assert(coord_system{coord_axis::forward, coord_axis::up,      coord_axis::left   }.is_left_handed(),  "forward-up-left coordinate system is left handed");
    static_assert(coord_system{coord_axis::forward, coord_axis::up,      coord_axis::right  }.is_right_handed(), "forward-up-right coordinate system is right handed");
    static_assert(coord_system{coord_axis::forward, coord_axis::down,    coord_axis::left   }.is_right_handed(), "forward-down-left coordinate system is right handed");
    static_assert(coord_system{coord_axis::forward, coord_axis::down,    coord_axis::right  }.is_left_handed(),  "forward-down-right coordinate system is left handed");
    static_assert(coord_system{coord_axis::back,    coord_axis::left,    coord_axis::up     }.is_left_handed(),  "back-left-up coordinate system is left handed");
    static_assert(coord_system{coord_axis::back,    coord_axis::left,    coord_axis::down   }.is_right_handed(), "back-left-down coordinate system is right handed");
    static_assert(coord_system{coord_axis::back,    coord_axis::right,   coord_axis::up     }.is_right_handed(), "back-right-up coordinate system is right handed");
    static_assert(coord_system{coord_axis::back,    coord_axis::right,   coord_axis::down   }.is_left_handed(),  "back-right-down coordinate system is left handed");
    static_assert(coord_system{coord_axis::back,    coord_axis::up,      coord_axis::left   }.is_right_handed(), "back-up-left coordinate system is right handed");
    static_assert(coord_system{coord_axis::back,    coord_axis::up,      coord_axis::right  }.is_left_handed(),  "back-up-right coordinate system is left handed");
    static_assert(coord_system{coord_axis::back,    coord_axis::down,    coord_axis::left   }.is_left_handed(),  "back-down-left coordinate system is left handed");
    static_assert(coord_system{coord_axis::back,    coord_axis::down,    coord_axis::right  }.is_right_handed(), "back-down-right coordinate system is right handed");
    static_assert(coord_system{coord_axis::left,    coord_axis::forward, coord_axis::up     }.is_left_handed(),  "left-forward-up coordinate system is left handed");
    static_assert(coord_system{coord_axis::left,    coord_axis::forward, coord_axis::down   }.is_right_handed(), "left-forward-down coordinate system is right handed");
    static_assert(coord_system{coord_axis::left,    coord_axis::back,    coord_axis::up     }.is_right_handed(), "left-back-up coordinate system is right handed");
    static_assert(coord_system{coord_axis::left,    coord_axis::back,    coord_axis::down   }.is_left_handed(),  "left-back-down coordinate system is left handed");
    static_assert(coord_system{coord_axis::left,    coord_axis::up,      coord_axis::forward}.is_right_handed(), "left-up-forward coordinate system is right handed");
    static_assert(coord_system{coord_axis::left,    coord_axis::up,      coord_axis::back   }.is_left_handed(),  "left-up-back coordinate system is left handed");
    static_assert(coord_system{coord_axis::left,    coord_axis::down,    coord_axis::forward}.is_left_handed(),  "left-down-forward coordinate system is left handed");
    static_assert(coord_system{coord_axis::left,    coord_axis::down,    coord_axis::back   }.is_right_handed(), "left-down-back coordinate system is right handed");                                                
    static_assert(coord_system{coord_axis::right,   coord_axis::forward, coord_axis::up     }.is_right_handed(), "right-forward-up coordinate system is right handed");
    static_assert(coord_system{coord_axis::right,   coord_axis::forward, coord_axis::down   }.is_left_handed(),  "right-forward-down coordinate system is left handed");
    static_assert(coord_system{coord_axis::right,   coord_axis::back,    coord_axis::up     }.is_left_handed(),  "right-back-up coordinate system is left handed");
    static_assert(coord_system{coord_axis::right,   coord_axis::back,    coord_axis::down   }.is_right_handed(), "right-back-down coordinate system is right handed");
    static_assert(coord_system{coord_axis::right,   coord_axis::up,      coord_axis::forward}.is_left_handed(),  "right-up-forward coordinate system is left handed");
    static_assert(coord_system{coord_axis::right,   coord_axis::up,      coord_axis::back   }.is_right_handed(), "right-up-back coordinate system is right handed");
    static_assert(coord_system{coord_axis::right,   coord_axis::down,    coord_axis::forward}.is_right_handed(), "right-down-forward coordinate system is right handed");
    static_assert(coord_system{coord_axis::right,   coord_axis::down,    coord_axis::back   }.is_left_handed(),  "right-down-back coordinate system is left handed");
    static_assert(coord_system{coord_axis::up,      coord_axis::left,    coord_axis::forward}.is_left_handed(),  "up-left-forward coordinate system is left handed");
    static_assert(coord_system{coord_axis::up,      coord_axis::left,    coord_axis::back   }.is_right_handed(), "up-left-back coordinate system is right handed");
    static_assert(coord_system{coord_axis::up,      coord_axis::right,   coord_axis::forward}.is_right_handed(), "up-right-forward coordinate system is right handed");
    static_assert(coord_system{coord_axis::up,      coord_axis::right,   coord_axis::back   }.is_left_handed(),  "up-right-back coordinate system is left handed");
    static_assert(coord_system{coord_axis::up,      coord_axis::forward, coord_axis::left   }.is_right_handed(), "up-forward-left coordinate system is right handed");
    static_assert(coord_system{coord_axis::up,      coord_axis::forward, coord_axis::right  }.is_left_handed(),  "up-forward-right coordinate system is left handed");
    static_assert(coord_system{coord_axis::up,      coord_axis::back,    coord_axis::left   }.is_left_handed(),  "up-back-left coordinate system is left handed");
    static_assert(coord_system{coord_axis::up,      coord_axis::back,    coord_axis::right  }.is_right_handed(), "up-back-right coordinate system is right handed");
    static_assert(coord_system{coord_axis::down,    coord_axis::left,    coord_axis::forward}.is_right_handed(), "down-left-forward coordinate system is right handed");
    static_assert(coord_system{coord_axis::down,    coord_axis::left,    coord_axis::back   }.is_left_handed(),  "down-left-back coordinate system is left handed");
    static_assert(coord_system{coord_axis::down,    coord_axis::right,   coord_axis::forward}.is_left_handed(),  "down-right-forward coordinate system is left handed");
    static_assert(coord_system{coord_axis::down,    coord_axis::right,   coord_axis::back   }.is_right_handed(), "down-right-back coordinate system is right handed");
    static_assert(coord_system{coord_axis::down,    coord_axis::forward, coord_axis::left   }.is_left_handed(),  "down-forward-left coordinate system is left handed");
    static_assert(coord_system{coord_axis::down,    coord_axis::forward, coord_axis::right  }.is_right_handed(), "down-forward-right coordinate system is right handed");
    static_assert(coord_system{coord_axis::down,    coord_axis::back,    coord_axis::left   }.is_right_handed(), "down-back-left coordinate system is right handed");
    static_assert(coord_system{coord_axis::down,    coord_axis::back,    coord_axis::right  }.is_left_handed(),  "down-back-right coordinate system is left handed");
}

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

TEST_CASE_TEMPLATE("Test linear and homogeneous transformations", T, float, double)
{
    random_number_generator rng;
    const linalg::coord_system ruf{linalg::coord_axis::right, linalg::coord_axis::up, linalg::coord_axis::forward};
    const linalg::coord_system rfu{linalg::coord_axis::right, linalg::coord_axis::forward, linalg::coord_axis::up};
    for(int i=0; i<reps; ++i)
    {
        const vec3<T> a = rng, b = rng, u = normalize(a);
        const vec3<T> c = cross(a,b), n = normalize(c);
        
        const T uniform_scale = rng;
        const vec3<T> scale_factors = rng;
        const quat<T> rotation_quat = normalize(rng.get<quat<T>>());
        const vec3<T> translation_vec = rng;

        // Check linear transformations of vectors
        using linear = linalg::linear_transformation<T,3>;
        check_approx_equal(transform_vector(linear::identity(), a),              a);
        check_approx_equal(transform_vector(linear::coord_change(rfu,ruf), a),   a.xzy());
        check_approx_equal(transform_vector(linear::scaling(uniform_scale), a),  a*uniform_scale);
        check_approx_equal(transform_vector(linear::scaling(scale_factors), a),  a*scale_factors);        
        check_approx_equal(transform_vector(linear::rotation(rotation_quat), a), qrot(rotation_quat,a));

        // Check linear transformation of directions
        check_approx_equal(transform_direction(linear::identity(), u),              normalize(a));
        check_approx_equal(transform_direction(linear::coord_change(rfu,ruf), u),   normalize(a.xzy()));
        check_approx_equal(transform_direction(linear::scaling(uniform_scale), u),  normalize(a*uniform_scale));
        check_approx_equal(transform_direction(linear::scaling(scale_factors), u),  normalize(a*scale_factors));        
        check_approx_equal(transform_direction(linear::rotation(rotation_quat), u), normalize(qrot(rotation_quat,a)));

        // Check linear transformation of bivectors
        check_approx_equal(transform_bivector(linear::identity(), c),              cross(a,b));
        check_approx_equal(transform_bivector(linear::coord_change(rfu,ruf), c),   cross(a.xzy(),b.xzy()));
        check_approx_equal(transform_bivector(linear::scaling(uniform_scale), c),  cross(a*uniform_scale,b*uniform_scale));
        check_approx_equal(transform_bivector(linear::scaling(scale_factors), c),  cross(a*scale_factors,b*scale_factors));        
        check_approx_equal(transform_bivector(linear::rotation(rotation_quat), c), cross(qrot(rotation_quat,a),qrot(rotation_quat,b)));

        // Check linear transformation of normals
        check_approx_equal(transform_normal(linear::identity(), n),              normalize(cross(a,b)));
        check_approx_equal(transform_normal(linear::coord_change(rfu,ruf), n),   normalize(cross(a.xzy(),b.xzy())));
        check_approx_equal(transform_normal(linear::scaling(uniform_scale), n),  normalize(cross(a*uniform_scale,b*uniform_scale)));
        check_approx_equal(transform_normal(linear::scaling(scale_factors), n),  normalize(cross(a*scale_factors,b*scale_factors)));        
        check_approx_equal(transform_normal(linear::rotation(rotation_quat), n), normalize(cross(qrot(rotation_quat,a),qrot(rotation_quat,b))));

        // Check linear transformations of points
        check_approx_equal(transform_point(linear::identity(), a),              a);
        check_approx_equal(transform_point(linear::coord_change(rfu,ruf), a),   a.xzy());
        check_approx_equal(transform_point(linear::scaling(uniform_scale), a),  a*uniform_scale);
        check_approx_equal(transform_point(linear::scaling(scale_factors), a),  a*scale_factors);        
        check_approx_equal(transform_point(linear::rotation(rotation_quat), a), qrot(rotation_quat,a));

        // Check homogeneous transformations of vectors
        using homog = linalg::homogeneous_transformation<T,3>;
        check_approx_equal(transform_vector(homog::identity(), a),                          a);
        check_approx_equal(transform_vector(homog::coord_change(rfu,ruf), a),               a.xzy());
        check_approx_equal(transform_vector(homog::scaling(uniform_scale), a),              a*uniform_scale);
        check_approx_equal(transform_vector(homog::scaling(scale_factors), a),              a*scale_factors);        
        check_approx_equal(transform_vector(homog::rotation(rotation_quat), a),             qrot(rotation_quat,a));
        check_approx_equal(transform_vector(homog::translation(translation_vec), a),        a);
        check_approx_equal(transform_vector(homog::pose(rotation_quat,translation_vec), a), qrot(rotation_quat,a));
        
        // Check homogeneous transformations of directions
        check_approx_equal(transform_direction(homog::identity(), u),                          normalize(a));
        check_approx_equal(transform_direction(homog::coord_change(rfu,ruf), u),               normalize(a.xzy()));
        check_approx_equal(transform_direction(homog::scaling(uniform_scale), u),              normalize(a*uniform_scale));
        check_approx_equal(transform_direction(homog::scaling(scale_factors), u),              normalize(a*scale_factors));
        check_approx_equal(transform_direction(homog::rotation(rotation_quat), u),             normalize(qrot(rotation_quat,a)));
        check_approx_equal(transform_direction(homog::translation(translation_vec), u),        normalize(a));
        check_approx_equal(transform_direction(homog::pose(rotation_quat,translation_vec), u), normalize(qrot(rotation_quat,a)));

        // Check homogeneous transformation of bivectors
        check_approx_equal(transform_bivector(homog::identity(), c),                          cross(a,b));
        check_approx_equal(transform_bivector(homog::coord_change(rfu,ruf), c),               cross(a.xzy(),b.xzy()));
        check_approx_equal(transform_bivector(homog::scaling(uniform_scale), c),              cross(a*uniform_scale,b*uniform_scale));
        check_approx_equal(transform_bivector(homog::scaling(scale_factors), c),              cross(a*scale_factors,b*scale_factors));        
        check_approx_equal(transform_bivector(homog::rotation(rotation_quat), c),             cross(qrot(rotation_quat,a),qrot(rotation_quat,b)));
        check_approx_equal(transform_bivector(homog::translation(translation_vec), c),        cross(a,b));
        check_approx_equal(transform_bivector(homog::pose(rotation_quat,translation_vec), c), cross(qrot(rotation_quat,a),qrot(rotation_quat,b)));

        // Check homogeneous transformation of normals
        check_approx_equal(transform_normal(homog::identity(), c),                          normalize(cross(a,b)));
        check_approx_equal(transform_normal(homog::coord_change(rfu,ruf), c),               normalize(cross(a.xzy(),b.xzy())));
        check_approx_equal(transform_normal(homog::scaling(uniform_scale), c),              normalize(cross(a*uniform_scale,b*uniform_scale)));
        check_approx_equal(transform_normal(homog::scaling(scale_factors), c),              normalize(cross(a*scale_factors,b*scale_factors)));
        check_approx_equal(transform_normal(homog::rotation(rotation_quat), c),             normalize(cross(qrot(rotation_quat,a),qrot(rotation_quat,b))));
        check_approx_equal(transform_normal(homog::translation(translation_vec), c),        normalize(cross(a,b)));
        check_approx_equal(transform_normal(homog::pose(rotation_quat,translation_vec), c), normalize(cross(qrot(rotation_quat,a),qrot(rotation_quat,b))));

        // Check homogeneous transformations of points
        check_approx_equal(transform_point(homog::identity(), a),                          a);
        check_approx_equal(transform_point(homog::coord_change(rfu,ruf), a),               a.xzy());
        check_approx_equal(transform_point(homog::scaling(uniform_scale), a),              a*uniform_scale);
        check_approx_equal(transform_point(homog::scaling(scale_factors), a),              a*scale_factors);        
        check_approx_equal(transform_point(homog::rotation(rotation_quat), a),             qrot(rotation_quat,a));
        check_approx_equal(transform_point(homog::translation(translation_vec), a),        a+translation_vec);
        check_approx_equal(transform_point(homog::pose(rotation_quat,translation_vec), a), qrot(rotation_quat,a)+translation_vec);
    }
}

TEST_CASE( "Projection matrices behave as intended" )
{
    const float n = 0.1f, f = 10.0f;
    const float nx0 = -0.9f*n, ny0 = -0.6f*n, nx1 = 0.8f*n, ny1 = 0.7f*n, ncx = (nx0+nx1)/2, ncy = (ny0+ny1)/2;
    const float fx0 = -0.9f*f, fy0 = -0.6f*f, fx1 = 0.8f*f, fy1 = 0.7f*f, fcx = (fx0+fx1)/2, fcy = (fy0+fy1)/2;

    using transform = linalg::homogeneous_transformation<float,3>;

    // Right handed OpenGL convention, x-right, y-up, z-back
    const float4x4 gl_rh = transform::frustum(nx0, nx1, ny0, ny1, n, f, linalg::neg_z, linalg::neg_one_to_one); 
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
    const float4x4 gl_lh = transform::frustum(nx0, nx1, ny0, ny1, n, f, linalg::pos_z, linalg::neg_one_to_one);
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
    const float4x4 vk_rh = transform::frustum(nx0, nx1, ny0, ny1, n, f, linalg::pos_z, linalg::zero_to_one);
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
    const float4x4 vk_lh = transform::frustum(nx0, nx1, ny0, ny1, n, f, linalg::neg_z, linalg::zero_to_one); 
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
