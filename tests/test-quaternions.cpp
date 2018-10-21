#include "test-linalg.h"


TEST_CASE_TEMPLATE("Quaternion exponent, logarithm, and power", T, float, double)
{
    constexpr T tau = T(6.2831853071795865), pi = tau/2;
    check_approx_equal( qexp(quat<T>{tau,0,0,0}), quat<T>{0,0,0,1} ); // e^tau*i = 1
    check_approx_equal( qexp(quat<T>{0,tau,0,0}), quat<T>{0,0,0,1} ); // e^tau*j = 1
    check_approx_equal( qexp(quat<T>{0,0,tau,0}), quat<T>{0,0,0,1} ); // e^tau*k = 1
    check_approx_equal( qexp(quat<T>{pi,0,0,0}), quat<T>{0,0,0,-1} ); // e^pi*i = -1
    check_approx_equal( qexp(quat<T>{0,pi,0,0}), quat<T>{0,0,0,-1} ); // e^pi*j = -1
    check_approx_equal( qexp(quat<T>{0,0,pi,0}), quat<T>{0,0,0,-1} ); // e^pi*k = -1

    random_number_generator rng;
    for(int i=0; i<reps; ++i)
    {
        // qexp, qlog, and qpow with a scalar should simply be the exp, log, pow of that scalar
        T s = rng, p = rng;
        check_approx_equal( qexp(quat<T>{0,0,0,s}  ), quat<T>{0,0,0,std::exp(s  )} );
        check_approx_equal( qlog(quat<T>{0,0,0,s}  ), quat<T>{0,0,0,std::log(s  )} );
        check_approx_equal( qpow(quat<T>{0,0,0,s},p), quat<T>{0,0,0,std::pow(s,p)} );

        // qpow with integral powers should have the effect of repeated multiplication
        quat<T> q = rng;
        check_approx_equal( qexp(qlog(q)), q );
        check_approx_equal( qpow(q, T(2)), q*q );
        check_approx_equal( qpow(q, T(3)), q*q*q );
        check_approx_equal( qpow(qpow(q, T(2)), T(3)), qpow(q, T(2*3)) );
    }   
}

TEST_CASE_TEMPLATE("qxdir(q) == qvq* where v={1,0,0,0}",  T, int, float, double)
{
    random_number_generator rng;
    for(int i=0; i<reps; ++i)
    {
        const quat<T> q = rng;
        const vec3<T> qx = qxdir(q), qx_ref = (q * quat<T>{1,0,0,0} * conjugate(q)).xyz();
        CHECK(qx.x == doctest::Approx(qx_ref.x));
        CHECK(qx.y == doctest::Approx(qx_ref.y));
        CHECK(qx.z == doctest::Approx(qx_ref.z));
    }
}

TEST_CASE_TEMPLATE("qzdir(q) == qvq* where v={0,1,0,0}",  T, int, float, double)
{
    random_number_generator rng;
    for(int i=0; i<reps; ++i)
    {
        const quat<T> q = rng;
        const vec3<T> qy = qydir(q), qy_ref = (q * quat<T>{0,1,0,0} * conjugate(q)).xyz();
        CHECK(qy.x == doctest::Approx(qy_ref.x));
        CHECK(qy.y == doctest::Approx(qy_ref.y));
        CHECK(qy.z == doctest::Approx(qy_ref.z));
    }
}

TEST_CASE_TEMPLATE("qzdir(q) == qvq* where v={0,0,1,0}",  T, int, float, double)
{
    random_number_generator rng;
    for(int i=0; i<reps; ++i)
    {
        const quat<T> q = rng;
        const vec3<T> qz = qzdir(q), qz_ref = (q * quat<T>{0,0,1,0} * conjugate(q)).xyz();
        CHECK(qz.x == doctest::Approx(qz_ref.x));
        CHECK(qz.y == doctest::Approx(qz_ref.y));
        CHECK(qz.z == doctest::Approx(qz_ref.z));
    }
}

TEST_CASE_TEMPLATE("qrot(q,v) == qvq*",  T, int, float, double)
{
    random_number_generator rng;
    for(int i=0; i<reps; ++i)
    {
        const quat<T> q = rng;
        for(int j=0; j<reps; ++j)
        {
            const vec3<T> v = rng, qr = qrot(q,v), qr_ref = (q * quat<T>{v,0} * conjugate(q)).xyz();
            CHECK(qr.x == doctest::Approx(qr_ref.x));
            CHECK(qr.y == doctest::Approx(qr_ref.y));
            CHECK(qr.z == doctest::Approx(qr_ref.z));
        }
    }
}

TEST_CASE_TEMPLATE( "rotation quaternions roundtrip with rotation matrices", T, float, double )
{
    std::mt19937 engine;
    std::normal_distribution<T> dist;

    for(int i=0; i<1000; ++i)
    {
        const quat<T> q = normalize(quat<T>(dist(engine), dist(engine), dist(engine), dist(engine)));
        const quat<T> q2 = rotation_quat(qmat(q));
        if(dot(q, q2) > 0) // q2 should either equal q or -q
        {
            REQUIRE( std::abs(q.x - q2.x) < 0.0001 );
            REQUIRE( std::abs(q.y - q2.y) < 0.0001 );
            REQUIRE( std::abs(q.z - q2.z) < 0.0001 );
            REQUIRE( std::abs(q.w - q2.w) < 0.0001 );
        }
        else
        {
            REQUIRE( std::abs(q.x + q2.x) < 0.0001 );
            REQUIRE( std::abs(q.y + q2.y) < 0.0001 );
            REQUIRE( std::abs(q.z + q2.z) < 0.0001 );
            REQUIRE( std::abs(q.w + q2.w) < 0.0001 );        
        }
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