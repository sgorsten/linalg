#include "test-linalg.h"

TEST_CASE_TEMPLATE("qxdir(q) == qvq* where v={1,0,0,0}",  T, floating_point_types)
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

TEST_CASE_TEMPLATE("qzdir(q) == qvq* where v={0,1,0,0}",  T, floating_point_types)
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

TEST_CASE_TEMPLATE("qzdir(q) == qvq* where v={0,0,1,0}",  T, floating_point_types)
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

TEST_CASE_TEMPLATE("qrot(q,v) == qvq*",  T, floating_point_types)
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