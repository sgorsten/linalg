#include "test-linalg.h"

TEST_CASE_TEMPLATE("Certain quaternion operations mirror the equivalent operations on 4D vectors", T, float, double)
{
    random_number_generator rng;
    for(int i=0; i<reps; ++i)
    {
        // Check that constructing quats from vecs retains the exact values
        const vec4<T> va = rng, vb = rng;
        const quat<T> qa = quat<T>{va}, qb = quat<T>{vb};
        CHECK( qa.x == va[0] );
        CHECK( qa.y == va[1] );
        CHECK( qa.z == va[2] );
        CHECK( qa.w == va[3] );
        CHECK( qb.x == vb[0] );
        CHECK( qb.y == vb[1] );
        CHECK( qb.z == vb[2] );
        CHECK( qb.w == vb[3] );

        // Check dot product and length
        CHECK( dot(va,vb) == dot(qa,qb) );
        CHECK( length2(va) == length2(qa) );
        CHECK( length2(vb) == length2(qb) );
        CHECK( length(va) == length(qa) );
        CHECK( length(vb) == length(qb) );

        // Check normalization and angle between unit length vecs/quats
        const vec4<T> nva = normalize(va), nvb = normalize(vb);
        const quat<T> nqa = normalize(qa), nqb = normalize(qb);
        CHECK( nqa.x == nva[0] );
        CHECK( nqa.y == nva[1] );
        CHECK( nqa.z == nva[2] );
        CHECK( nqa.w == nva[3] );
        CHECK( nqb.x == nvb[0] );
        CHECK( nqb.y == nvb[1] );
        CHECK( nqb.z == nvb[2] );
        CHECK( nqb.w == nvb[3] );
        CHECK( uangle(nva,nvb) == uangle(nqa,nqb) );

        // Check lerp, nlerp, and slerp
        const T t = rng;
        const vec4<T> lv = lerp(va, vb, t), nlv = nlerp(va, vb, t), slv = slerp(nva, nvb, t);
        const quat<T> lq = lerp(qa, qb, t), nlq = nlerp(qa, qb, t), slq = slerp(nqa, nqb, t);
        CHECK( lq.x == lv[0] );
        CHECK( lq.y == lv[1] );
        CHECK( lq.z == lv[2] );
        CHECK( lq.w == lv[3] );
        CHECK( nlq.x == nlv[0] );
        CHECK( nlq.y == nlv[1] );
        CHECK( nlq.z == nlv[2] );
        CHECK( nlq.w == nlv[3] );
        CHECK( slq.x == slv[0] );
        CHECK( slq.y == slv[1] );
        CHECK( slq.z == slv[2] );
        CHECK( slq.w == slv[3] );
    }   
}

TEST_CASE_TEMPLATE("Quaternion exponent, logarithm, and power", T, float, double)
{
    constexpr T tau = T(6.2831853071795865), pi = tau/2;
    check_approx_equal( exp(quat<T>{tau,0,0,0}), quat<T>{0,0,0,1} ); // e^tau*i = 1
    check_approx_equal( exp(quat<T>{0,tau,0,0}), quat<T>{0,0,0,1} ); // e^tau*j = 1
    check_approx_equal( exp(quat<T>{0,0,tau,0}), quat<T>{0,0,0,1} ); // e^tau*k = 1
    check_approx_equal( exp(quat<T>{pi,0,0,0}), quat<T>{0,0,0,-1} ); // e^pi*i = -1
    check_approx_equal( exp(quat<T>{0,pi,0,0}), quat<T>{0,0,0,-1} ); // e^pi*j = -1
    check_approx_equal( exp(quat<T>{0,0,pi,0}), quat<T>{0,0,0,-1} ); // e^pi*k = -1

    random_number_generator rng;
    for(int i=0; i<reps; ++i)
    {
        // qexp, qlog, and qpow with a scalar should simply be the exp, log, pow of that scalar
        T s = rng, a = std::abs(s), p = rng;
        check_approx_equal( exp(quat<T>{0,0,0,s}  ), quat<T>{0,0,0,std::exp(s  )} );
        check_approx_equal( log(quat<T>{0,0,0,a}  ), quat<T>{0,0,0,std::log(a  )} );
        check_approx_equal( pow(quat<T>{0,0,0,a},p), quat<T>{0,0,0,std::pow(a,p)} );

        // qpow with integral powers should have the effect of repeated multiplication
        quat<T> q = rng;
        check_approx_equal( exp(log(q)), q );
        check_approx_equal( pow(q, T(2)), q*q );
        check_approx_equal( pow(q, T(3)), q*q*q );
        check_approx_equal( pow(pow(q, T(2)), T(3)), pow(q, T(2*3)) );
        check_approx_equal( pow(q, p+1), q*pow(q, p) );
    }   
}

TEST_CASE_TEMPLATE("qxdir(q) == qvq* where v={1,0,0,0}",  T, int, float, double)
{
    random_number_generator rng;
    for(int i=0; i<reps; ++i)
    {
        const quat<T> q = rng;
        check_approx_equal(quat<T>{qxdir(q),0}, q * quat<T>{1,0,0,0} * conjugate(q));
    }
}

TEST_CASE_TEMPLATE("qzdir(q) == qvq* where v={0,1,0,0}",  T, int, float, double)
{
    random_number_generator rng;
    for(int i=0; i<reps; ++i)
    {
        const quat<T> q = rng;
        check_approx_equal(quat<T>{qydir(q),0}, q * quat<T>{0,1,0,0} * conjugate(q));
    }
}

TEST_CASE_TEMPLATE("qzdir(q) == qvq* where v={0,0,1,0}",  T, int, float, double)
{
    random_number_generator rng;
    for(int i=0; i<reps; ++i)
    {
        const quat<T> q = rng;
        check_approx_equal(quat<T>{qzdir(q),0}, q * quat<T>{0,0,1,0} * conjugate(q));
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
            const vec3<T> v = rng;
            check_approx_equal(quat<T>{qrot(q,v),0}, q * quat<T>{v,0} * conjugate(q));
            check_approx_equal(qrot(q,v), qmat(q)*v);
        }
    }
}
