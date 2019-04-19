#include "test-linalg.h"

TEST_CASE_TEMPLATE("Quaternion exponent, logarithm, and power", T, float, double)
{
    constexpr T tau = T(6.2831853071795865), pi = tau/2;
    check_approx_equal( qexp(vec4<T>{tau,0,0,0}), vec4<T>{0,0,0,1} ); // e^tau*i = 1
    check_approx_equal( qexp(vec4<T>{0,tau,0,0}), vec4<T>{0,0,0,1} ); // e^tau*j = 1
    check_approx_equal( qexp(vec4<T>{0,0,tau,0}), vec4<T>{0,0,0,1} ); // e^tau*k = 1
    check_approx_equal( qexp(vec4<T>{pi,0,0,0}), vec4<T>{0,0,0,-1} ); // e^pi*i = -1
    check_approx_equal( qexp(vec4<T>{0,pi,0,0}), vec4<T>{0,0,0,-1} ); // e^pi*j = -1
    check_approx_equal( qexp(vec4<T>{0,0,pi,0}), vec4<T>{0,0,0,-1} ); // e^pi*k = -1

    random_number_generator rng;
    for(int i=0; i<reps; ++i)
    {
        // qexp, qlog, and qpow with a scalar should simply be the exp, log, pow of that scalar
        T s = rng, a = std::abs(s), p = rng;
        check_approx_equal( qexp(vec4<T>{0,0,0,s}  ), vec4<T>{0,0,0,std::exp(s  )} );
        check_approx_equal( qlog(vec4<T>{0,0,0,a}  ), vec4<T>{0,0,0,std::log(a  )} );
        check_approx_equal( qpow(vec4<T>{0,0,0,a},p), vec4<T>{0,0,0,std::pow(a,p)} );

        // qpow with integral powers should have the effect of repeated multiplication
        vec4<T> q = rng;
        check_approx_equal( qexp(qlog(q)), q );
        check_approx_equal( qpow(q, T(2)), qmul(q,q) );
        check_approx_equal( qpow(q, T(3)), qmul(q,q,q) );
        check_approx_equal( qpow(q, p+1), qmul(q, qpow(q, p)) );
    } 
}

TEST_CASE_TEMPLATE("qxdir(q) == qvq* where v={1,0,0,0}",  T, int, float, double)
{
    random_number_generator rng;
    for(int i=0; i<reps; ++i)
    {
        const vec4<T> q = rng;
        check_approx_equal(vec4<T>{qxdir(q),0}, qmul(q, vec4<T>{1,0,0,0}, qconj(q)));
    }
}

TEST_CASE_TEMPLATE("qzdir(q) == qvq* where v={0,1,0,0}",  T, int, float, double)
{
    random_number_generator rng;
    for(int i=0; i<reps; ++i)
    {
        const vec4<T> q = rng;
        check_approx_equal(vec4<T>{qydir(q),0}, qmul(q, vec4<T>{0,1,0,0}, qconj(q)));
    }
}

TEST_CASE_TEMPLATE("qzdir(q) == qvq* where v={0,0,1,0}",  T, int, float, double)
{
    random_number_generator rng;
    for(int i=0; i<reps; ++i)
    {
        const vec4<T> q = rng;
        check_approx_equal(vec4<T>{qzdir(q),0}, qmul(q, vec4<T>{0,0,1,0}, qconj(q)));
    }
}

TEST_CASE_TEMPLATE("qrot(q,v) == qvq*",  T, int, float, double)
{
    random_number_generator rng;
    for(int i=0; i<reps; ++i)
    {
        const vec4<T> q = rng;
        for(int j=0; j<reps; ++j)
        {
            const vec3<T> v = rng;
            check_approx_equal(vec4<T>{qrot(q,v),0}, qmul(q, vec4<T>{v,0}, qconj(q)));
            check_approx_equal(qrot(q,v), mul(qmat(q),v));
        }
    }
}
