#include "test-linalg.h"
#include <tuple>

struct make_tuple { template<class... T> std::tuple<T...> operator() (const T & ... args) const { return std::tuple<T...>(args...); } };

TEST_CASE("Test linalg::apply(...)")
{
    // We will use the apply(...) function to form tuples from our arguments, then check that the correct arguments were passed in
    const float  fs=1.f, f0=2.f, f1=3.f, f2=4.f, f3=5.f;
    const double ds=1.1, d0=2.2, d1=3.3, d2=4.4, d3=5.5;
    const int    is=111, i0=222, i1=333, i2=444, i3=555;
    const float2  fv {f0,f1}; const float2x2  fm {{f0,f1}, {f2,f3}}; const quatf fq {f0,f1,f2,f3};
    const double2 dv {d0,d1}; const double2x2 dm {{d0,d1}, {d2,d3}}; const quatd dq {d0,d1,d2,d3};
    const int2    iv {i0,i1}; const int2x2    im {{i0,i1}, {i2,i3}};
    auto tup = make_tuple{};

    // Check unary application with vectors/matrices/quaternions/scalars
    CHECK(linalg::apply(tup, fv) == linalg::vec<std::tuple<float>,2>{tup(f0), tup(f1)});
    CHECK(linalg::apply(tup, fm) == linalg::mat<std::tuple<float>,2,2>{{tup(f0), tup(f1)}, {tup(f2), tup(f3)}});
    CHECK(linalg::apply(tup, fq) == linalg::quat<std::tuple<float>>{tup(f0), tup(f1), tup(f2), tup(f3)});
    CHECK(linalg::apply(tup, fs) == tup(fs));

    // Check binary application with vectors/matrices/quaternions/scalars
    CHECK(linalg::apply(tup, fv, dv) == linalg::vec<std::tuple<float,double>,2>{tup(f0,d0), tup(f1,d1)});
    CHECK(linalg::apply(tup, fv, ds) == linalg::vec<std::tuple<float,double>,2>{tup(f0,ds), tup(f1,ds)});
    CHECK(linalg::apply(tup, fs, dv) == linalg::vec<std::tuple<float,double>,2>{tup(fs,d0), tup(fs,d1)});
    CHECK(linalg::apply(tup, fm, dm) == linalg::mat<std::tuple<float,double>,2,2>{{tup(f0,d0), tup(f1,d1)}, {tup(f2,d2), tup(f3,d3)}});
    CHECK(linalg::apply(tup, fm, ds) == linalg::mat<std::tuple<float,double>,2,2>{{tup(f0,ds), tup(f1,ds)}, {tup(f2,ds), tup(f3,ds)}});
    CHECK(linalg::apply(tup, fs, dm) == linalg::mat<std::tuple<float,double>,2,2>{{tup(fs,d0), tup(fs,d1)}, {tup(fs,d2), tup(fs,d3)}});
    CHECK(linalg::apply(tup, fq, dq) == linalg::quat<std::tuple<float,double>>{tup(f0,d0), tup(f1,d1), tup(f2,d2), tup(f3,d3)});
    CHECK(linalg::apply(tup, fq, ds) == linalg::quat<std::tuple<float,double>>{tup(f0,ds), tup(f1,ds), tup(f2,ds), tup(f3,ds)});
    CHECK(linalg::apply(tup, fs, dq) == linalg::quat<std::tuple<float,double>>{tup(fs,d0), tup(fs,d1), tup(fs,d2), tup(fs,d3)});
    CHECK(linalg::apply(tup, fs, ds) == tup(fs,ds));

    // Check ternary application with vectors/scalars
    CHECK(linalg::apply(tup, fv, dv, iv) == linalg::vec<std::tuple<float,double,int>,2>{tup(f0,d0,i0), tup(f1,d1,i1)});
    CHECK(linalg::apply(tup, fv, dv, is) == linalg::vec<std::tuple<float,double,int>,2>{tup(f0,d0,is), tup(f1,d1,is)});
    CHECK(linalg::apply(tup, fv, ds, iv) == linalg::vec<std::tuple<float,double,int>,2>{tup(f0,ds,i0), tup(f1,ds,i1)});
    CHECK(linalg::apply(tup, fv, ds, is) == linalg::vec<std::tuple<float,double,int>,2>{tup(f0,ds,is), tup(f1,ds,is)});
    CHECK(linalg::apply(tup, fs, dv, iv) == linalg::vec<std::tuple<float,double,int>,2>{tup(fs,d0,i0), tup(fs,d1,i1)});
    CHECK(linalg::apply(tup, fs, dv, is) == linalg::vec<std::tuple<float,double,int>,2>{tup(fs,d0,is), tup(fs,d1,is)});
    CHECK(linalg::apply(tup, fs, ds, iv) == linalg::vec<std::tuple<float,double,int>,2>{tup(fs,ds,i0), tup(fs,ds,i1)});
    CHECK(linalg::apply(tup, fs, ds, is) == tup(fs,ds,is));

    // Check quaternary application with scalars
    CHECK(linalg::apply(tup, fs, ds, is, true) == tup(fs,ds,is,true));

    // Check legacy map(...)
    CHECK(linalg::map(fv, tup) == linalg::vec<std::tuple<float>,2>{tup(f0), tup(f1)});
    CHECK(linalg::map(fm, tup) == linalg::mat<std::tuple<float>,2,2>{{tup(f0), tup(f1)}, {tup(f2), tup(f3)}});
    CHECK(linalg::map(fq, tup) == linalg::quat<std::tuple<float>>{tup(f0), tup(f1), tup(f2), tup(f3)});
    CHECK(linalg::map(fs, tup) == tup(fs));

    // Check legacy zip(...)
    CHECK(linalg::zip(fv, dv, tup) == linalg::vec<std::tuple<float,double>,2>{tup(f0,d0), tup(f1,d1)});
    CHECK(linalg::zip(fv, ds, tup) == linalg::vec<std::tuple<float,double>,2>{tup(f0,ds), tup(f1,ds)});
    CHECK(linalg::zip(fs, dv, tup) == linalg::vec<std::tuple<float,double>,2>{tup(fs,d0), tup(fs,d1)});
    CHECK(linalg::zip(fm, dm, tup) == linalg::mat<std::tuple<float,double>,2,2>{{tup(f0,d0), tup(f1,d1)}, {tup(f2,d2), tup(f3,d3)}});
    CHECK(linalg::zip(fm, ds, tup) == linalg::mat<std::tuple<float,double>,2,2>{{tup(f0,ds), tup(f1,ds)}, {tup(f2,ds), tup(f3,ds)}});
    CHECK(linalg::zip(fs, dm, tup) == linalg::mat<std::tuple<float,double>,2,2>{{tup(fs,d0), tup(fs,d1)}, {tup(fs,d2), tup(fs,d3)}});
    CHECK(linalg::zip(fq, dq, tup) == linalg::quat<std::tuple<float,double>>{tup(f0,d0), tup(f1,d1), tup(f2,d2), tup(f3,d3)});
    CHECK(linalg::zip(fq, ds, tup) == linalg::quat<std::tuple<float,double>>{tup(f0,ds), tup(f1,ds), tup(f2,ds), tup(f3,ds)});
    CHECK(linalg::zip(fs, dq, tup) == linalg::quat<std::tuple<float,double>>{tup(fs,d0), tup(fs,d1), tup(fs,d2), tup(fs,d3)});
    CHECK(linalg::zip(fs, ds, tup) == tup(fs,ds));
}