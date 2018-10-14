#include "test-linalg.h"

TEST_CASE("Test linalg::apply(...)")
{
    // We will use the apply(...) function to form tuples from our arguments, then check that the correct arguments were passed in
    const float  fs=1.f, f0=2.f, f1=3.f, f2=4.f, f3=5.f;
    const double ds=1.1, d0=2.2, d1=3.3, d2=4.4, d3=5.5;
    const int    is=111, i0=222, i1=333, i2=444, i3=555;
    const float2  fv {f0,f1}; const float2x2  fm {{f0,f1},{f2,f3}}; const quatf fq {f0,f1,f2,f3};
    const double2 dv {d0,d1}; const double2x2 dm {{d0,d1},{d2,d3}}; const quatd dq {d0,d1,d2,d3};
    const int2    iv {i0,i1}; const int2x2    im {{i0,i1},{i2,i3}};
    const auto make_tuple = [](const auto & ... args) { return std::make_tuple(args...); };

    // Check unary application with vectors/matrices/quaternions/scalars
    CHECK(linalg::apply(make_tuple, fv) == linalg::vec<std::tuple<float>,2>{{f0},{f1}});
    CHECK(linalg::apply(make_tuple, fm) == linalg::mat<std::tuple<float>,2,2>{{{f0},{f1}},{{f2},{f3}}});
    CHECK(linalg::apply(make_tuple, fq) == linalg::quat<std::tuple<float>>{{f0},{f1},{f2},{f3}});
    CHECK(linalg::apply(make_tuple, fs) == std::tuple<float>{fs});

    // Check binary application with vectors/matrices/quaternions/scalars
    CHECK(linalg::apply(make_tuple, fv, dv) == linalg::vec<std::tuple<float,double>,2>{{f0,d0},{f1,d1}});
    CHECK(linalg::apply(make_tuple, fv, ds) == linalg::vec<std::tuple<float,double>,2>{{f0,ds},{f1,ds}});
    CHECK(linalg::apply(make_tuple, fs, dv) == linalg::vec<std::tuple<float,double>,2>{{fs,d0},{fs,d1}});
    CHECK(linalg::apply(make_tuple, fm, dm) == linalg::mat<std::tuple<float,double>,2,2>{{{f0,d0},{f1,d1}},{{f2,d2},{f3,d3}}});
    CHECK(linalg::apply(make_tuple, fm, ds) == linalg::mat<std::tuple<float,double>,2,2>{{{f0,ds},{f1,ds}},{{f2,ds},{f3,ds}}});
    CHECK(linalg::apply(make_tuple, fs, dm) == linalg::mat<std::tuple<float,double>,2,2>{{{fs,d0},{fs,d1}},{{fs,d2},{fs,d3}}});
    CHECK(linalg::apply(make_tuple, fq, dq) == linalg::quat<std::tuple<float,double>>{{f0,d0},{f1,d1},{f2,d2},{f3,d3}});
    CHECK(linalg::apply(make_tuple, fq, ds) == linalg::quat<std::tuple<float,double>>{{f0,ds},{f1,ds},{f2,ds},{f3,ds}});
    CHECK(linalg::apply(make_tuple, fs, dq) == linalg::quat<std::tuple<float,double>>{{fs,d0},{fs,d1},{fs,d2},{fs,d3}});
    CHECK(linalg::apply(make_tuple, fs, ds) == std::tuple<float,double>{fs,ds});

    // Check ternary application with vectors/scalars
    CHECK(linalg::apply(make_tuple, fv, dv, iv) == linalg::vec<std::tuple<float,double,int>,2>{{f0,d0,i0},{f1,d1,i1}});
    CHECK(linalg::apply(make_tuple, fv, dv, is) == linalg::vec<std::tuple<float,double,int>,2>{{f0,d0,is},{f1,d1,is}});
    CHECK(linalg::apply(make_tuple, fv, ds, iv) == linalg::vec<std::tuple<float,double,int>,2>{{f0,ds,i0},{f1,ds,i1}});
    CHECK(linalg::apply(make_tuple, fv, ds, is) == linalg::vec<std::tuple<float,double,int>,2>{{f0,ds,is},{f1,ds,is}});
    CHECK(linalg::apply(make_tuple, fs, dv, iv) == linalg::vec<std::tuple<float,double,int>,2>{{fs,d0,i0},{fs,d1,i1}});
    CHECK(linalg::apply(make_tuple, fs, dv, is) == linalg::vec<std::tuple<float,double,int>,2>{{fs,d0,is},{fs,d1,is}});
    CHECK(linalg::apply(make_tuple, fs, ds, iv) == linalg::vec<std::tuple<float,double,int>,2>{{fs,ds,i0},{fs,ds,i1}});
    CHECK(linalg::apply(make_tuple, fs, ds, is) == std::tuple<float,double,int>{fs,ds,is});

    // Check quaternary application with scalars
    CHECK(linalg::apply(make_tuple, fs, ds, is, true) == std::tuple<float,double,int,bool>{fs,ds,is,true});

    // Check legacy map(...)
    CHECK(linalg::map(fv, make_tuple) == linalg::vec<std::tuple<float>,2>{{f0},{f1}});
    CHECK(linalg::map(fm, make_tuple) == linalg::mat<std::tuple<float>,2,2>{{{f0},{f1}},{{f2},{f3}}});
    CHECK(linalg::map(fq, make_tuple) == linalg::quat<std::tuple<float>>{{f0},{f1},{f2},{f3}});
    CHECK(linalg::map(fs, make_tuple) == std::tuple<float>{fs});

    // Check legacy zip(...)
    CHECK(linalg::zip(fv, dv, make_tuple) == linalg::vec<std::tuple<float,double>,2>{{f0,d0},{f1,d1}});
    CHECK(linalg::zip(fv, ds, make_tuple) == linalg::vec<std::tuple<float,double>,2>{{f0,ds},{f1,ds}});
    CHECK(linalg::zip(fs, dv, make_tuple) == linalg::vec<std::tuple<float,double>,2>{{fs,d0},{fs,d1}});
    CHECK(linalg::zip(fm, dm, make_tuple) == linalg::mat<std::tuple<float,double>,2,2>{{{f0,d0},{f1,d1}},{{f2,d2},{f3,d3}}});
    CHECK(linalg::zip(fm, ds, make_tuple) == linalg::mat<std::tuple<float,double>,2,2>{{{f0,ds},{f1,ds}},{{f2,ds},{f3,ds}}});
    CHECK(linalg::zip(fs, dm, make_tuple) == linalg::mat<std::tuple<float,double>,2,2>{{{fs,d0},{fs,d1}},{{fs,d2},{fs,d3}}});
    CHECK(linalg::zip(fq, dq, make_tuple) == linalg::quat<std::tuple<float,double>>{{f0,d0},{f1,d1},{f2,d2},{f3,d3}});
    CHECK(linalg::zip(fq, ds, make_tuple) == linalg::quat<std::tuple<float,double>>{{f0,ds},{f1,ds},{f2,ds},{f3,ds}});
    CHECK(linalg::zip(fs, dq, make_tuple) == linalg::quat<std::tuple<float,double>>{{fs,d0},{fs,d1},{fs,d2},{fs,d3}});
    CHECK(linalg::zip(fs, ds, make_tuple) == std::tuple<float,double>{fs,ds});
}