#include "test-linalg.h"
#include <vector>
#include <tuple>

// Function object which repeatedly appends to a std::vector<T>
struct { template<class T> std::vector<T> operator() (std::vector<T> a, T b) const { a.push_back(b); return a; } } flatten;

TEST_CASE("linalg::fold(...)")
{
    std::vector<float> values {0};

    SUBCASE("vec<T,1>")
    {
        values = fold(flatten, values, float1{1});
        CHECK(values == std::vector<float>{0,1});
    }

    SUBCASE("vec<T,2>")
    {
        values = fold(flatten, values, float2{1,2});
        CHECK(values == std::vector<float>{0,1,2});
    }

    SUBCASE("vec<T,3>")
    {
        values = fold(flatten, values, float3{1,2,3});
        CHECK(values == std::vector<float>{0,1,2,3});
    }

    SUBCASE("vec<T,4>")
    {
        values = fold(flatten, values, float4{1,2,3,4});
        CHECK(values == std::vector<float>{0,1,2,3,4});
    }

    SUBCASE("mat<T,3,2>")
    {
        values = fold(flatten, values, float3x2{{1,2,3},{4,5,6}});
        CHECK(values == std::vector<float>{0,1,2,3,4,5,6});
    }

    SUBCASE("mat<T,3,3>")
    {
        values = fold(flatten, values, float3x3{{1,2,3},{4,5,6},{7,8,9}});
        CHECK(values == std::vector<float>{0,1,2,3,4,5,6,7,8,9});
    }

    SUBCASE("mat<T,2,4>")
    {
        values = fold(flatten, values, float2x4{{1,2},{3,4},{5,6},{7,8}});
        CHECK(values == std::vector<float>{0,1,2,3,4,5,6,7,8});
    }

    SUBCASE("mat<T,1,4>")
    {
        values = fold(flatten, values, float1x4{{1},{3},{5},{7}});
        CHECK(values == std::vector<float>{0,1,3,5,7});
    }
}

// Function object equivalent to the overload set of std::make_tuple(...)
struct { template<class... T> std::tuple<T...> operator() (const T & ... args) const { return std::tuple<T...>(args...); } } tup;

TEST_CASE("linalg::apply(...)")
{
    // We will use the apply(...) function to form tuples from our arguments, then check that the correct arguments were passed in
    const float  fs=1.f, f0=2.f, f1=3.f, f2=4.f, f3=5.f;
    const double ds=1.1, d0=2.2, d1=3.3, d2=4.4, d3=5.5;
    const int    is=111, i0=222, i1=333, i2=444, i3=555;
    const float3  fv {f0,f1,f2}; const float2x2  fm {{f0,f1}, {f2,f3}};
    const double3 dv {d0,d1,d2}; const double2x2 dm {{d0,d1}, {d2,d3}};
    const int3    iv {i0,i1,i2}; const int2x2    im {{i0,i1}, {i2,i3}};

    SUBCASE("Unary application with vectors/matrices/scalars")
    {
        CHECK(apply(tup, fv) == linalg::vec<std::tuple<float>,3>{tup(f0), tup(f1), tup(f2)});
        CHECK(apply(tup, dv) == linalg::vec<std::tuple<double>,3>{tup(d0), tup(d1), tup(d2)});
        CHECK(apply(tup, fm) == linalg::mat<std::tuple<float>,2,2>{{tup(f0), tup(f1)}, {tup(f2), tup(f3)}});
        CHECK(apply(tup, dm) == linalg::mat<std::tuple<double>,2,2>{{tup(d0), tup(d1)}, {tup(d2), tup(d3)}});
        CHECK(linalg::apply(tup, fs) == tup(fs));
        CHECK(linalg::apply(tup, ds) == tup(ds));
    }

    SUBCASE("Binary application with vectors/matrices/scalars")
    {
        CHECK(apply(tup, fv, dv) == linalg::vec<std::tuple<float,double>,3>{tup(f0,d0), tup(f1,d1), tup(f2,d2)});
        CHECK(apply(tup, fv, ds) == linalg::vec<std::tuple<float,double>,3>{tup(f0,ds), tup(f1,ds), tup(f2,ds)});
        CHECK(apply(tup, fs, dv) == linalg::vec<std::tuple<float,double>,3>{tup(fs,d0), tup(fs,d1), tup(fs,d2)});
        CHECK(apply(tup, dv, fv) == linalg::vec<std::tuple<double,float>,3>{tup(d0,f0), tup(d1,f1), tup(d2,f2)});
        CHECK(apply(tup, dv, fs) == linalg::vec<std::tuple<double,float>,3>{tup(d0,fs), tup(d1,fs), tup(d2,fs)});
        CHECK(apply(tup, ds, fv) == linalg::vec<std::tuple<double,float>,3>{tup(ds,f0), tup(ds,f1), tup(ds,f2)});
        CHECK(apply(tup, fm, dm) == linalg::mat<std::tuple<float,double>,2,2>{{tup(f0,d0), tup(f1,d1)}, {tup(f2,d2), tup(f3,d3)}});
        CHECK(apply(tup, fm, ds) == linalg::mat<std::tuple<float,double>,2,2>{{tup(f0,ds), tup(f1,ds)}, {tup(f2,ds), tup(f3,ds)}});
        CHECK(apply(tup, fs, dm) == linalg::mat<std::tuple<float,double>,2,2>{{tup(fs,d0), tup(fs,d1)}, {tup(fs,d2), tup(fs,d3)}});
        CHECK(apply(tup, dm, fm) == linalg::mat<std::tuple<double,float>,2,2>{{tup(d0,f0), tup(d1,f1)}, {tup(d2,f2), tup(d3,f3)}});
        CHECK(apply(tup, dm, fs) == linalg::mat<std::tuple<double,float>,2,2>{{tup(d0,fs), tup(d1,fs)}, {tup(d2,fs), tup(d3,fs)}});
        CHECK(apply(tup, ds, fm) == linalg::mat<std::tuple<double,float>,2,2>{{tup(ds,f0), tup(ds,f1)}, {tup(ds,f2), tup(ds,f3)}});
        CHECK(linalg::apply(tup, fs, ds) == std::tuple<float,double>{tup(fs,ds)});
        CHECK(linalg::apply(tup, ds, fs) == std::tuple<double,float>{tup(ds,fs)});
    }

    SUBCASE("Ternary application with vectors/scalars")
    {
        CHECK(apply(tup, fv, dv, iv) == linalg::vec<std::tuple<float,double,int>,3>{tup(f0,d0,i0), tup(f1,d1,i1), tup(f2,d2,i2)});
        CHECK(apply(tup, fv, dv, is) == linalg::vec<std::tuple<float,double,int>,3>{tup(f0,d0,is), tup(f1,d1,is), tup(f2,d2,is)});
        CHECK(apply(tup, fv, ds, iv) == linalg::vec<std::tuple<float,double,int>,3>{tup(f0,ds,i0), tup(f1,ds,i1), tup(f2,ds,i2)});
        CHECK(apply(tup, fv, ds, is) == linalg::vec<std::tuple<float,double,int>,3>{tup(f0,ds,is), tup(f1,ds,is), tup(f2,ds,is)});
        CHECK(apply(tup, fs, dv, iv) == linalg::vec<std::tuple<float,double,int>,3>{tup(fs,d0,i0), tup(fs,d1,i1), tup(fs,d2,i2)});
        CHECK(apply(tup, fs, dv, is) == linalg::vec<std::tuple<float,double,int>,3>{tup(fs,d0,is), tup(fs,d1,is), tup(fs,d2,is)});
        CHECK(apply(tup, fs, ds, iv) == linalg::vec<std::tuple<float,double,int>,3>{tup(fs,ds,i0), tup(fs,ds,i1), tup(fs,ds,i2)});
        CHECK(apply(tup, iv, fv, dv) == linalg::vec<std::tuple<int,float,double>,3>{tup(i0,f0,d0), tup(i1,f1,d1), tup(i2,f2,d2)});
        CHECK(apply(tup, is, fv, dv) == linalg::vec<std::tuple<int,float,double>,3>{tup(is,f0,d0), tup(is,f1,d1), tup(is,f2,d2)});
        CHECK(apply(tup, iv, fv, ds) == linalg::vec<std::tuple<int,float,double>,3>{tup(i0,f0,ds), tup(i1,f1,ds), tup(i2,f2,ds)});
        CHECK(apply(tup, is, fv, ds) == linalg::vec<std::tuple<int,float,double>,3>{tup(is,f0,ds), tup(is,f1,ds), tup(is,f2,ds)});
        CHECK(apply(tup, iv, fs, dv) == linalg::vec<std::tuple<int,float,double>,3>{tup(i0,fs,d0), tup(i1,fs,d1), tup(i2,fs,d2)});
        CHECK(apply(tup, is, fs, dv) == linalg::vec<std::tuple<int,float,double>,3>{tup(is,fs,d0), tup(is,fs,d1), tup(is,fs,d2)});
        CHECK(apply(tup, iv, fs, ds) == linalg::vec<std::tuple<int,float,double>,3>{tup(i0,fs,ds), tup(i1,fs,ds), tup(i2,fs,ds)});
        CHECK(linalg::apply(tup, fs, ds, is) == std::tuple<float,double,int>{tup(fs,ds,is)});
        CHECK(linalg::apply(tup, is, fs, ds) == std::tuple<int,float,double>{tup(is,fs,ds)});
    }

    SUBCASE("Quaternary application with scalars")
    {
        CHECK(linalg::apply(tup, fs, ds, is, true) == std::tuple<float,double,int,bool>{tup(fs,ds,is,true)});
        CHECK(linalg::apply(tup, true, is, ds, fs) == std::tuple<bool,int,double,float>{tup(true,is,ds,fs)});
    }

    SUBCASE("map(a,f) == apply(f,a)")
    {
        CHECK(map(fv, tup) == linalg::vec<std::tuple<float>,3>{tup(f0), tup(f1), tup(f2)});
        CHECK(map(dv, tup) == linalg::vec<std::tuple<double>,3>{tup(d0), tup(d1), tup(d2)});
        CHECK(map(fm, tup) == linalg::mat<std::tuple<float>,2,2>{{tup(f0), tup(f1)}, {tup(f2), tup(f3)}});
        CHECK(map(dm, tup) == linalg::mat<std::tuple<double>,2,2>{{tup(d0), tup(d1)}, {tup(d2), tup(d3)}});
        CHECK(linalg::map(fs, tup) == tup(fs));
        CHECK(linalg::map(ds, tup) == tup(ds));
    }

    SUBCASE("zip(a,b,f) == apply(f,a,b)")
    {
        CHECK(zip(fv, dv, tup) == linalg::vec<std::tuple<float,double>,3>{tup(f0,d0), tup(f1,d1), tup(f2,d2)});
        CHECK(zip(fv, ds, tup) == linalg::vec<std::tuple<float,double>,3>{tup(f0,ds), tup(f1,ds), tup(f2,ds)});
        CHECK(zip(fs, dv, tup) == linalg::vec<std::tuple<float,double>,3>{tup(fs,d0), tup(fs,d1), tup(fs,d2)});
        CHECK(zip(dv, fv, tup) == linalg::vec<std::tuple<double,float>,3>{tup(d0,f0), tup(d1,f1), tup(d2,f2)});
        CHECK(zip(dv, fs, tup) == linalg::vec<std::tuple<double,float>,3>{tup(d0,fs), tup(d1,fs), tup(d2,fs)});
        CHECK(zip(ds, fv, tup) == linalg::vec<std::tuple<double,float>,3>{tup(ds,f0), tup(ds,f1), tup(ds,f2)});
        CHECK(zip(fm, dm, tup) == linalg::mat<std::tuple<float,double>,2,2>{{tup(f0,d0), tup(f1,d1)}, {tup(f2,d2), tup(f3,d3)}});
        CHECK(zip(fm, ds, tup) == linalg::mat<std::tuple<float,double>,2,2>{{tup(f0,ds), tup(f1,ds)}, {tup(f2,ds), tup(f3,ds)}});
        CHECK(zip(fs, dm, tup) == linalg::mat<std::tuple<float,double>,2,2>{{tup(fs,d0), tup(fs,d1)}, {tup(fs,d2), tup(fs,d3)}});
        CHECK(zip(dm, fm, tup) == linalg::mat<std::tuple<double,float>,2,2>{{tup(d0,f0), tup(d1,f1)}, {tup(d2,f2), tup(d3,f3)}});
        CHECK(zip(dm, fs, tup) == linalg::mat<std::tuple<double,float>,2,2>{{tup(d0,fs), tup(d1,fs)}, {tup(d2,fs), tup(d3,fs)}});
        CHECK(zip(ds, fm, tup) == linalg::mat<std::tuple<double,float>,2,2>{{tup(ds,f0), tup(ds,f1)}, {tup(ds,f2), tup(ds,f3)}});
        CHECK(linalg::zip(fs, ds, tup) == std::tuple<float,double>{tup(fs,ds)});
        CHECK(linalg::zip(ds, fs, tup) == std::tuple<double,float>{tup(ds,fs)});
    }
}