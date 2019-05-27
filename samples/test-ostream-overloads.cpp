#include "test-linalg.h"
using namespace linalg::ostream_overloads;

#include <sstream>
template<class T> std::string stringify(const T & value) { std::ostringstream ss; ss << value; return ss.str(); }
template<class T> std::wstring wstringify(const T & value) { std::wostringstream ss; ss << value; return ss.str(); }

TEST_CASE("Test ostream << vec")
{
    CHECK(stringify(float2{1,2}) == "{1,2}");
    CHECK(stringify(float3{1,2,3}) == "{1,2,3}");
    CHECK(stringify(float4{1,2,3,4}) == "{1,2,3,4}");
    CHECK(wstringify(float2{1,2}) == L"{1,2}");
    CHECK(wstringify(float3{1,2,3}) == L"{1,2,3}");
    CHECK(wstringify(float4{1,2,3,4}) == L"{1,2,3,4}");
}

TEST_CASE("Test ostream << _scalar")
{
    const float4 v {1,2,3,4};
    CHECK(stringify(v.x) == "1");
    CHECK(stringify(v.y) == "2");
    CHECK(stringify(v.z) == "3");
    CHECK(stringify(v.w) == "4");
    CHECK(wstringify(v.x) == L"1");
    CHECK(wstringify(v.y) == L"2");
    CHECK(wstringify(v.z) == L"3");
    CHECK(wstringify(v.w) == L"4");
}

TEST_CASE("Test ostream << mat")
{
    CHECK(stringify(float2x4{{1,2},{3,4},{5,6},{7,8}}) == "{{1,2},{3,4},{5,6},{7,8}}");
    CHECK(stringify(float3x3{{1,2,3},{4,5,6},{7,8,9}}) == "{{1,2,3},{4,5,6},{7,8,9}}");
    CHECK(stringify(float4x2{{1,2,3,4},{5,6,7,8}}) == "{{1,2,3,4},{5,6,7,8}}");
    CHECK(wstringify(float2x4{{1,2},{3,4},{5,6},{7,8}}) == L"{{1,2},{3,4},{5,6},{7,8}}");
    CHECK(wstringify(float3x3{{1,2,3},{4,5,6},{7,8,9}}) == L"{{1,2,3},{4,5,6},{7,8,9}}");
    CHECK(wstringify(float4x2{{1,2,3,4},{5,6,7,8}}) == L"{{1,2,3,4},{5,6,7,8}}");
}
