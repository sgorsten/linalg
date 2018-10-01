#include "../linalg.h"
using namespace linalg::aliases;

#include "thirdparty/doctest.h"

namespace detail
{
    // These function overloads exist if specific the expression T{} op U{} is valid, and return true_type
    template<class T, class U> static decltype(std::declval<T>() + std::declval<U>(), std::true_type{}) has_op_add(int);
    template<class T, class U> static decltype(std::declval<T>() - std::declval<U>(), std::true_type{}) has_op_sub(int);  
    template<class T, class U> static decltype(std::declval<T>() * std::declval<U>(), std::true_type{}) has_op_mul(int);
    template<class T, class U> static decltype(std::declval<T>() / std::declval<U>(), std::true_type{}) has_op_div(int);
    
    // These function overloads always exist, but have lowest selection priority, and return false_type
    template<class T, class U> static std::false_type has_op_add(...);
    template<class T, class U> static std::false_type has_op_sub(...);  
    template<class T, class U> static std::false_type has_op_mul(...);
    template<class T, class U> static std::false_type has_op_div(...);
}

template<class T, class U> static constexpr bool has_op_add = decltype(detail::has_op_add<T,U>(0))::value;
template<class T, class U> static constexpr bool has_op_sub = decltype(detail::has_op_sub<T,U>(0))::value;
template<class T, class U> static constexpr bool has_op_mul = decltype(detail::has_op_mul<T,U>(0))::value;
template<class T, class U> static constexpr bool has_op_div = decltype(detail::has_op_div<T,U>(0))::value;

TEST_CASE("Sums between matrices of the same dimension are well formed")
{
    CHECK(has_op_add<float2x2, float2x2>);
    CHECK(has_op_add<float2x3, float2x3>);
    CHECK(has_op_add<float2x4, float2x4>);
    CHECK(has_op_add<float3x2, float3x2>);
    CHECK(has_op_add<float3x3, float3x3>);
    CHECK(has_op_add<float3x4, float3x4>);
    CHECK(has_op_add<float4x2, float4x2>);
    CHECK(has_op_add<float4x3, float4x3>);
    CHECK(has_op_add<float4x4, float4x4>);
}

TEST_CASE("Differences between matrices of the same dimension are well formed")
{
    CHECK(has_op_sub<float2x2, float2x2>);
    CHECK(has_op_sub<float2x3, float2x3>);
    CHECK(has_op_sub<float2x4, float2x4>);
    CHECK(has_op_sub<float3x2, float3x2>);
    CHECK(has_op_sub<float3x3, float3x3>);
    CHECK(has_op_sub<float3x4, float3x4>);
    CHECK(has_op_sub<float4x2, float4x2>);
    CHECK(has_op_sub<float4x3, float4x3>);
    CHECK(has_op_sub<float4x4, float4x4>);
}

TEST_CASE("Matrices may be multiplied by a scalar on the right or on the left")
{
    CHECK(has_op_mul<float2x2, float>);
    CHECK(has_op_mul<float2x3, float>);
    CHECK(has_op_mul<float2x4, float>);
    CHECK(has_op_mul<float3x2, float>);
    CHECK(has_op_mul<float3x3, float>);
    CHECK(has_op_mul<float3x4, float>);
    CHECK(has_op_mul<float4x2, float>);
    CHECK(has_op_mul<float4x3, float>);
    CHECK(has_op_mul<float4x4, float>);

    CHECK(has_op_mul<float, float2x2>);
    CHECK(has_op_mul<float, float2x3>);
    CHECK(has_op_mul<float, float2x4>);
    CHECK(has_op_mul<float, float3x2>);
    CHECK(has_op_mul<float, float3x3>);
    CHECK(has_op_mul<float, float3x4>);
    CHECK(has_op_mul<float, float4x2>);
    CHECK(has_op_mul<float, float4x3>);
    CHECK(has_op_mul<float, float4x4>);
}

TEST_CASE("Matrices may be divided by a scalar")
{
    CHECK(has_op_div<float2x2, float>);
    CHECK(has_op_div<float2x3, float>);
    CHECK(has_op_div<float2x4, float>);
    CHECK(has_op_div<float3x2, float>);
    CHECK(has_op_div<float3x3, float>);
    CHECK(has_op_div<float3x4, float>);
    CHECK(has_op_div<float4x2, float>);
    CHECK(has_op_div<float4x3, float>);
    CHECK(has_op_div<float4x4, float>);
}

TEST_CASE("Scalars may not be divided by a matrix")
{
    CHECK_FALSE(has_op_div<float, float2x2>);
    CHECK_FALSE(has_op_div<float, float2x3>);
    CHECK_FALSE(has_op_div<float, float2x4>);
    CHECK_FALSE(has_op_div<float, float3x2>);
    CHECK_FALSE(has_op_div<float, float3x3>);
    CHECK_FALSE(has_op_div<float, float3x4>);
    CHECK_FALSE(has_op_div<float, float4x2>);
    CHECK_FALSE(has_op_div<float, float4x3>);
    CHECK_FALSE(has_op_div<float, float4x4>);
}

TEST_CASE("Sums between objects of the same size are well formed")
{
    CHECK(has_op_add<float2x2, float2x2>);
    CHECK(has_op_add<float2x3, float2x3>);
    CHECK(has_op_add<float2x4, float2x4>);
    CHECK(has_op_add<float3x2, float3x2>);
    CHECK(has_op_add<float3x3, float3x3>);
    CHECK(has_op_add<float3x4, float3x4>);
    CHECK(has_op_add<float4x2, float4x2>);
    CHECK(has_op_add<float4x3, float4x3>);
    CHECK(has_op_add<float4x4, float4x4>);

    CHECK(has_op_sub<float2x2, float2x2>);
    CHECK(has_op_sub<float2x3, float2x3>);
    CHECK(has_op_sub<float2x4, float2x4>);
    CHECK(has_op_sub<float3x2, float3x2>);
    CHECK(has_op_sub<float3x3, float3x3>);
    CHECK(has_op_sub<float3x4, float3x4>);
    CHECK(has_op_sub<float4x2, float4x2>);
    CHECK(has_op_sub<float4x3, float4x3>);
    CHECK(has_op_sub<float4x4, float4x4>);
}

TEST_CASE("Matrix products are well-formed if the left-hand number of columns matches the right-hand number of rows") 
{
    CHECK(has_op_mul<float2x2, float2x2>);
    CHECK(has_op_mul<float2x2, float2x3>);
    CHECK(has_op_mul<float2x2, float2x4>);
    CHECK(has_op_mul<float2x3, float3x2>);
    CHECK(has_op_mul<float2x3, float3x3>);
    CHECK(has_op_mul<float2x3, float3x4>);
    CHECK(has_op_mul<float2x4, float4x2>);
    CHECK(has_op_mul<float2x4, float4x3>);
    CHECK(has_op_mul<float2x4, float4x4>);

    CHECK(has_op_mul<float3x2, float2x2>);
    CHECK(has_op_mul<float3x2, float2x3>);
    CHECK(has_op_mul<float3x2, float2x4>);
    CHECK(has_op_mul<float3x3, float3x2>);
    CHECK(has_op_mul<float3x3, float3x3>);
    CHECK(has_op_mul<float3x3, float3x4>);
    CHECK(has_op_mul<float3x4, float4x2>);
    CHECK(has_op_mul<float3x4, float4x3>);
    CHECK(has_op_mul<float3x4, float4x4>);

    CHECK(has_op_mul<float4x2, float2x2>);
    CHECK(has_op_mul<float4x2, float2x3>);
    CHECK(has_op_mul<float4x2, float2x4>);
    CHECK(has_op_mul<float4x3, float3x2>);
    CHECK(has_op_mul<float4x3, float3x3>);
    CHECK(has_op_mul<float4x3, float3x4>);
    CHECK(has_op_mul<float4x4, float4x2>);
    CHECK(has_op_mul<float4x4, float4x3>);
    CHECK(has_op_mul<float4x4, float4x4>);
}

TEST_CASE("Matrix products are ill-formed if the left-hand number of columns does not match the right-hand number of rows") 
{
    CHECK_FALSE(has_op_mul<float2x2, float3x2>);
    CHECK_FALSE(has_op_mul<float2x2, float3x3>);
    CHECK_FALSE(has_op_mul<float2x2, float3x4>);
    CHECK_FALSE(has_op_mul<float2x2, float4x2>);
    CHECK_FALSE(has_op_mul<float2x2, float4x3>);
    CHECK_FALSE(has_op_mul<float2x2, float4x4>);
    CHECK_FALSE(has_op_mul<float2x3, float2x2>);
    CHECK_FALSE(has_op_mul<float2x3, float2x3>);
    CHECK_FALSE(has_op_mul<float2x3, float2x4>);
    CHECK_FALSE(has_op_mul<float2x3, float4x2>);
    CHECK_FALSE(has_op_mul<float2x3, float4x3>);
    CHECK_FALSE(has_op_mul<float2x3, float4x4>);
    CHECK_FALSE(has_op_mul<float2x4, float2x2>);
    CHECK_FALSE(has_op_mul<float2x4, float2x3>);
    CHECK_FALSE(has_op_mul<float2x4, float2x4>);
    CHECK_FALSE(has_op_mul<float2x4, float3x2>);
    CHECK_FALSE(has_op_mul<float2x4, float3x3>);
    CHECK_FALSE(has_op_mul<float2x4, float3x4>);

    CHECK_FALSE(has_op_mul<float3x2, float3x2>);
    CHECK_FALSE(has_op_mul<float3x2, float3x3>);
    CHECK_FALSE(has_op_mul<float3x2, float3x4>);
    CHECK_FALSE(has_op_mul<float3x2, float4x2>);
    CHECK_FALSE(has_op_mul<float3x2, float4x3>);
    CHECK_FALSE(has_op_mul<float3x2, float4x4>);
    CHECK_FALSE(has_op_mul<float3x3, float2x2>);
    CHECK_FALSE(has_op_mul<float3x3, float2x3>);
    CHECK_FALSE(has_op_mul<float3x3, float2x4>);
    CHECK_FALSE(has_op_mul<float3x3, float4x2>);
    CHECK_FALSE(has_op_mul<float3x3, float4x3>);
    CHECK_FALSE(has_op_mul<float3x3, float4x4>);
    CHECK_FALSE(has_op_mul<float3x4, float2x2>);
    CHECK_FALSE(has_op_mul<float3x4, float2x3>);
    CHECK_FALSE(has_op_mul<float3x4, float2x4>);
    CHECK_FALSE(has_op_mul<float3x4, float3x2>);
    CHECK_FALSE(has_op_mul<float3x4, float3x3>);
    CHECK_FALSE(has_op_mul<float3x4, float3x4>);

    CHECK_FALSE(has_op_mul<float4x2, float3x2>);
    CHECK_FALSE(has_op_mul<float4x2, float3x3>);
    CHECK_FALSE(has_op_mul<float4x2, float3x4>);
    CHECK_FALSE(has_op_mul<float4x2, float4x2>);
    CHECK_FALSE(has_op_mul<float4x2, float4x3>);
    CHECK_FALSE(has_op_mul<float4x2, float4x4>);
    CHECK_FALSE(has_op_mul<float4x3, float2x2>);
    CHECK_FALSE(has_op_mul<float4x3, float2x3>);
    CHECK_FALSE(has_op_mul<float4x3, float2x4>);
    CHECK_FALSE(has_op_mul<float4x3, float4x2>);
    CHECK_FALSE(has_op_mul<float4x3, float4x3>);
    CHECK_FALSE(has_op_mul<float4x3, float4x4>);
    CHECK_FALSE(has_op_mul<float4x4, float2x2>);
    CHECK_FALSE(has_op_mul<float4x4, float2x3>);
    CHECK_FALSE(has_op_mul<float4x4, float2x4>);
    CHECK_FALSE(has_op_mul<float4x4, float3x2>);
    CHECK_FALSE(has_op_mul<float4x4, float3x3>);
    CHECK_FALSE(has_op_mul<float4x4, float3x4>);
}