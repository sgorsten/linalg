#include "../linalg.h"
using namespace linalg::aliases;

#include "thirdparty/doctest.h"

// SFINAE helper for determining if T + U is a valid expression
template<class T, class U> static decltype(std::declval<T>() + std::declval<U>(), std::true_type{}) has_operator_plus_helper(int);
template<class T, class U> static std::false_type has_operator_plus_helper(...);
template<class T, class U> static constexpr bool has_operator_plus = decltype(has_operator_plus_helper<T,U>(0))::value;

// SFINAE helper for determining if T * U is a valid expression
template<class T, class U> static decltype(std::declval<T>() * std::declval<U>(), std::true_type{}) has_operator_times_helper(int);
template<class T, class U> static std::false_type has_operator_times_helper(...);
template<class T, class U> static constexpr bool has_operator_times = decltype(has_operator_times_helper<T,U>(0))::value;

TEST_CASE("Well-formed matrix products compile") 
{
    // Enforce deprecation of component-wise matrix multiply
    CHECK_FALSE(has_operator_times<float2x2, float2x2>);
    CHECK_FALSE(has_operator_times<float2x3, float2x3>);
    CHECK_FALSE(has_operator_times<float2x4, float2x4>);
    CHECK_FALSE(has_operator_times<float3x2, float3x2>);
    CHECK_FALSE(has_operator_times<float3x3, float3x3>);
    CHECK_FALSE(has_operator_times<float3x4, float3x4>);
    CHECK_FALSE(has_operator_times<float4x2, float4x2>);
    CHECK_FALSE(has_operator_times<float4x3, float4x3>);
    CHECK_FALSE(has_operator_times<float4x4, float4x4>);
}