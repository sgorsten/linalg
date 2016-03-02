// This file exists to assist in testing that the abstractions used in linalg.h do not interfere with the optimizer's ability
// to inline math code. For now, we are looking to see that there are no function calls to linalg.h functions, and that the
// number of math operations (addss, subss, mulss, divss) does not exceed the amount required to compute the solution using
// scalar arithmetic. In Visual Studio, this can be done by setting a breakpoint, running with the debugger attached, and
// hitting Alt-8. I'd eventually like to automate this process and have it run as part of our test suite.

#include "../linalg.h"
using namespace linalg::aliases;

// A small number of functions, such as matrix multiplication with 4x4 matrices, are SIMD friendly and could conceivably be
// done with fewer operations. The compiler might be smart enough to auto-vectorize. If not, the following trick can be done
// to hint to the compiler which operations can be done with SSE instructions.

#if 0 // Change to 1 to hint to the compiler that it can do ops on float4 (and by proxy, float4xN) using SSE ops
#include <xmmintrin.h>
namespace linalg
{
    float4 operator + (const float4 & a, const float4 & b) { return (const float4 &)_mm_add_ps((const __m128 &)a, (const __m128 &)b); }
    float4 operator - (const float4 & a, const float4 & b) { return (const float4 &)_mm_sub_ps((const __m128 &)a, (const __m128 &)b); }
    float4 operator * (const float4 & a, const float4 & b) { return (const float4 &)_mm_mul_ps((const __m128 &)a, (const __m128 &)b); }
    float4 operator / (const float4 & a, const float4 & b) { return (const float4 &)_mm_div_ps((const __m128 &)a, (const __m128 &)b); }
    float4 operator + (const float4 & a, const float &  b) { return (const float4 &)_mm_add_ps((const __m128 &)a, _mm_set1_ps(b)); }
    float4 operator - (const float4 & a, const float &  b) { return (const float4 &)_mm_sub_ps((const __m128 &)a, _mm_set1_ps(b)); }
    float4 operator * (const float4 & a, const float &  b) { return (const float4 &)_mm_mul_ps((const __m128 &)a, _mm_set1_ps(b)); }
    float4 operator / (const float4 & a, const float &  b) { return (const float4 &)_mm_div_ps((const __m128 &)a, _mm_set1_ps(b)); }
    float4 operator + (const float &  a, const float4 & b) { return (const float4 &)_mm_add_ps(_mm_set1_ps(a), (const __m128 &)b); }
    float4 operator - (const float &  a, const float4 & b) { return (const float4 &)_mm_sub_ps(_mm_set1_ps(a), (const __m128 &)b); }
    float4 operator * (const float &  a, const float4 & b) { return (const float4 &)_mm_mul_ps(_mm_set1_ps(a), (const __m128 &)b); }
    float4 operator / (const float &  a, const float4 & b) { return (const float4 &)_mm_div_ps(_mm_set1_ps(a), (const __m128 &)b); }
}
#endif

volatile float v;

// Release Win32 codegen should contain 3 mulss and 2 addss, no function calls
__declspec(noinline) void test_dot_float3_float3()
{
    v = dot(float3(v,v,v), float3(v,v,v)); 
}

// Release Win32 codegen should contain 9 mulss, 3 divss, 5 addss, 1 subss, and 1 call to __libm_sse2_sqrt_precise
__declspec(noinline) void test_nlerp_float3_float3()
{
    auto r = nlerp(float3(v,v,v), float3(v,v,v), v);
    v = r.x; v = r.y; v = r.z;
}

// Release Win32 codegen should contain no more than 16 mulss and 12 addss, no function calls
// If SSE definitions above are brought into play, should contain no more than 4 mulss and 3 addss, no function calls
__declspec(noinline) void test_mul_float4x4_float4()
{
    auto r = linalg::mul(float4x4({v,v,v,v},{v,v,v,v},{v,v,v,v},{v,v,v,v}), float4(v,v,v,v));
    v = r.x; v = r.y; v = r.z; v = r.w;
}

int main()
{
    test_dot_float3_float3();
    test_nlerp_float3_float3();
    test_mul_float4x4_float4();
    return 0;
}