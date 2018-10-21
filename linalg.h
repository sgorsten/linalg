// linalg.h - v3.0 alpha - Single-header public domain linear algebra library
// This is a prerelease version and should be considered unstable.
//
// The intent of this library is to provide the bulk of the functionality
// you need to write programs that frequently use small, fixed-size vectors
// and matrices, in domains such as computational geometry or computer
// graphics. It strives for terse, readable source code.
//
// The original author of this software is Sterling Orsten, and its permanent
// home is <http://github.com/sgorsten/linalg/>. If you find this software
// useful, an acknowledgement in your source text and/or product documentation
// is appreciated, but not required.
//
// The author acknowledges significant insights and contributions by:
//     Stan Melax <http://github.com/melax/>
//     Dimitri Diakopoulos <http://github.com/ddiakopoulos/>



// This is free and unencumbered software released into the public domain.
// 
// Anyone is free to copy, modify, publish, use, compile, sell, or
// distribute this software, either in source code form or as a compiled
// binary, for any purpose, commercial or non-commercial, and by any
// means.
// 
// In jurisdictions that recognize copyright laws, the author or authors
// of this software dedicate any and all copyright interest in the
// software to the public domain. We make this dedication for the benefit
// of the public at large and to the detriment of our heirs and
// successors. We intend this dedication to be an overt act of
// relinquishment in perpetuity of all present and future rights to this
// software under copyright law.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
// IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
// OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
// ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
// OTHER DEALINGS IN THE SOFTWARE.
// 
// For more information, please refer to <http://unlicense.org/>



#pragma once
#ifndef LINALG_H
#define LINALG_H

#include <cmath>        // For std::sqrt, std::sin, std::cos, etc.
#include <cstdlib>      // For std::abs
#include <cstddef>      // For std::nullptr_t
#include <cstdint>      // For std::uint8_t, std::uint16_t, std::int16_t, etc.
#include <type_traits>  // For std::is_arithmetic, std::is_same, std::enable_if, std::conditional, std::remove_reference
#include <functional>   // For std::hash
#include <iosfwd>       // For std::basic_ostream

// In Visual Studio 2015, `constexpr` applied to a member function implies `const`, which causes ambiguous overload resolution
#if _MSC_VER <= 1900
#define LINALG_CONSTEXPR14
#else
#define LINALG_CONSTEXPR14 constexpr
#endif

namespace linalg
{
    //////////////////////
    // Type definitions //
    //////////////////////

    // Small, fixed-length vector type, consisting of exactly M elements of type T, and presumed to be a column-vector unless otherwise noted.
    template<class T, int M> struct vec;

    // Small, fixed-size matrix type, consisting of exactly M rows and N columns of type T, stored in column-major order.
    template<class T, int M, int N> struct mat;

    // Quaternion type, consisting of four elements of type T, representing the quaternion xi + yj + zk + w.
    template<class T> struct quat;

    // Specialize converter<T,U> with a function application operator that converts type U to type T to enable implicit conversions
    template<class T, class U> struct converter {};

    // Define a type which will convert to the multiplicative identity of any given algebraic object
    struct identity_t { constexpr explicit identity_t(int) {} };
    template<class T> struct converter<mat<T,1,1>, identity_t> { mat<T,1,1> operator() (identity_t) const { return {vec<T,1>{1}}; } };
    template<class T> struct converter<mat<T,2,2>, identity_t> { mat<T,2,2> operator() (identity_t) const { return {{1,0},{0,1}}; } };
    template<class T> struct converter<mat<T,3,3>, identity_t> { mat<T,3,3> operator() (identity_t) const { return {{1,0,0},{0,1,0},{0,0,1}}; } };
    template<class T> struct converter<mat<T,4,4>, identity_t> { mat<T,4,4> operator() (identity_t) const { return {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}}; } };
    template<class T> struct converter<quat<T>, identity_t> { quat<T> operator() (identity_t) const { return {0,0,0,1}; } };
    constexpr identity_t identity {1};

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Implementation details. Do not make use of the contents of this namespace from outside the library. //
    /////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Support for named scalar types, as in https://t0rakka.silvrback.com/simd-scalar-accessor, cannot live in detail namespace for ADL reasons, but should not be used by user code
    template<class T, class A, int I> struct _scalar;
    template<class T, class A, int... I> struct _lswizzle;
    template<class T, class A, int... I> struct _rswizzle;

    namespace detail
    {
        template<class A> struct element_storage { A elems; };

        // Support for user defined conversions
        template<class T, class U> using conv_t = typename std::enable_if<!std::is_same<T,U>::value, decltype(converter<T,U>()(std::declval<U>()))>::type;
        template<class T, class U> constexpr T convert(const U & u) { return converter<T,U>()(u); }

        // Helper which reduces _scalar, _lswizzleN, _rswizzleN types to their intended targets
        template<class T> struct unpack { using type=T; };
        template<class T, class A, int I> struct unpack<_scalar<T,A,I>> { using type=T; };
        template<class T, class A, int... I> struct unpack<_lswizzle<T,A,I...>> { using type=vec<T,sizeof...(I)>; };
        template<class T, class A, int... I> struct unpack<_rswizzle<T,A,I...>> { using type=vec<T,sizeof...(I)>; };
        template<class T> using unpack_t = typename unpack<T>::type;

        // Type returned by the compare(...) function which supports all six comparison operators against 0
        template<class T> struct ord { T a,b; };
        template<class T> constexpr bool operator == (const ord<T> & o, std::nullptr_t) { return o.a == o.b; }
        template<class T> constexpr bool operator != (const ord<T> & o, std::nullptr_t) { return !(o.a == o.b); }
        template<class T> constexpr bool operator < (const ord<T> & o, std::nullptr_t) { return o.a < o.b; }
        template<class T> constexpr bool operator > (const ord<T> & o, std::nullptr_t) { return o.b < o.a; }
        template<class T> constexpr bool operator <= (const ord<T> & o, std::nullptr_t) { return !(o.b < o.a); }
        template<class T> constexpr bool operator >= (const ord<T> & o, std::nullptr_t) { return !(o.a < o.b); }

        // Patterns which can be used with the compare(...) function
        template<class A, class B> struct any_compare {};
        template<class T> struct any_compare<vec<T,1>,vec<T,1>> { using type=ord<T>; constexpr ord<T> operator() (const vec<T,1> & a, const vec<T,1> & b) const { return ord<T>{a[0],b[0]}; } };
        template<class T> struct any_compare<vec<T,2>,vec<T,2>> { using type=ord<T>; constexpr ord<T> operator() (const vec<T,2> & a, const vec<T,2> & b) const { return !(a[0]==b[0]) ? ord<T>{a[0],b[0]} : ord<T>{a[1],b[1]}; } };
        template<class T> struct any_compare<vec<T,3>,vec<T,3>> { using type=ord<T>; constexpr ord<T> operator() (const vec<T,3> & a, const vec<T,3> & b) const { return !(a[0]==b[0]) ? ord<T>{a[0],b[0]} : !(a[1]==b[1]) ? ord<T>{a[1],b[1]} : ord<T>{a[2],b[2]}; } };
        template<class T> struct any_compare<vec<T,4>,vec<T,4>> { using type=ord<T>; constexpr ord<T> operator() (const vec<T,4> & a, const vec<T,4> & b) const { return !(a[0]==b[0]) ? ord<T>{a[0],b[0]} : !(a[1]==b[1]) ? ord<T>{a[1],b[1]} : !(a[2]==b[2]) ? ord<T>{a[2],b[2]} : ord<T>{a[3],b[3]}; } };
        template<class T, int M> struct any_compare<mat<T,M,1>,mat<T,M,1>> { using type=ord<T>; constexpr ord<T> operator() (const mat<T,M,1> & a, const mat<T,M,1> & b) const { return compare(a[0],b[0]); } };
        template<class T, int M> struct any_compare<mat<T,M,2>,mat<T,M,2>> { using type=ord<T>; constexpr ord<T> operator() (const mat<T,M,2> & a, const mat<T,M,2> & b) const { return a[0]!=b[0] ? compare(a[0],b[0]) : compare(a[1],b[1]); } };
        template<class T, int M> struct any_compare<mat<T,M,3>,mat<T,M,3>> { using type=ord<T>; constexpr ord<T> operator() (const mat<T,M,3> & a, const mat<T,M,3> & b) const { return a[0]!=b[0] ? compare(a[0],b[0]) : a[1]!=b[1] ? compare(a[1],b[1]) : compare(a[2],b[2]); } };
        template<class T, int M> struct any_compare<mat<T,M,4>,mat<T,M,4>> { using type=ord<T>; constexpr ord<T> operator() (const mat<T,M,4> & a, const mat<T,M,4> & b) const { return a[0]!=b[0] ? compare(a[0],b[0]) : a[1]!=b[1] ? compare(a[1],b[1]) : a[2]!=b[2] ? compare(a[2],b[2]) : compare(a[3],b[3]); } };
        template<class T> struct any_compare<quat<T>,quat<T>> { using type=ord<T>; constexpr ord<T> operator() (const quat<T> & a, const quat<T> & b) const { return !(a.x==b.x) ? ord<T>{a.x,b.x} : !(a.y==b.y) ? ord<T>{a.y,b.y} : !(a.z==b.z) ? ord<T>{a.z,b.z} : ord<T>{a.w,b.w}; } };

        // Stand-in for std::integer_sequence/std::make_integer_sequence
        template<int... I> struct seq {};
        template<int N> struct make_seq_impl;
        template<> struct make_seq_impl<0> { using type=seq<>; };
        template<> struct make_seq_impl<1> { using type=seq<0>; };
        template<> struct make_seq_impl<2> { using type=seq<0,1>; };
        template<> struct make_seq_impl<3> { using type=seq<0,1,2>; };
        template<> struct make_seq_impl<4> { using type=seq<0,1,2,3>; };
        template<int N> using make_seq = typename make_seq_impl<N>::type;

        // SFINAE helpers to determine result of function application
        template<class F, class... T> using ret_t = decltype(std::declval<F>()(std::declval<T>()...));

        // SFINAE helper which is defined if all provided types are scalars
        struct empty {};
        template<class... T> struct scalars;
        template<> struct scalars<> { using type=void; };
        template<class T, class... U> struct scalars<T,U...> : std::conditional<std::is_arithmetic<T>::value, scalars<U...>, empty>::type {};
        template<class... T> using scalars_t = typename scalars<T...>::type;

        // Define various patterns for apply(...)
        template<class F, class Void, class... T> struct vec_apply {}; // Patterns which contain only vectors or scalars
        template<class F, class Void, class... T> struct axa_apply {}; // Patterns of the form: algebraic, algebraic
        template<class F, class Void, class... T> struct axs_apply {}; // Patterns of the form: algebraic, scalar
        template<class F, class Void, class... T> struct sxa_apply {}; // Patterns of the form: scalar, algebraic

        template<class F, int M, class A                  > struct vec_apply<F, scalars_t<   >, vec<A,M>                    > { using type=vec<ret_t<F,A    >,M>; enum {size=M}; template<int... I> static constexpr type impl(seq<I...>, F f, const vec<A,M> & a                                        ) { return {f(a[I]            )...}; } };
        template<class F, int M, class A, class B         > struct vec_apply<F, scalars_t<   >, vec<A,M>, vec<B,M>          > { using type=vec<ret_t<F,A,B  >,M>; enum {size=M}; template<int... I> static constexpr type impl(seq<I...>, F f, const vec<A,M> & a, const vec<B,M> & b                    ) { return {f(a[I], b[I]      )...}; } };
        template<class F, int M, class A, class B         > struct vec_apply<F, scalars_t<B  >, vec<A,M>, B                 > { using type=vec<ret_t<F,A,B  >,M>; enum {size=M}; template<int... I> static constexpr type impl(seq<I...>, F f, const vec<A,M> & a, B                b                    ) { return {f(a[I], b         )...}; } };
        template<class F, int M, class A, class B         > struct vec_apply<F, scalars_t<A  >, A,        vec<B,M>          > { using type=vec<ret_t<F,A,B  >,M>; enum {size=M}; template<int... I> static constexpr type impl(seq<I...>, F f, A                a, const vec<B,M> & b                    ) { return {f(a,    b[I]      )...}; } };
        template<class F, int M, class A, class B, class C> struct vec_apply<F, scalars_t<   >, vec<A,M>, vec<B,M>, vec<C,M>> { using type=vec<ret_t<F,A,B,C>,M>; enum {size=M}; template<int... I> static constexpr type impl(seq<I...>, F f, const vec<A,M> & a, const vec<B,M> & b, const vec<C,M> & c) { return {f(a[I], b[I], c[I])...}; } };
        template<class F, int M, class A, class B, class C> struct vec_apply<F, scalars_t<C  >, vec<A,M>, vec<B,M>, C       > { using type=vec<ret_t<F,A,B,C>,M>; enum {size=M}; template<int... I> static constexpr type impl(seq<I...>, F f, const vec<A,M> & a, const vec<B,M> & b, C                c) { return {f(a[I], b[I], c   )...}; } };
        template<class F, int M, class A, class B, class C> struct vec_apply<F, scalars_t<B  >, vec<A,M>, B,        vec<C,M>> { using type=vec<ret_t<F,A,B,C>,M>; enum {size=M}; template<int... I> static constexpr type impl(seq<I...>, F f, const vec<A,M> & a, B                b, const vec<C,M> & c) { return {f(a[I], b,    c[I])...}; } };
        template<class F, int M, class A, class B, class C> struct vec_apply<F, scalars_t<B,C>, vec<A,M>, B,        C       > { using type=vec<ret_t<F,A,B,C>,M>; enum {size=M}; template<int... I> static constexpr type impl(seq<I...>, F f, const vec<A,M> & a, B                b, C                c) { return {f(a[I], b,    c   )...}; } };
        template<class F, int M, class A, class B, class C> struct vec_apply<F, scalars_t<A  >, A,        vec<B,M>, vec<C,M>> { using type=vec<ret_t<F,A,B,C>,M>; enum {size=M}; template<int... I> static constexpr type impl(seq<I...>, F f, A                a, const vec<B,M> & b, const vec<C,M> & c) { return {f(a,    b[I], c[I])...}; } };
        template<class F, int M, class A, class B, class C> struct vec_apply<F, scalars_t<A,C>, A,        vec<B,M>, C       > { using type=vec<ret_t<F,A,B,C>,M>; enum {size=M}; template<int... I> static constexpr type impl(seq<I...>, F f, A                a, const vec<B,M> & b, C                c) { return {f(a,    b[I], c   )...}; } };
        template<class F, int M, class A, class B, class C> struct vec_apply<F, scalars_t<A,B>, A,        B,        vec<C,M>> { using type=vec<ret_t<F,A,B,C>,M>; enum {size=M}; template<int... I> static constexpr type impl(seq<I...>, F f, A                a, B                b, const vec<C,M> & c) { return {f(a,    b,    c[I])...}; } };

        template<class F, int M, int N, class A         > struct axa_apply<F, scalars_t< >, mat<A,M,N>            > { using type=mat<ret_t<F,A  >,M,N>; enum {size=N}; template<int... J> static constexpr type impl(seq<J...>, F f, const mat<A,M,N> & a                      ) { return {vec_apply<F, void, vec<A,M>          >::impl(make_seq<M>{}, f, a[J]      )...}; } };
        template<class F, int M, int N, class A, class B> struct axa_apply<F, scalars_t< >, mat<A,M,N>, mat<B,M,N>> { using type=mat<ret_t<F,A,B>,M,N>; enum {size=N}; template<int... J> static constexpr type impl(seq<J...>, F f, const mat<A,M,N> & a, const mat<B,M,N> & b) { return {vec_apply<F, void, vec<A,M>, vec<B,M>>::impl(make_seq<M>{}, f, a[J], b[J])...}; } };
        template<class F, int M, int N, class A, class B> struct axs_apply<F, scalars_t<B>, mat<A,M,N>, B         > { using type=mat<ret_t<F,A,B>,M,N>; enum {size=N}; template<int... J> static constexpr type impl(seq<J...>, F f, const mat<A,M,N> & a, B                  b) { return {vec_apply<F, void, vec<A,M>, B       >::impl(make_seq<M>{}, f, a[J], b   )...}; } };
        template<class F, int M, int N, class A, class B> struct sxa_apply<F, scalars_t<A>, A,          mat<B,M,N>> { using type=mat<ret_t<F,A,B>,M,N>; enum {size=N}; template<int... J> static constexpr type impl(seq<J...>, F f, A                  a, const mat<B,M,N> & b) { return {vec_apply<F, void, A,        vec<B,M>>::impl(make_seq<M>{}, f, a,    b[J])...}; } };

        template<class F, class A          > struct axa_apply<F, scalars_t< >, quat<A>         > { using type=quat<ret_t<F,A  >>; enum {size=0}; static constexpr type impl(seq<>, F f, const quat<A> & a                   ) { return {f(a.x     ), f(a.y     ), f(a.z     ), f(a.w     )}; } };
        template<class F, class A, class B > struct axa_apply<F, scalars_t< >, quat<A>, quat<B>> { using type=quat<ret_t<F,A,B>>; enum {size=0}; static constexpr type impl(seq<>, F f, const quat<A> & a, const quat<B> & b) { return {f(a.x, b.x), f(a.y, b.y), f(a.z, b.z), f(a.w, b.w)}; } };
        template<class F, class A, class B > struct axs_apply<F, scalars_t<B>, quat<A>, B      > { using type=quat<ret_t<F,A,B>>; enum {size=0}; static constexpr type impl(seq<>, F f, const quat<A> & a, B               b) { return {f(a.x, b  ), f(a.y, b  ), f(a.z, b  ), f(a.w, b  )}; } };
        template<class F, class A, class B > struct sxa_apply<F, scalars_t<A>, A,       quat<B>> { using type=quat<ret_t<F,A,B>>; enum {size=0}; static constexpr type impl(seq<>, F f, A               a, const quat<B> & b) { return {f(a,   b.x), f(a,   b.y), f(a,   b.z), f(a,   b.w)}; } };

        template<class F, class... A                      > struct vec_apply<F, scalars_t<A...>, A...> { using type = ret_t<F,A...>; enum {size=0}; static constexpr type impl(seq<>, F f, A... a) { return f(a...); } };
        template<class F, class... T> struct any_apply : vec_apply<F,void,T...>, axa_apply<F,void,T...> , axs_apply<F,void,T...>, sxa_apply<F,void,T...>{};

        // Function objects for selecting between alternatives
        struct min    { template<class A, class B> constexpr auto operator() (A a, B b) const -> typename std::remove_reference<decltype(a<b ? a : b)>::type { return a<b ? a : b; } };
        struct max    { template<class A, class B> constexpr auto operator() (A a, B b) const -> typename std::remove_reference<decltype(a<b ? b : a)>::type { return a<b ? b : a; } };
        struct clamp  { template<class A, class B, class C> constexpr auto operator() (A a, B b, C c) const -> typename std::remove_reference<decltype(a<b ? b : a<c ? a : c)>::type { return a<b ? b : a<c ? a : c; } };
        struct select { template<class A, class B, class C> constexpr auto operator() (A a, B b, C c) const -> typename std::remove_reference<decltype(a ? b : c)>::type             { return a ? b : c; } };
        struct lerp   { template<class A, class B, class C> constexpr auto operator() (A a, B b, C c) const -> decltype(a*(1-c) + b*c)                                               { return a*(1-c) + b*c; } };

        // Function objects for applying operators
        struct op_pos { template<class A> constexpr auto operator() (A a) const -> decltype(+a) { return +a; } };
        struct op_neg { template<class A> constexpr auto operator() (A a) const -> decltype(-a) { return -a; } };
        struct op_not { template<class A> constexpr auto operator() (A a) const -> decltype(!a) { return !a; } };
        struct op_cmp { template<class A> constexpr auto operator() (A a) const -> decltype(~(a)) { return ~a; } };
        struct op_mul { template<class A, class B> constexpr auto operator() (A a, B b) const -> decltype(a * b)  { return a * b; } };
        struct op_div { template<class A, class B> constexpr auto operator() (A a, B b) const -> decltype(a / b)  { return a / b; } };
        struct op_mod { template<class A, class B> constexpr auto operator() (A a, B b) const -> decltype(a % b)  { return a % b; } };
        struct op_add { template<class A, class B> constexpr auto operator() (A a, B b) const -> decltype(a + b)  { return a + b; } };
        struct op_sub { template<class A, class B> constexpr auto operator() (A a, B b) const -> decltype(a - b)  { return a - b; } };
        struct op_lsh { template<class A, class B> constexpr auto operator() (A a, B b) const -> decltype(a << b) { return a << b; } };
        struct op_rsh { template<class A, class B> constexpr auto operator() (A a, B b) const -> decltype(a >> b) { return a >> b; } };
        struct op_lt  { template<class A, class B> constexpr auto operator() (A a, B b) const -> decltype(a < b)  { return a < b; } };
        struct op_gt  { template<class A, class B> constexpr auto operator() (A a, B b) const -> decltype(a > b)  { return a > b; } };
        struct op_le  { template<class A, class B> constexpr auto operator() (A a, B b) const -> decltype(a <= b) { return a <= b; } };
        struct op_ge  { template<class A, class B> constexpr auto operator() (A a, B b) const -> decltype(a >= b) { return a >= b; } };
        struct op_eq  { template<class A, class B> constexpr auto operator() (A a, B b) const -> decltype(a == b) { return a == b; } };
        struct op_ne  { template<class A, class B> constexpr auto operator() (A a, B b) const -> decltype(a != b) { return a != b; } };
        struct op_int { template<class A, class B> constexpr auto operator() (A a, B b) const -> decltype(a & b)  { return a & b; } };        
        struct op_xor { template<class A, class B> constexpr auto operator() (A a, B b) const -> decltype(a ^ b)  { return a ^ b; } };
        struct op_un  { template<class A, class B> constexpr auto operator() (A a, B b) const -> decltype(a | b)  { return a | b; } };
        struct op_and { template<class A, class B> constexpr auto operator() (A a, B b) const -> decltype(a && b) { return a && b; } };
        struct op_or  { template<class A, class B> constexpr auto operator() (A a, B b) const -> decltype(a || b) { return a || b; } };

        // Function objects for applying standard library math functions
        struct std_abs      { template<class A> auto operator() (A a) const -> decltype(std::abs  (a)) { return std::abs  (a); } };
        struct std_floor    { template<class A> auto operator() (A a) const -> decltype(std::floor(a)) { return std::floor(a); } };
        struct std_ceil     { template<class A> auto operator() (A a) const -> decltype(std::ceil (a)) { return std::ceil (a); } };
        struct std_exp      { template<class A> auto operator() (A a) const -> decltype(std::exp  (a)) { return std::exp  (a); } };
        struct std_log      { template<class A> auto operator() (A a) const -> decltype(std::log  (a)) { return std::log  (a); } };
        struct std_log10    { template<class A> auto operator() (A a) const -> decltype(std::log10(a)) { return std::log10(a); } };
        struct std_sqrt     { template<class A> auto operator() (A a) const -> decltype(std::sqrt (a)) { return std::sqrt (a); } };
        struct std_sin      { template<class A> auto operator() (A a) const -> decltype(std::sin  (a)) { return std::sin  (a); } };
        struct std_cos      { template<class A> auto operator() (A a) const -> decltype(std::cos  (a)) { return std::cos  (a); } };
        struct std_tan      { template<class A> auto operator() (A a) const -> decltype(std::tan  (a)) { return std::tan  (a); } };
        struct std_asin     { template<class A> auto operator() (A a) const -> decltype(std::asin (a)) { return std::asin (a); } };
        struct std_acos     { template<class A> auto operator() (A a) const -> decltype(std::acos (a)) { return std::acos (a); } };
        struct std_atan     { template<class A> auto operator() (A a) const -> decltype(std::atan (a)) { return std::atan (a); } };
        struct std_sinh     { template<class A> auto operator() (A a) const -> decltype(std::sinh (a)) { return std::sinh (a); } };
        struct std_cosh     { template<class A> auto operator() (A a) const -> decltype(std::cosh (a)) { return std::cosh (a); } };
        struct std_tanh     { template<class A> auto operator() (A a) const -> decltype(std::tanh (a)) { return std::tanh (a); } };
        struct std_round    { template<class A> auto operator() (A a) const -> decltype(std::round(a)) { return std::round(a); } };
        struct std_fmod     { template<class A, class B> auto operator() (A a, B b) const -> decltype(std::fmod    (a, b)) { return std::fmod    (a, b); } };
        struct std_pow      { template<class A, class B> auto operator() (A a, B b) const -> decltype(std::pow     (a, b)) { return std::pow     (a, b); } };
        struct std_atan2    { template<class A, class B> auto operator() (A a, B b) const -> decltype(std::atan2   (a, b)) { return std::atan2   (a, b); } };
        struct std_copysign { template<class A, class B> auto operator() (A a, B b) const -> decltype(std::copysign(a, b)) { return std::copysign(a, b); } };
    }

    //////////////////////////////////////////////
    // Implementation of named scalar accessors //
    //////////////////////////////////////////////

    template<class T, class A, int I> struct _scalar 
    { 
        A elems;
        _scalar()                          = delete;
        T operator () () const             { return elems[I]; }
        operator const T & () const        { return elems[I]; }
        operator T & ()                    { return elems[I]; }
        const T * operator & () const      { return &elems[I]; }
        T * operator & ()                  { return &elems[I]; }
        T & operator = (const T & value)   { return elems[I] = value; }
        T & operator = (const _scalar & r) { return elems[I] = r; } // Default copy operator has incorrect semantics, force interpretation as scalar
    };
    template<class T, class A, int I0, int I1> struct _lswizzle<T,A,I0,I1>
    {
        A                           e;
        constexpr                   _lswizzle()                              : e{} {} // This constructor makes _lswizzle a literal type, allowing it to live inside unions
                                    _lswizzle(T e0, T e1)                    { e[I0]=e0; e[I1]=e1; }
                                    _lswizzle(const vec<T,2> & r)            : _lswizzle(r[0], r[1]) {}
        template<class B, int... J> _lswizzle(const _lswizzle<T,B,J...> & r) : _lswizzle(vec<T,2>(r)) {}
        vec<T,2>                    operator () () const                     { return {e[I0], e[I1]}; }
                                    operator vec<T,2> () const               { return {e[I0], e[I1]}; }           
        _lswizzle &                 operator = (const _lswizzle & r)         { e[I0]=r.e[I0]; e[I1]=r.e[I1]; return *this; } // Default copy operator has incorrect semantics
    };
    template<class T, class A, int I0, int I1, int I2> struct _lswizzle<T,A,I0,I1,I2>
    {
        A                           e;
        constexpr                   _lswizzle()                              : e{} {} // This constructor makes _lswizzle a literal type, allowing it to live inside unions
                                    _lswizzle(T e0, T e1, T e2)              { e[I0]=e0; e[I1]=e1; e[I2]=e2; }
                                    _lswizzle(const vec<T,3> & r)            : _lswizzle(r[0], r[1], r[2]) {}
        template<class B, int... J> _lswizzle(const _lswizzle<T,B,J...> & r) : _lswizzle(vec<T,3>(r)) {}
        vec<T,3>                    operator () () const                     { return {e[I0], e[I1], e[I2]}; }
                                    operator vec<T,3> () const               { return {e[I0], e[I1], e[I2]}; }           
        _lswizzle &                 operator = (const _lswizzle & r)         { e[I0]=r.e[I0]; e[I1]=r.e[I1]; e[I2]=r.e[I2]; return *this; } // Default copy operator has incorrect semantics
    };
    template<class T, class A, int I0, int I1, int I2, int I3> struct _lswizzle<T,A,I0,I1,I2,I3>
    {
        A                           e;
        constexpr                   _lswizzle()                              : e{} {} // This constructor makes _lswizzle a literal type, allowing it to live inside unions
                                    _lswizzle(T e0, T e1, T e2, T e3)        { e[I0]=e0; e[I1]=e1; e[I2]=e2; e[I3]=e3; }
                                    _lswizzle(const vec<T,4> & r)            : _lswizzle(r[0], r[1], r[2], r[3]) {}
        template<class B, int... J> _lswizzle(const _lswizzle<T,B,J...> & r) : _lswizzle(vec<T,4>(r)) {}
        vec<T,4>                    operator () () const                     { return {e[I0], e[I1], e[I2], e[I3]}; }
                                    operator vec<T,4> () const               { return {e[I0], e[I1], e[I2], e[I3]}; }           
        _lswizzle &                 operator = (const _lswizzle & r)         { e[I0]=r.e[I0]; e[I1]=r.e[I1]; e[I2]=r.e[I2]; e[I3]=r.e[I3]; return *this; } // Default copy operator has incorrect semantics
    };
    template<class T, class A, int... I> struct _rswizzle
    {
        A                           e;
                                    operator vec<T,sizeof...(I)> () const    { return {e[I]...}; }
        _rswizzle &                 operator = (const _rswizzle & r)         = delete;
    };

    /////////////////////////////////////////////////////////////////
    // vec<T,M> specializations for 1, 2, 3, and 4 element vectors //
    /////////////////////////////////////////////////////////////////

    template<class T> struct vec<T,1>
    {
        union
        {
            detail::element_storage<T[1]> elems;
            _scalar<T,T[1],0> x,r,s;
            _rswizzle<T,T[1],0,0> xx,rr,ss;
            _rswizzle<T,T[1],0,0,0> xxx,rrr,sss;
            _rswizzle<T,T[1],0,0,0,0> xxxx,rrrr,ssss;
        };
        constexpr                            vec()                                                       : elems{} {}
        constexpr                            vec(const vec & v)                                          : elems{v[0]} {}
        constexpr                            vec(const T & e0)                                           : elems{e0} {}
        template<class U> constexpr explicit vec(const vec<U,1> & v)                                     : elems{static_cast<T>(v[0])} {}
        LINALG_CONSTEXPR14 vec &             operator = (const vec & r)                                  { elems = r.elems; return *this; }
        constexpr const T &                  operator[] (int i) const                                    { return elems.elems[i]; }
        LINALG_CONSTEXPR14 T &               operator[] (int i)                                          { return elems.elems[i]; }
        constexpr const T *                  data() const                                                { return elems.elems; }
        LINALG_CONSTEXPR14 T *               data()                                                      { return elems.elems; }

        template<class U, class=detail::conv_t<vec,U>> constexpr vec(const U & u)                        : vec(detail::convert<vec>(u)) {}
        template<class U, class=detail::conv_t<U,vec>> constexpr operator U () const                     { return detail::convert<U>(*this); }
    };
    template<class T> struct vec<T,2>
    {
        union
        {
            detail::element_storage<T[2]> elems;
            _scalar<T,T[2],0> x,r,s; _scalar<T,T[2],1> y,g,t;
            _rswizzle<T,T[2],0,0> xx,rr,ss; _lswizzle<T,T[2],0,1> xy,rg,st;
            _lswizzle<T,T[2],1,0> yx,gr,ts; _rswizzle<T,T[2],1,1> yy,gg,tt;
            _rswizzle<T,T[2],0,0,0> xxx,rrr,sss; _rswizzle<T,T[2],0,0,1> xxy,rrg,sst;
            _rswizzle<T,T[2],0,1,0> xyx,rgr,sts; _rswizzle<T,T[2],0,1,1> xyy,rgg,stt;
            _rswizzle<T,T[2],1,0,0> yxx,grr,tss; _rswizzle<T,T[2],1,0,1> yxy,grg,tst;
            _rswizzle<T,T[2],1,1,0> yyx,ggr,tts; _rswizzle<T,T[2],1,1,1> yyy,ggg,ttt;
            _rswizzle<T,T[2],0,0,0,0> xxxx,rrrr,ssss; _rswizzle<T,T[2],0,0,0,1> xxxy,rrrg,ssst;
            _rswizzle<T,T[2],0,0,1,0> xxyx,rrgr,ssts; _rswizzle<T,T[2],0,0,1,1> xxyy,rrgg,sstt;
            _rswizzle<T,T[2],0,1,0,0> xyxx,rgrr,stss; _rswizzle<T,T[2],0,1,0,1> xyxy,rgrg,stst;
            _rswizzle<T,T[2],0,1,1,0> xyyx,rggr,stts; _rswizzle<T,T[2],0,1,1,1> xyyy,rggg,sttt;
            _rswizzle<T,T[2],1,0,0,0> yxxx,grrr,tsss; _rswizzle<T,T[2],1,0,0,1> yxxy,grrg,tsst;
            _rswizzle<T,T[2],1,0,1,0> yxyx,grgr,tsts; _rswizzle<T,T[2],1,0,1,1> yxyy,grgg,tstt;
            _rswizzle<T,T[2],1,1,0,0> yyxx,ggrr,ttss; _rswizzle<T,T[2],1,1,0,1> yyxy,ggrg,ttst;
            _rswizzle<T,T[2],1,1,1,0> yyyx,gggr,ttts; _rswizzle<T,T[2],1,1,1,1> yyyy,gggg,tttt;
        };
        constexpr                            vec()                                                       : elems{} {}
        constexpr                            vec(const vec & v)                                          : elems{v[0], v[1]} {}
        constexpr                            vec(const T & e0, const T & e1)                             : elems{e0, e1} {}
        constexpr explicit                   vec(const T & s)                                            : elems{s, s} {}
        template<class U> constexpr explicit vec(const vec<U,2> & v)                                     : elems{static_cast<T>(v[0]), static_cast<T>(v[1])} {}
        LINALG_CONSTEXPR14 vec &             operator = (const vec & r)                                  { elems = r.elems; return *this; }
        constexpr const T &                  operator[] (int i) const                                    { return elems.elems[i]; }
        LINALG_CONSTEXPR14 T &               operator[] (int i)                                          { return elems.elems[i]; }
        constexpr const T *                  data() const                                                { return elems.elems; }
        LINALG_CONSTEXPR14 T *               data()                                                      { return elems.elems; }

        template<class U, class=detail::conv_t<vec,U>> constexpr vec(const U & u)                        : vec(detail::convert<vec>(u)) {}
        template<class U, class=detail::conv_t<U,vec>> constexpr operator U () const                     { return detail::convert<U>(*this); }
    };
    template<class T> struct vec<T,3>
    {
        union
        {
            detail::element_storage<T[3]> elems;
            _scalar<T,T[3],0> x,r,s; _scalar<T,T[3],1> y,g,t; _scalar<T,T[3],2> z,b,p;
            _rswizzle<T,T[3],0,0> xx,rr,ss; _lswizzle<T,T[3],0,1> xy,rg,st; _lswizzle<T,T[3],0,2> xz,rb,sp;
            _lswizzle<T,T[3],1,0> yx,gr,ts; _rswizzle<T,T[3],1,1> yy,gg,tt; _lswizzle<T,T[3],1,2> yz,gb,tp;
            _lswizzle<T,T[3],2,0> zx,br,ps; _lswizzle<T,T[3],2,1> zy,bg,pt; _rswizzle<T,T[3],2,2> zz,bb,pp;
            _rswizzle<T,T[3],0,0,0> xxx,rrr,sss; _rswizzle<T,T[3],0,0,1> xxy,rrg,sst; _rswizzle<T,T[3],0,0,2> xxz,rrb,ssp;
            _rswizzle<T,T[3],0,1,0> xyx,rgr,sts; _rswizzle<T,T[3],0,1,1> xyy,rgg,stt; _lswizzle<T,T[3],0,1,2> xyz,rgb,stp;
            _rswizzle<T,T[3],0,2,0> xzx,rbr,sps; _lswizzle<T,T[3],0,2,1> xzy,rbg,spt; _rswizzle<T,T[3],0,2,2> xzz,rbb,spp;
            _rswizzle<T,T[3],1,0,0> yxx,grr,tss; _rswizzle<T,T[3],1,0,1> yxy,grg,tst; _lswizzle<T,T[3],1,0,2> yxz,grb,tsp;
            _rswizzle<T,T[3],1,1,0> yyx,ggr,tts; _rswizzle<T,T[3],1,1,1> yyy,ggg,ttt; _rswizzle<T,T[3],1,1,2> yyz,ggb,ttp;
            _lswizzle<T,T[3],1,2,0> yzx,gbr,tps; _rswizzle<T,T[3],1,2,1> yzy,gbg,tpt; _rswizzle<T,T[3],1,2,2> yzz,gbb,tpp;
            _rswizzle<T,T[3],2,0,0> zxx,brr,pss; _lswizzle<T,T[3],2,0,1> zxy,brg,pst; _rswizzle<T,T[3],2,0,2> zxz,brb,psp;
            _lswizzle<T,T[3],2,1,0> zyx,bgr,pts; _rswizzle<T,T[3],2,1,1> zyy,bgg,ptt; _rswizzle<T,T[3],2,1,2> zyz,bgb,ptp;
            _rswizzle<T,T[3],2,2,0> zzx,bbr,pps; _rswizzle<T,T[3],2,2,1> zzy,bbg,ppt; _rswizzle<T,T[3],2,2,2> zzz,bbb,ppp;
            _rswizzle<T,T[3],0,0,0,0> xxxx,rrrr,ssss; _rswizzle<T,T[3],0,0,0,1> xxxy,rrrg,ssst; _rswizzle<T,T[3],0,0,0,2> xxxz,rrrb,sssp;
            _rswizzle<T,T[3],0,0,1,0> xxyx,rrgr,ssts; _rswizzle<T,T[3],0,0,1,1> xxyy,rrgg,sstt; _rswizzle<T,T[3],0,0,1,2> xxyz,rrgb,sstp;
            _rswizzle<T,T[3],0,0,2,0> xxzx,rrbr,ssps; _rswizzle<T,T[3],0,0,2,1> xxzy,rrbg,sspt; _rswizzle<T,T[3],0,0,2,2> xxzz,rrbb,sspp;
            _rswizzle<T,T[3],0,1,0,0> xyxx,rgrr,stss; _rswizzle<T,T[3],0,1,0,1> xyxy,rgrg,stst; _rswizzle<T,T[3],0,1,0,2> xyxz,rgrb,stsp;
            _rswizzle<T,T[3],0,1,1,0> xyyx,rggr,stts; _rswizzle<T,T[3],0,1,1,1> xyyy,rggg,sttt; _rswizzle<T,T[3],0,1,1,2> xyyz,rggb,sttp;
            _rswizzle<T,T[3],0,1,2,0> xyzx,rgbr,stps; _rswizzle<T,T[3],0,1,2,1> xyzy,rgbg,stpt; _rswizzle<T,T[3],0,1,2,2> xyzz,rgbb,stpp;
            _rswizzle<T,T[3],0,2,0,0> xzxx,rbrr,spss; _rswizzle<T,T[3],0,2,0,1> xzxy,rbrg,spst; _rswizzle<T,T[3],0,2,0,2> xzxz,rbrb,spsp;
            _rswizzle<T,T[3],0,2,1,0> xzyx,rbgr,spts; _rswizzle<T,T[3],0,2,1,1> xzyy,rbgg,sptt; _rswizzle<T,T[3],0,2,1,2> xzyz,rbgb,sptp;
            _rswizzle<T,T[3],0,2,2,0> xzzx,rbbr,spps; _rswizzle<T,T[3],0,2,2,1> xzzy,rbbg,sppt; _rswizzle<T,T[3],0,2,2,2> xzzz,rbbb,sppp;
            _rswizzle<T,T[3],1,0,0,0> yxxx,grrr,tsss; _rswizzle<T,T[3],1,0,0,1> yxxy,grrg,tsst; _rswizzle<T,T[3],1,0,0,2> yxxz,grrb,tssp;
            _rswizzle<T,T[3],1,0,1,0> yxyx,grgr,tsts; _rswizzle<T,T[3],1,0,1,1> yxyy,grgg,tstt; _rswizzle<T,T[3],1,0,1,2> yxyz,grgb,tstp;
            _rswizzle<T,T[3],1,0,2,0> yxzx,grbr,tsps; _rswizzle<T,T[3],1,0,2,1> yxzy,grbg,tspt; _rswizzle<T,T[3],1,0,2,2> yxzz,grbb,tspp;
            _rswizzle<T,T[3],1,1,0,0> yyxx,ggrr,ttss; _rswizzle<T,T[3],1,1,0,1> yyxy,ggrg,ttst; _rswizzle<T,T[3],1,1,0,2> yyxz,ggrb,ttsp;
            _rswizzle<T,T[3],1,1,1,0> yyyx,gggr,ttts; _rswizzle<T,T[3],1,1,1,1> yyyy,gggg,tttt; _rswizzle<T,T[3],1,1,1,2> yyyz,gggb,tttp;
            _rswizzle<T,T[3],1,1,2,0> yyzx,ggbr,ttps; _rswizzle<T,T[3],1,1,2,1> yyzy,ggbg,ttpt; _rswizzle<T,T[3],1,1,2,2> yyzz,ggbb,ttpp;
            _rswizzle<T,T[3],1,2,0,0> yzxx,gbrr,tpss; _rswizzle<T,T[3],1,2,0,1> yzxy,gbrg,tpst; _rswizzle<T,T[3],1,2,0,2> yzxz,gbrb,tpsp;
            _rswizzle<T,T[3],1,2,1,0> yzyx,gbgr,tpts; _rswizzle<T,T[3],1,2,1,1> yzyy,gbgg,tptt; _rswizzle<T,T[3],1,2,1,2> yzyz,gbgb,tptp;
            _rswizzle<T,T[3],1,2,2,0> yzzx,gbbr,tpps; _rswizzle<T,T[3],1,2,2,1> yzzy,gbbg,tppt; _rswizzle<T,T[3],1,2,2,2> yzzz,gbbb,tppp;
            _rswizzle<T,T[3],2,0,0,0> zxxx,brrr,psss; _rswizzle<T,T[3],2,0,0,1> zxxy,brrg,psst; _rswizzle<T,T[3],2,0,0,2> zxxz,brrb,pssp;
            _rswizzle<T,T[3],2,0,1,0> zxyx,brgr,psts; _rswizzle<T,T[3],2,0,1,1> zxyy,brgg,pstt; _rswizzle<T,T[3],2,0,1,2> zxyz,brgb,pstp;
            _rswizzle<T,T[3],2,0,2,0> zxzx,brbr,psps; _rswizzle<T,T[3],2,0,2,1> zxzy,brbg,pspt; _rswizzle<T,T[3],2,0,2,2> zxzz,brbb,pspp;
            _rswizzle<T,T[3],2,1,0,0> zyxx,bgrr,ptss; _rswizzle<T,T[3],2,1,0,1> zyxy,bgrg,ptst; _rswizzle<T,T[3],2,1,0,2> zyxz,bgrb,ptsp;
            _rswizzle<T,T[3],2,1,1,0> zyyx,bggr,ptts; _rswizzle<T,T[3],2,1,1,1> zyyy,bggg,pttt; _rswizzle<T,T[3],2,1,1,2> zyyz,bggb,pttp;
            _rswizzle<T,T[3],2,1,2,0> zyzx,bgbr,ptps; _rswizzle<T,T[3],2,1,2,1> zyzy,bgbg,ptpt; _rswizzle<T,T[3],2,1,2,2> zyzz,bgbb,ptpp;
            _rswizzle<T,T[3],2,2,0,0> zzxx,bbrr,ppss; _rswizzle<T,T[3],2,2,0,1> zzxy,bbrg,ppst; _rswizzle<T,T[3],2,2,0,2> zzxz,bbrb,ppsp;
            _rswizzle<T,T[3],2,2,1,0> zzyx,bbgr,ppts; _rswizzle<T,T[3],2,2,1,1> zzyy,bbgg,pptt; _rswizzle<T,T[3],2,2,1,2> zzyz,bbgb,pptp;
            _rswizzle<T,T[3],2,2,2,0> zzzx,bbbr,ppps; _rswizzle<T,T[3],2,2,2,1> zzzy,bbbg,pppt; _rswizzle<T,T[3],2,2,2,2> zzzz,bbbb,pppp;
        };
        constexpr                            vec()                                                       : elems{} {}
        constexpr                            vec(const vec & v)                                          : elems{v[0], v[1], v[2]} {}
        constexpr                            vec(const T & e0, const T & e1, const T & e2)               : elems{e0, e1, e2} {}
        constexpr                            vec(const vec<T,2> & e01, const T & e2)                     : elems{e01[0], e01[1], e2} {}
        constexpr explicit                   vec(const T & s)                                            : elems{s, s, s} {}
        template<class U> constexpr explicit vec(const vec<U,3> & v)                                     : elems{static_cast<T>(v[0]), static_cast<T>(v[1]), static_cast<T>(v[2])} {}
        LINALG_CONSTEXPR14 vec &             operator = (const vec & r)                                  { elems = r.elems; return *this; }
        constexpr const T &                  operator[] (int i) const                                    { return elems.elems[i]; }
        LINALG_CONSTEXPR14 T &               operator[] (int i)                                          { return elems.elems[i]; }
        constexpr const T *                  data() const                                                { return elems.elems; }
        LINALG_CONSTEXPR14 T *               data()                                                      { return elems.elems; }

        template<class U, class=detail::conv_t<vec,U>> constexpr vec(const U & u)                        : vec(detail::convert<vec>(u)) {}
        template<class U, class=detail::conv_t<U,vec>> constexpr operator U () const                     { return detail::convert<U>(*this); }
    };
    template<class T> struct vec<T,4>
    {
        union
        {
            detail::element_storage<T[4]> elems;
            _scalar<T,T[4],0> x,r,s; _scalar<T,T[4],1> y,g,t; _scalar<T,T[4],2> z,b,p; _scalar<T,T[4],3> w,a,q;
            _rswizzle<T,T[4],0,0> xx,rr,ss; _lswizzle<T,T[4],0,1> xy,rg,st; _lswizzle<T,T[4],0,2> xz,rb,sp; _lswizzle<T,T[4],0,3> xw,ra,sq;
            _lswizzle<T,T[4],1,0> yx,gr,ts; _rswizzle<T,T[4],1,1> yy,gg,tt; _lswizzle<T,T[4],1,2> yz,gb,tp; _lswizzle<T,T[4],1,3> yw,ga,tq;
            _lswizzle<T,T[4],2,0> zx,br,ps; _lswizzle<T,T[4],2,1> zy,bg,pt; _rswizzle<T,T[4],2,2> zz,bb,pp; _lswizzle<T,T[4],2,3> zw,ba,pq;
            _lswizzle<T,T[4],3,0> wx,ar,qs; _lswizzle<T,T[4],3,1> wy,ag,qt; _lswizzle<T,T[4],3,2> wz,ab,qp; _rswizzle<T,T[4],3,3> ww,aa,qq;
            _rswizzle<T,T[4],0,0,0> xxx,rrr,sss; _rswizzle<T,T[4],0,0,1> xxy,rrg,sst; _rswizzle<T,T[4],0,0,2> xxz,rrb,ssp; _rswizzle<T,T[4],0,0,3> xxw,rra,ssq;
            _rswizzle<T,T[4],0,1,0> xyx,rgr,sts; _rswizzle<T,T[4],0,1,1> xyy,rgg,stt; _lswizzle<T,T[4],0,1,2> xyz,rgb,stp; _lswizzle<T,T[4],0,1,3> xyw,rga,stq;
            _rswizzle<T,T[4],0,2,0> xzx,rbr,sps; _lswizzle<T,T[4],0,2,1> xzy,rbg,spt; _rswizzle<T,T[4],0,2,2> xzz,rbb,spp; _lswizzle<T,T[4],0,2,3> xzw,rba,spq;
            _rswizzle<T,T[4],0,3,0> xwx,rar,sqs; _lswizzle<T,T[4],0,3,1> xwy,rag,sqt; _lswizzle<T,T[4],0,3,2> xwz,rab,sqp; _rswizzle<T,T[4],0,3,3> xww,raa,sqq;
            _rswizzle<T,T[4],1,0,0> yxx,grr,tss; _rswizzle<T,T[4],1,0,1> yxy,grg,tst; _lswizzle<T,T[4],1,0,2> yxz,grb,tsp; _lswizzle<T,T[4],1,0,3> yxw,gra,tsq;
            _rswizzle<T,T[4],1,1,0> yyx,ggr,tts; _rswizzle<T,T[4],1,1,1> yyy,ggg,ttt; _rswizzle<T,T[4],1,1,2> yyz,ggb,ttp; _rswizzle<T,T[4],1,1,3> yyw,gga,ttq;
            _lswizzle<T,T[4],1,2,0> yzx,gbr,tps; _rswizzle<T,T[4],1,2,1> yzy,gbg,tpt; _rswizzle<T,T[4],1,2,2> yzz,gbb,tpp; _lswizzle<T,T[4],1,2,3> yzw,gba,tpq;
            _lswizzle<T,T[4],1,3,0> ywx,gar,tqs; _rswizzle<T,T[4],1,3,1> ywy,gag,tqt; _lswizzle<T,T[4],1,3,2> ywz,gab,tqp; _rswizzle<T,T[4],1,3,3> yww,gaa,tqq;
            _rswizzle<T,T[4],2,0,0> zxx,brr,pss; _lswizzle<T,T[4],2,0,1> zxy,brg,pst; _rswizzle<T,T[4],2,0,2> zxz,brb,psp; _lswizzle<T,T[4],2,0,3> zxw,bra,psq;
            _lswizzle<T,T[4],2,1,0> zyx,bgr,pts; _rswizzle<T,T[4],2,1,1> zyy,bgg,ptt; _rswizzle<T,T[4],2,1,2> zyz,bgb,ptp; _lswizzle<T,T[4],2,1,3> zyw,bga,ptq;
            _rswizzle<T,T[4],2,2,0> zzx,bbr,pps; _rswizzle<T,T[4],2,2,1> zzy,bbg,ppt; _rswizzle<T,T[4],2,2,2> zzz,bbb,ppp; _rswizzle<T,T[4],2,2,3> zzw,bba,ppq;
            _lswizzle<T,T[4],2,3,0> zwx,bar,pqs; _lswizzle<T,T[4],2,3,1> zwy,bag,pqt; _rswizzle<T,T[4],2,3,2> zwz,bab,pqp; _rswizzle<T,T[4],2,3,3> zww,baa,pqq;
            _rswizzle<T,T[4],3,0,0> wxx,arr,qss; _lswizzle<T,T[4],3,0,1> wxy,arg,qst; _lswizzle<T,T[4],3,0,2> wxz,arb,qsp; _rswizzle<T,T[4],3,0,3> wxw,ara,qsq;
            _lswizzle<T,T[4],3,1,0> wyx,agr,qts; _rswizzle<T,T[4],3,1,1> wyy,agg,qtt; _lswizzle<T,T[4],3,1,2> wyz,agb,qtp; _rswizzle<T,T[4],3,1,3> wyw,aga,qtq;
            _lswizzle<T,T[4],3,2,0> wzx,abr,qps; _lswizzle<T,T[4],3,2,1> wzy,abg,qpt; _rswizzle<T,T[4],3,2,2> wzz,abb,qpp; _rswizzle<T,T[4],3,2,3> wzw,aba,qpq;
            _rswizzle<T,T[4],3,3,0> wwx,aar,qqs; _rswizzle<T,T[4],3,3,1> wwy,aag,qqt; _rswizzle<T,T[4],3,3,2> wwz,aab,qqp; _rswizzle<T,T[4],3,3,3> www,aaa,qqq;
            _rswizzle<T,T[4],0,0,0,0> xxxx,rrrr,ssss; _rswizzle<T,T[4],0,0,0,1> xxxy,rrrg,ssst; _rswizzle<T,T[4],0,0,0,2> xxxz,rrrb,sssp; _rswizzle<T,T[4],0,0,0,3> xxxw,rrra,sssq;
            _rswizzle<T,T[4],0,0,1,0> xxyx,rrgr,ssts; _rswizzle<T,T[4],0,0,1,1> xxyy,rrgg,sstt; _rswizzle<T,T[4],0,0,1,2> xxyz,rrgb,sstp; _rswizzle<T,T[4],0,0,1,3> xxyw,rrga,sstq;
            _rswizzle<T,T[4],0,0,2,0> xxzx,rrbr,ssps; _rswizzle<T,T[4],0,0,2,1> xxzy,rrbg,sspt; _rswizzle<T,T[4],0,0,2,2> xxzz,rrbb,sspp; _rswizzle<T,T[4],0,0,2,3> xxzw,rrba,sspq;
            _rswizzle<T,T[4],0,0,3,0> xxwx,rrar,ssqs; _rswizzle<T,T[4],0,0,3,1> xxwy,rrag,ssqt; _rswizzle<T,T[4],0,0,3,2> xxwz,rrab,ssqp; _rswizzle<T,T[4],0,0,3,3> xxww,rraa,ssqq;
            _rswizzle<T,T[4],0,1,0,0> xyxx,rgrr,stss; _rswizzle<T,T[4],0,1,0,1> xyxy,rgrg,stst; _rswizzle<T,T[4],0,1,0,2> xyxz,rgrb,stsp; _rswizzle<T,T[4],0,1,0,3> xyxw,rgra,stsq;
            _rswizzle<T,T[4],0,1,1,0> xyyx,rggr,stts; _rswizzle<T,T[4],0,1,1,1> xyyy,rggg,sttt; _rswizzle<T,T[4],0,1,1,2> xyyz,rggb,sttp; _rswizzle<T,T[4],0,1,1,3> xyyw,rgga,sttq;
            _rswizzle<T,T[4],0,1,2,0> xyzx,rgbr,stps; _rswizzle<T,T[4],0,1,2,1> xyzy,rgbg,stpt; _rswizzle<T,T[4],0,1,2,2> xyzz,rgbb,stpp; _lswizzle<T,T[4],0,1,2,3> xyzw,rgba,stpq;
            _rswizzle<T,T[4],0,1,3,0> xywx,rgar,stqs; _rswizzle<T,T[4],0,1,3,1> xywy,rgag,stqt; _lswizzle<T,T[4],0,1,3,2> xywz,rgab,stqp; _rswizzle<T,T[4],0,1,3,3> xyww,rgaa,stqq;
            _rswizzle<T,T[4],0,2,0,0> xzxx,rbrr,spss; _rswizzle<T,T[4],0,2,0,1> xzxy,rbrg,spst; _rswizzle<T,T[4],0,2,0,2> xzxz,rbrb,spsp; _rswizzle<T,T[4],0,2,0,3> xzxw,rbra,spsq;
            _rswizzle<T,T[4],0,2,1,0> xzyx,rbgr,spts; _rswizzle<T,T[4],0,2,1,1> xzyy,rbgg,sptt; _rswizzle<T,T[4],0,2,1,2> xzyz,rbgb,sptp; _lswizzle<T,T[4],0,2,1,3> xzyw,rbga,sptq;
            _rswizzle<T,T[4],0,2,2,0> xzzx,rbbr,spps; _rswizzle<T,T[4],0,2,2,1> xzzy,rbbg,sppt; _rswizzle<T,T[4],0,2,2,2> xzzz,rbbb,sppp; _rswizzle<T,T[4],0,2,2,3> xzzw,rbba,sppq;
            _rswizzle<T,T[4],0,2,3,0> xzwx,rbar,spqs; _lswizzle<T,T[4],0,2,3,1> xzwy,rbag,spqt; _rswizzle<T,T[4],0,2,3,2> xzwz,rbab,spqp; _rswizzle<T,T[4],0,2,3,3> xzww,rbaa,spqq;
            _rswizzle<T,T[4],0,3,0,0> xwxx,rarr,sqss; _rswizzle<T,T[4],0,3,0,1> xwxy,rarg,sqst; _rswizzle<T,T[4],0,3,0,2> xwxz,rarb,sqsp; _rswizzle<T,T[4],0,3,0,3> xwxw,rara,sqsq;
            _rswizzle<T,T[4],0,3,1,0> xwyx,ragr,sqts; _rswizzle<T,T[4],0,3,1,1> xwyy,ragg,sqtt; _lswizzle<T,T[4],0,3,1,2> xwyz,ragb,sqtp; _rswizzle<T,T[4],0,3,1,3> xwyw,raga,sqtq;
            _rswizzle<T,T[4],0,3,2,0> xwzx,rabr,sqps; _lswizzle<T,T[4],0,3,2,1> xwzy,rabg,sqpt; _rswizzle<T,T[4],0,3,2,2> xwzz,rabb,sqpp; _rswizzle<T,T[4],0,3,2,3> xwzw,raba,sqpq;
            _rswizzle<T,T[4],0,3,3,0> xwwx,raar,sqqs; _rswizzle<T,T[4],0,3,3,1> xwwy,raag,sqqt; _rswizzle<T,T[4],0,3,3,2> xwwz,raab,sqqp; _rswizzle<T,T[4],0,3,3,3> xwww,raaa,sqqq;
            _rswizzle<T,T[4],1,0,0,0> yxxx,grrr,tsss; _rswizzle<T,T[4],1,0,0,1> yxxy,grrg,tsst; _rswizzle<T,T[4],1,0,0,2> yxxz,grrb,tssp; _rswizzle<T,T[4],1,0,0,3> yxxw,grra,tssq;
            _rswizzle<T,T[4],1,0,1,0> yxyx,grgr,tsts; _rswizzle<T,T[4],1,0,1,1> yxyy,grgg,tstt; _rswizzle<T,T[4],1,0,1,2> yxyz,grgb,tstp; _rswizzle<T,T[4],1,0,1,3> yxyw,grga,tstq;
            _rswizzle<T,T[4],1,0,2,0> yxzx,grbr,tsps; _rswizzle<T,T[4],1,0,2,1> yxzy,grbg,tspt; _rswizzle<T,T[4],1,0,2,2> yxzz,grbb,tspp; _lswizzle<T,T[4],1,0,2,3> yxzw,grba,tspq;
            _rswizzle<T,T[4],1,0,3,0> yxwx,grar,tsqs; _rswizzle<T,T[4],1,0,3,1> yxwy,grag,tsqt; _lswizzle<T,T[4],1,0,3,2> yxwz,grab,tsqp; _rswizzle<T,T[4],1,0,3,3> yxww,graa,tsqq;
            _rswizzle<T,T[4],1,1,0,0> yyxx,ggrr,ttss; _rswizzle<T,T[4],1,1,0,1> yyxy,ggrg,ttst; _rswizzle<T,T[4],1,1,0,2> yyxz,ggrb,ttsp; _rswizzle<T,T[4],1,1,0,3> yyxw,ggra,ttsq;
            _rswizzle<T,T[4],1,1,1,0> yyyx,gggr,ttts; _rswizzle<T,T[4],1,1,1,1> yyyy,gggg,tttt; _rswizzle<T,T[4],1,1,1,2> yyyz,gggb,tttp; _rswizzle<T,T[4],1,1,1,3> yyyw,ggga,tttq;
            _rswizzle<T,T[4],1,1,2,0> yyzx,ggbr,ttps; _rswizzle<T,T[4],1,1,2,1> yyzy,ggbg,ttpt; _rswizzle<T,T[4],1,1,2,2> yyzz,ggbb,ttpp; _rswizzle<T,T[4],1,1,2,3> yyzw,ggba,ttpq;
            _rswizzle<T,T[4],1,1,3,0> yywx,ggar,ttqs; _rswizzle<T,T[4],1,1,3,1> yywy,ggag,ttqt; _rswizzle<T,T[4],1,1,3,2> yywz,ggab,ttqp; _rswizzle<T,T[4],1,1,3,3> yyww,ggaa,ttqq;
            _rswizzle<T,T[4],1,2,0,0> yzxx,gbrr,tpss; _rswizzle<T,T[4],1,2,0,1> yzxy,gbrg,tpst; _rswizzle<T,T[4],1,2,0,2> yzxz,gbrb,tpsp; _lswizzle<T,T[4],1,2,0,3> yzxw,gbra,tpsq;
            _rswizzle<T,T[4],1,2,1,0> yzyx,gbgr,tpts; _rswizzle<T,T[4],1,2,1,1> yzyy,gbgg,tptt; _rswizzle<T,T[4],1,2,1,2> yzyz,gbgb,tptp; _rswizzle<T,T[4],1,2,1,3> yzyw,gbga,tptq;
            _rswizzle<T,T[4],1,2,2,0> yzzx,gbbr,tpps; _rswizzle<T,T[4],1,2,2,1> yzzy,gbbg,tppt; _rswizzle<T,T[4],1,2,2,2> yzzz,gbbb,tppp; _rswizzle<T,T[4],1,2,2,3> yzzw,gbba,tppq;
            _lswizzle<T,T[4],1,2,3,0> yzwx,gbar,tpqs; _rswizzle<T,T[4],1,2,3,1> yzwy,gbag,tpqt; _rswizzle<T,T[4],1,2,3,2> yzwz,gbab,tpqp; _rswizzle<T,T[4],1,2,3,3> yzww,gbaa,tpqq;
            _rswizzle<T,T[4],1,3,0,0> ywxx,garr,tqss; _rswizzle<T,T[4],1,3,0,1> ywxy,garg,tqst; _lswizzle<T,T[4],1,3,0,2> ywxz,garb,tqsp; _rswizzle<T,T[4],1,3,0,3> ywxw,gara,tqsq;
            _rswizzle<T,T[4],1,3,1,0> ywyx,gagr,tqts; _rswizzle<T,T[4],1,3,1,1> ywyy,gagg,tqtt; _rswizzle<T,T[4],1,3,1,2> ywyz,gagb,tqtp; _rswizzle<T,T[4],1,3,1,3> ywyw,gaga,tqtq;
            _lswizzle<T,T[4],1,3,2,0> ywzx,gabr,tqps; _rswizzle<T,T[4],1,3,2,1> ywzy,gabg,tqpt; _rswizzle<T,T[4],1,3,2,2> ywzz,gabb,tqpp; _rswizzle<T,T[4],1,3,2,3> ywzw,gaba,tqpq;
            _rswizzle<T,T[4],1,3,3,0> ywwx,gaar,tqqs; _rswizzle<T,T[4],1,3,3,1> ywwy,gaag,tqqt; _rswizzle<T,T[4],1,3,3,2> ywwz,gaab,tqqp; _rswizzle<T,T[4],1,3,3,3> ywww,gaaa,tqqq;
            _rswizzle<T,T[4],2,0,0,0> zxxx,brrr,psss; _rswizzle<T,T[4],2,0,0,1> zxxy,brrg,psst; _rswizzle<T,T[4],2,0,0,2> zxxz,brrb,pssp; _rswizzle<T,T[4],2,0,0,3> zxxw,brra,pssq;
            _rswizzle<T,T[4],2,0,1,0> zxyx,brgr,psts; _rswizzle<T,T[4],2,0,1,1> zxyy,brgg,pstt; _rswizzle<T,T[4],2,0,1,2> zxyz,brgb,pstp; _lswizzle<T,T[4],2,0,1,3> zxyw,brga,pstq;
            _rswizzle<T,T[4],2,0,2,0> zxzx,brbr,psps; _rswizzle<T,T[4],2,0,2,1> zxzy,brbg,pspt; _rswizzle<T,T[4],2,0,2,2> zxzz,brbb,pspp; _rswizzle<T,T[4],2,0,2,3> zxzw,brba,pspq;
            _rswizzle<T,T[4],2,0,3,0> zxwx,brar,psqs; _lswizzle<T,T[4],2,0,3,1> zxwy,brag,psqt; _rswizzle<T,T[4],2,0,3,2> zxwz,brab,psqp; _rswizzle<T,T[4],2,0,3,3> zxww,braa,psqq;
            _rswizzle<T,T[4],2,1,0,0> zyxx,bgrr,ptss; _rswizzle<T,T[4],2,1,0,1> zyxy,bgrg,ptst; _rswizzle<T,T[4],2,1,0,2> zyxz,bgrb,ptsp; _lswizzle<T,T[4],2,1,0,3> zyxw,bgra,ptsq;
            _rswizzle<T,T[4],2,1,1,0> zyyx,bggr,ptts; _rswizzle<T,T[4],2,1,1,1> zyyy,bggg,pttt; _rswizzle<T,T[4],2,1,1,2> zyyz,bggb,pttp; _rswizzle<T,T[4],2,1,1,3> zyyw,bgga,pttq;
            _rswizzle<T,T[4],2,1,2,0> zyzx,bgbr,ptps; _rswizzle<T,T[4],2,1,2,1> zyzy,bgbg,ptpt; _rswizzle<T,T[4],2,1,2,2> zyzz,bgbb,ptpp; _rswizzle<T,T[4],2,1,2,3> zyzw,bgba,ptpq;
            _lswizzle<T,T[4],2,1,3,0> zywx,bgar,ptqs; _rswizzle<T,T[4],2,1,3,1> zywy,bgag,ptqt; _rswizzle<T,T[4],2,1,3,2> zywz,bgab,ptqp; _rswizzle<T,T[4],2,1,3,3> zyww,bgaa,ptqq;
            _rswizzle<T,T[4],2,2,0,0> zzxx,bbrr,ppss; _rswizzle<T,T[4],2,2,0,1> zzxy,bbrg,ppst; _rswizzle<T,T[4],2,2,0,2> zzxz,bbrb,ppsp; _rswizzle<T,T[4],2,2,0,3> zzxw,bbra,ppsq;
            _rswizzle<T,T[4],2,2,1,0> zzyx,bbgr,ppts; _rswizzle<T,T[4],2,2,1,1> zzyy,bbgg,pptt; _rswizzle<T,T[4],2,2,1,2> zzyz,bbgb,pptp; _rswizzle<T,T[4],2,2,1,3> zzyw,bbga,pptq;
            _rswizzle<T,T[4],2,2,2,0> zzzx,bbbr,ppps; _rswizzle<T,T[4],2,2,2,1> zzzy,bbbg,pppt; _rswizzle<T,T[4],2,2,2,2> zzzz,bbbb,pppp; _rswizzle<T,T[4],2,2,2,3> zzzw,bbba,pppq;
            _rswizzle<T,T[4],2,2,3,0> zzwx,bbar,ppqs; _rswizzle<T,T[4],2,2,3,1> zzwy,bbag,ppqt; _rswizzle<T,T[4],2,2,3,2> zzwz,bbab,ppqp; _rswizzle<T,T[4],2,2,3,3> zzww,bbaa,ppqq;
            _rswizzle<T,T[4],2,3,0,0> zwxx,barr,pqss; _lswizzle<T,T[4],2,3,0,1> zwxy,barg,pqst; _rswizzle<T,T[4],2,3,0,2> zwxz,barb,pqsp; _rswizzle<T,T[4],2,3,0,3> zwxw,bara,pqsq;
            _lswizzle<T,T[4],2,3,1,0> zwyx,bagr,pqts; _rswizzle<T,T[4],2,3,1,1> zwyy,bagg,pqtt; _rswizzle<T,T[4],2,3,1,2> zwyz,bagb,pqtp; _rswizzle<T,T[4],2,3,1,3> zwyw,baga,pqtq;
            _rswizzle<T,T[4],2,3,2,0> zwzx,babr,pqps; _rswizzle<T,T[4],2,3,2,1> zwzy,babg,pqpt; _rswizzle<T,T[4],2,3,2,2> zwzz,babb,pqpp; _rswizzle<T,T[4],2,3,2,3> zwzw,baba,pqpq;
            _rswizzle<T,T[4],2,3,3,0> zwwx,baar,pqqs; _rswizzle<T,T[4],2,3,3,1> zwwy,baag,pqqt; _rswizzle<T,T[4],2,3,3,2> zwwz,baab,pqqp; _rswizzle<T,T[4],2,3,3,3> zwww,baaa,pqqq;
            _rswizzle<T,T[4],3,0,0,0> wxxx,arrr,qsss; _rswizzle<T,T[4],3,0,0,1> wxxy,arrg,qsst; _rswizzle<T,T[4],3,0,0,2> wxxz,arrb,qssp; _rswizzle<T,T[4],3,0,0,3> wxxw,arra,qssq;
            _rswizzle<T,T[4],3,0,1,0> wxyx,argr,qsts; _rswizzle<T,T[4],3,0,1,1> wxyy,argg,qstt; _lswizzle<T,T[4],3,0,1,2> wxyz,argb,qstp; _rswizzle<T,T[4],3,0,1,3> wxyw,arga,qstq;
            _rswizzle<T,T[4],3,0,2,0> wxzx,arbr,qsps; _lswizzle<T,T[4],3,0,2,1> wxzy,arbg,qspt; _rswizzle<T,T[4],3,0,2,2> wxzz,arbb,qspp; _rswizzle<T,T[4],3,0,2,3> wxzw,arba,qspq;
            _rswizzle<T,T[4],3,0,3,0> wxwx,arar,qsqs; _rswizzle<T,T[4],3,0,3,1> wxwy,arag,qsqt; _rswizzle<T,T[4],3,0,3,2> wxwz,arab,qsqp; _rswizzle<T,T[4],3,0,3,3> wxww,araa,qsqq;
            _rswizzle<T,T[4],3,1,0,0> wyxx,agrr,qtss; _rswizzle<T,T[4],3,1,0,1> wyxy,agrg,qtst; _lswizzle<T,T[4],3,1,0,2> wyxz,agrb,qtsp; _rswizzle<T,T[4],3,1,0,3> wyxw,agra,qtsq;
            _rswizzle<T,T[4],3,1,1,0> wyyx,aggr,qtts; _rswizzle<T,T[4],3,1,1,1> wyyy,aggg,qttt; _rswizzle<T,T[4],3,1,1,2> wyyz,aggb,qttp; _rswizzle<T,T[4],3,1,1,3> wyyw,agga,qttq;
            _lswizzle<T,T[4],3,1,2,0> wyzx,agbr,qtps; _rswizzle<T,T[4],3,1,2,1> wyzy,agbg,qtpt; _rswizzle<T,T[4],3,1,2,2> wyzz,agbb,qtpp; _rswizzle<T,T[4],3,1,2,3> wyzw,agba,qtpq;
            _rswizzle<T,T[4],3,1,3,0> wywx,agar,qtqs; _rswizzle<T,T[4],3,1,3,1> wywy,agag,qtqt; _rswizzle<T,T[4],3,1,3,2> wywz,agab,qtqp; _rswizzle<T,T[4],3,1,3,3> wyww,agaa,qtqq;
            _rswizzle<T,T[4],3,2,0,0> wzxx,abrr,qpss; _lswizzle<T,T[4],3,2,0,1> wzxy,abrg,qpst; _rswizzle<T,T[4],3,2,0,2> wzxz,abrb,qpsp; _rswizzle<T,T[4],3,2,0,3> wzxw,abra,qpsq;
            _lswizzle<T,T[4],3,2,1,0> wzyx,abgr,qpts; _rswizzle<T,T[4],3,2,1,1> wzyy,abgg,qptt; _rswizzle<T,T[4],3,2,1,2> wzyz,abgb,qptp; _rswizzle<T,T[4],3,2,1,3> wzyw,abga,qptq;
            _rswizzle<T,T[4],3,2,2,0> wzzx,abbr,qpps; _rswizzle<T,T[4],3,2,2,1> wzzy,abbg,qppt; _rswizzle<T,T[4],3,2,2,2> wzzz,abbb,qppp; _rswizzle<T,T[4],3,2,2,3> wzzw,abba,qppq;
            _rswizzle<T,T[4],3,2,3,0> wzwx,abar,qpqs; _rswizzle<T,T[4],3,2,3,1> wzwy,abag,qpqt; _rswizzle<T,T[4],3,2,3,2> wzwz,abab,qpqp; _rswizzle<T,T[4],3,2,3,3> wzww,abaa,qpqq;
            _rswizzle<T,T[4],3,3,0,0> wwxx,aarr,qqss; _rswizzle<T,T[4],3,3,0,1> wwxy,aarg,qqst; _rswizzle<T,T[4],3,3,0,2> wwxz,aarb,qqsp; _rswizzle<T,T[4],3,3,0,3> wwxw,aara,qqsq;
            _rswizzle<T,T[4],3,3,1,0> wwyx,aagr,qqts; _rswizzle<T,T[4],3,3,1,1> wwyy,aagg,qqtt; _rswizzle<T,T[4],3,3,1,2> wwyz,aagb,qqtp; _rswizzle<T,T[4],3,3,1,3> wwyw,aaga,qqtq;
            _rswizzle<T,T[4],3,3,2,0> wwzx,aabr,qqps; _rswizzle<T,T[4],3,3,2,1> wwzy,aabg,qqpt; _rswizzle<T,T[4],3,3,2,2> wwzz,aabb,qqpp; _rswizzle<T,T[4],3,3,2,3> wwzw,aaba,qqpq;
            _rswizzle<T,T[4],3,3,3,0> wwwx,aaar,qqqs; _rswizzle<T,T[4],3,3,3,1> wwwy,aaag,qqqt; _rswizzle<T,T[4],3,3,3,2> wwwz,aaab,qqqp; _rswizzle<T,T[4],3,3,3,3> wwww,aaaa,qqqq;
        };
        constexpr                            vec()                                                       : elems{} {}
        constexpr                            vec(const vec & v)                                          : elems{v[0], v[1], v[2], v[3]} {}
        constexpr                            vec(const T & e0, const T & e1, const T & e2, const T & e3) : elems{e0, e1, e2, e3} {}
        constexpr                            vec(const vec<T,2> & e01, const T & e2, const T & e3)       : elems{e01[0], e01[1], e2, e3} {}
        constexpr                            vec(const vec<T,3> & e012, const T & e3)                    : elems{e012[0], e012[1], e012[2], e3} {}
        constexpr explicit                   vec(const T & s)                                            : elems{s, s, s, s} {}
        template<class U> constexpr explicit vec(const vec<U,4> & v)                                     : elems{static_cast<T>(v[0]), static_cast<T>(v[1]), static_cast<T>(v[2]), static_cast<T>(v[3])} {}
        constexpr explicit                   vec(const quat<T> & q)                                      : elems{q.x, q.y, q.z, q.w} {}
        LINALG_CONSTEXPR14 vec &             operator = (const vec & r)                                  { elems = r.elems; return *this; }
        constexpr const T &                  operator[] (int i) const                                    { return elems.elems[i]; }
        LINALG_CONSTEXPR14 T &               operator[] (int i)                                          { return elems.elems[i]; }
        constexpr const T *                  data() const                                                { return elems.elems; }
        LINALG_CONSTEXPR14 T *               data()                                                      { return elems.elems; }

        template<class U, class=detail::conv_t<vec,U>> constexpr vec(const U & u)                        : vec(detail::convert<vec>(u)) {}
        template<class U, class=detail::conv_t<U,vec>> constexpr operator U () const                     { return detail::convert<U>(*this); }
    };

    ///////////////////////////////////////////////////////////////////
    // mat<T,M,N> specializations for 1, 2, 3, and 4 column matrices //
    ///////////////////////////////////////////////////////////////////
   
    template<class T, int M> struct mat<T,M,1>
    {
        using V=vec<T,M>;
        V                                    cols[1];
        constexpr                            mat()                                                       : cols{} {}
        constexpr                            mat(const V & x_)                                           : cols{x_} {}
        constexpr explicit                   mat(const T & s)                                            : cols{V(s)} {}
        template<class U> constexpr explicit mat(const mat<U,M,1> & m)                                   : cols{V(m[0])} {}
        constexpr const V &                  operator[] (int j) const                                    { return cols[j]; }
        LINALG_CONSTEXPR14 V &               operator[] (int j)                                          { return cols[j]; }
        constexpr vec<T,1>                   row(int i) const                                            { return {cols[0][i]}; }
        constexpr const T *                  data() const                                                { return cols[0].data(); }
        LINALG_CONSTEXPR14 T *               data()                                                      { return cols[0].data(); }

        template<class U, class=detail::conv_t<mat,U>> constexpr mat(const U & u)                        : mat(detail::convert<mat>(u)) {}
        template<class U, class=detail::conv_t<U,mat>> constexpr operator U () const                     { return detail::convert<U>(*this); }
    };
    template<class T, int M> struct mat<T,M,2>
    {
        using V=vec<T,M>;
        V                                    cols[2];
        constexpr                            mat()                                                       : cols{} {}
        constexpr                            mat(const V & x_, const V & y_)                             : cols{x_, y_} {}
        constexpr explicit                   mat(const T & s)                                            : cols{V(s), V(s)} {}
        template<class U> constexpr explicit mat(const mat<U,M,2> & m)                                   : cols{V(m[0]), V(m[1])} {}
        constexpr const V &                  operator[] (int j) const                                    { return cols[j]; }
        LINALG_CONSTEXPR14 V &               operator[] (int j)                                          { return cols[j]; }
        constexpr vec<T,2>                   row(int i) const                                            { return {cols[0][i], cols[1][i]}; }
        constexpr const T *                  data() const                                                { return cols[0].data(); }
        LINALG_CONSTEXPR14 T *               data()                                                      { return cols[0].data(); }

        template<class U, class=detail::conv_t<mat,U>> constexpr mat(const U & u)                        : mat(detail::convert<mat>(u)) {}
        template<class U, class=detail::conv_t<U,mat>> constexpr operator U () const                     { return detail::convert<U>(*this); }
    };
    template<class T, int M> struct mat<T,M,3>
    {
        using V=vec<T,M>;
        V                                    cols[3];
        constexpr                            mat()                                                       : cols{} {}
        constexpr                            mat(const V & x_, const V & y_, const V & z_)               : cols{x_, y_, z_} {}
        constexpr explicit                   mat(const T & s)                                            : cols{V(s), V(s), V(s)} {}
        template<class U> constexpr explicit mat(const mat<U,M,3> & m)                                   : cols{V(m[0]), V(m[1]), V(m[2])} {}
        constexpr const V &                  operator[] (int j) const                                    { return cols[j]; }
        LINALG_CONSTEXPR14 V &               operator[] (int j)                                          { return cols[j]; }
        constexpr vec<T,3>                   row(int i) const                                            { return {cols[0][i], cols[1][i], cols[2][i]}; }
        constexpr const T *                  data() const                                                { return cols[0].data(); }
        LINALG_CONSTEXPR14 T *               data()                                                      { return cols[0].data(); }
                                                                                                         
        template<class U, class=detail::conv_t<mat,U>> constexpr mat(const U & u)                        : mat(detail::convert<mat>(u)) {}
        template<class U, class=detail::conv_t<U,mat>> constexpr operator U () const                     { return detail::convert<U>(*this); }
    };
    template<class T, int M> struct mat<T,M,4>
    {
        using V=vec<T,M>;
        V                                    cols[4];
        constexpr                            mat()                                                       : cols{} {}
        constexpr                            mat(const V & x_, const V & y_, const V & z_, const V & w_) : cols{x_, y_, z_, w_} {}
        constexpr explicit                   mat(const T & s)                                            : cols{V(s), V(s), V(s), V(s)} {}
        template<class U> constexpr explicit mat(const mat<U,M,4> & m)                                   : cols{V(m[0]), V(m[1]), V(m[2]), V(m[3])} {}
        constexpr const V &                  operator[] (int j) const                                    { return cols[j]; }
        LINALG_CONSTEXPR14 V &               operator[] (int j)                                          { return cols[j]; }
        constexpr vec<T,4>                   row(int i) const                                            { return {cols[0][i], cols[1][i], cols[2][i], cols[3][i]}; }
        constexpr const T *                  data() const                                                { return cols[0].data(); }
        LINALG_CONSTEXPR14 T *               data()                                                      { return cols[0].data(); }

        template<class U, class=detail::conv_t<mat,U>> constexpr mat(const U & u)                        : mat(detail::convert<mat>(u)) {}
        template<class U, class=detail::conv_t<U,mat>> constexpr operator U () const                     { return detail::convert<U>(*this); }
    };

    ////////////////////////////
    // quat<T> implementation //
    ////////////////////////////

    template<class T> struct quat
    {
        T x,y,z,w;
        constexpr                            quat()                                                       : x(), y(), z(), w() {}
        constexpr                            quat(const T & x_, const T & y_, const T & z_, const T & w_) : x(x_), y(y_), z(z_), w(w_) {}
        constexpr                            quat(const vec<T,3> & xyz, const T & w_)                     : quat(xyz[0], xyz[1], xyz[2], w_) {}
        constexpr explicit                   quat(const vec<T,4> & xyzw)                                  : quat(xyzw[0], xyzw[1], xyzw[2], xyzw[3]) {}
        template<class U> constexpr explicit quat(const quat<U> & q)                                      : quat(static_cast<T>(q.x), static_cast<T>(q.y), static_cast<T>(q.z), static_cast<T>(q.w)) {}
        constexpr vec<T,3>                   xyz() const                                                  { return {x,y,z}; }

        template<class U, class=detail::conv_t<quat,U>> constexpr quat(const U & u)                       : quat(detail::convert<quat>(u)) {}
        template<class U, class=detail::conv_t<U,quat>> constexpr operator U () const                     { return detail::convert<U>(*this); }
    };

    //////////////////////////
    // Relational operators //
    //////////////////////////

    template<class A, class B> using compare_t = typename detail::any_compare<detail::unpack_t<A>, detail::unpack_t<B>>::type;
    template<class A, class B> constexpr compare_t<A,B> compare(const A & a, const B & b) { return detail::any_compare<detail::unpack_t<A>, detail::unpack_t<B>>()(a,b); }

    template<class A, class B> constexpr auto operator == (const A & a, const B & b) -> decltype(compare(a,b) == 0) { return compare(a,b) == 0; }
    template<class A, class B> constexpr auto operator != (const A & a, const B & b) -> decltype(compare(a,b) != 0) { return compare(a,b) != 0; }
    template<class A, class B> constexpr auto operator <  (const A & a, const B & b) -> decltype(compare(a,b) <  0) { return compare(a,b) <  0; }
    template<class A, class B> constexpr auto operator >  (const A & a, const B & b) -> decltype(compare(a,b) >  0) { return compare(a,b) >  0; }
    template<class A, class B> constexpr auto operator <= (const A & a, const B & b) -> decltype(compare(a,b) <= 0) { return compare(a,b) <= 0; }
    template<class A, class B> constexpr auto operator >= (const A & a, const B & b) -> decltype(compare(a,b) >= 0) { return compare(a,b) >= 0; }

    ///////////////////
    // Range support //
    ///////////////////

    template<class T, int M> constexpr const T * begin(const vec<T,M> & a) { return a.data(); }
    template<class T, int M> constexpr const T * end  (const vec<T,M> & a) { return a.data() + M; }
    template<class T, int M, int N> constexpr const vec<T,M> * begin(const mat<T,M,N> & a) { return a.cols; }
    template<class T, int M, int N> constexpr const vec<T,M> * end  (const mat<T,M,N> & a) { return a.cols + N; }

    template<class T, int M> constexpr T * begin(vec<T,M> & a) { return a.data(); }
    template<class T, int M> constexpr T * end  (vec<T,M> & a) { return a.data() + M; }
    template<class T, int M, int N> constexpr vec<T,M> * begin(mat<T,M,N> & a) { return a.cols; }
    template<class T, int M, int N> constexpr vec<T,M> * end  (mat<T,M,N> & a) { return a.cols + N; }

    ////////////////////////////
    // Higher-order functions //
    ////////////////////////////

    // Produce a scalar by applying f(T,T) -> T to adjacent pairs of elements from vector/matrix a in left-to-right order (matching the associativity of arithmetic and logical operators)
    template<class T, class F> constexpr T fold(const vec<T,1> & a, F f) { return a[0]; }
    template<class T, class F> constexpr T fold(const vec<T,2> & a, F f) { return f(a[0],a[1]); }
    template<class T, class F> constexpr T fold(const vec<T,3> & a, F f) { return f(f(a[0],a[1]),a[2]); }
    template<class T, class F> constexpr T fold(const vec<T,4> & a, F f) { return f(f(f(a[0],a[1]),a[2]),a[3]); }

    template<class T, int M, class F> constexpr T fold(const mat<T,M,1> & a, F f) { return fold(a[0],f); }
    template<class T, int M, class F> constexpr T fold(const mat<T,M,2> & a, F f) { return f(fold(a[0],f),fold(a[1],f)); }
    template<class T, int M, class F> constexpr T fold(const mat<T,M,3> & a, F f) { return f(f(fold(a[0],f),fold(a[1],f)),fold(a[2],f)); }
    template<class T, int M, class F> constexpr T fold(const mat<T,M,4> & a, F f) { return f(f(f(fold(a[0],f),fold(a[1],f)),fold(a[2],f)),fold(a[3],f)); }   

    template<class T, class F> constexpr T fold(const quat<T> & a, F f) { return f(f(f(a.x,a.y),a.z),a.w); }

    // Type aliases for the result of calling apply(...) with various arguments, can be used with return type SFINAE to constrian overload sets
    template<class F, class... A> using apply_t = typename detail::any_apply<F,detail::unpack_t<A>...>::type;
    template<class F, class... A> using vec_apply_t = typename detail::vec_apply<F,void,detail::unpack_t<A>...>::type;
    template<class F, class... A> using axa_apply_t = typename detail::axa_apply<F,void,detail::unpack_t<A>...>::type;
    template<class F, class... A> using axs_apply_t = typename detail::axs_apply<F,void,detail::unpack_t<A>...>::type;
    template<class F, class... A> using sxa_apply_t = typename detail::sxa_apply<F,void,detail::unpack_t<A>...>::type;

    // apply(f,...) applies the provided function in an elementwise fashion to its arguments, producing an object of the same dimensions
    template<class F, class... A> constexpr apply_t<F,A...> apply(F func, const A & ... args) { return detail::any_apply<F,detail::unpack_t<A>...>::impl(detail::make_seq<detail::any_apply<F,detail::unpack_t<A>...>::size>{}, func, args...); }

    // map(a,f) is equivalent to apply(f,a)
    template<class A, class F> constexpr apply_t<F,A> map(const A & a, F func) { return apply(func, a); }

    // zip(a,b,f) is equivalent to apply(f,a,b)
    template<class A, class B, class F> constexpr apply_t<F,A,B> zip(const A & a, const B & b, F func) { return apply(func, a, b); }

    ////////////////////////////////////
    // Vector operators and functions //
    ////////////////////////////////////

    // Unary and binary + and - are defined component-wise for all linalg types
    template<class A> constexpr apply_t<detail::op_pos, A> operator + (const A & a) { return apply(detail::op_pos{}, a); }
    template<class A> constexpr apply_t<detail::op_neg, A> operator - (const A & a) { return apply(detail::op_neg{}, a); }

    // Remaining operators are defined componentwise for vector and scalar types
    template<class A> constexpr vec_apply_t<detail::op_cmp, A> operator ~ (const A & a) { return apply(detail::op_cmp{}, a); }
    template<class A> constexpr vec_apply_t<detail::op_not, A> operator ! (const A & a) { return apply(detail::op_not{}, a); }
    template<class A, class B> constexpr vec_apply_t<detail::op_add, A, B> operator +  (const A & a, const B & b) { return apply(detail::op_add{}, a, b); }
    template<class A, class B> constexpr vec_apply_t<detail::op_sub, A, B> operator -  (const A & a, const B & b) { return apply(detail::op_sub{}, a, b); }
    template<class A, class B> constexpr vec_apply_t<detail::op_mul, A, B> operator *  (const A & a, const B & b) { return apply(detail::op_mul{}, a, b); }
    template<class A, class B> constexpr vec_apply_t<detail::op_div, A, B> operator /  (const A & a, const B & b) { return apply(detail::op_div{}, a, b); }
    template<class A, class B> constexpr vec_apply_t<detail::op_mod, A, B> operator %  (const A & a, const B & b) { return apply(detail::op_mod{}, a, b); }
    template<class A, class B> constexpr vec_apply_t<detail::op_un,  A, B> operator |  (const A & a, const B & b) { return apply(detail::op_un{},  a, b); }
    template<class A, class B> constexpr vec_apply_t<detail::op_xor, A, B> operator ^  (const A & a, const B & b) { return apply(detail::op_xor{}, a, b); }
    template<class A, class B> constexpr vec_apply_t<detail::op_int, A, B> operator &  (const A & a, const B & b) { return apply(detail::op_int{}, a, b); }
    template<class A, class B> constexpr vec_apply_t<detail::op_lsh, A, B> operator << (const A & a, const B & b) { return apply(detail::op_lsh{}, a, b); }
    template<class A, class B> constexpr vec_apply_t<detail::op_rsh, A, B> operator >> (const A & a, const B & b) { return apply(detail::op_rsh{}, a, b); }

    // Addition, substraction, and multiplication-by-scalar are defined componentwise for matrix and quaternion types
    template<class A, class B> constexpr axa_apply_t<detail::op_add, A, B> operator + (const A & a, const B & b) { return apply(detail::op_add{}, a, b); }
    template<class A, class B> constexpr axa_apply_t<detail::op_sub, A, B> operator - (const A & a, const B & b) { return apply(detail::op_sub{}, a, b); }
    template<class A, class B> constexpr axs_apply_t<detail::op_mul, A, B> operator * (const A & a, const B & b) { return apply(detail::op_mul{}, a, b); }
    template<class A, class B> constexpr sxa_apply_t<detail::op_mul, A, B> operator * (const A & a, const B & b) { return apply(detail::op_mul{}, a, b); }
    template<class A, class B> constexpr axs_apply_t<detail::op_div, A, B> operator / (const A & a, const B & b) { return apply(detail::op_div{}, a, b); }

    // Binary assignment operators a $= b is always defined as though it were explicitly written a = a $ b
    template<class A, class B> constexpr auto operator +=  (A & a, const B & b) -> decltype(a = a + b) { return a = a + b; }
    template<class A, class B> constexpr auto operator -=  (A & a, const B & b) -> decltype(a = a - b) { return a = a - b; }
    template<class A, class B> constexpr auto operator *=  (A & a, const B & b) -> decltype(a = a * b) { return a = a * b; }
    template<class A, class B> constexpr auto operator /=  (A & a, const B & b) -> decltype(a = a / b) { return a = a / b; }
    template<class A, class B> constexpr auto operator %=  (A & a, const B & b) -> decltype(a = a % b) { return a = a % b; }
    template<class A, class B> constexpr auto operator |=  (A & a, const B & b) -> decltype(a = a | b) { return a = a | b; }
    template<class A, class B> constexpr auto operator ^=  (A & a, const B & b) -> decltype(a = a ^ b) { return a = a ^ b; }
    template<class A, class B> constexpr auto operator &=  (A & a, const B & b) -> decltype(a = a & b) { return a = a & b; }
    template<class A, class B> constexpr auto operator <<= (A & a, const B & b) -> decltype(a = a << b) { return a = a << b; }
    template<class A, class B> constexpr auto operator >>= (A & a, const B & b) -> decltype(a = a >> b) { return a = a >> b; }

    // Reduction functions on vectors
    template<class T, int M> constexpr bool any (const vec<T,M> & a) { return fold(a, detail::op_or{}); }
    template<class T, int M> constexpr bool all (const vec<T,M> & a) { return fold(a, detail::op_and{}); }
    template<class T, int M> constexpr T sum    (const vec<T,M> & a) { return fold(a, detail::op_add{}); }
    template<class T, int M> constexpr T product(const vec<T,M> & a) { return fold(a, detail::op_mul{}); }
    template<class T, int M> constexpr T minelem(const vec<T,M> & a) { return fold(a, detail::min{}); }
    template<class T, int M> constexpr T maxelem(const vec<T,M> & a) { return fold(a, detail::max{}); }

    // Search functions on vectors
    template<class T, int M> int argmin(const vec<T,M> & a) { int j=0; for(int i=1; i<M; ++i) if(a[i] < a[j]) j = i; return j; }
    template<class T, int M> int argmax(const vec<T,M> & a) { int j=0; for(int i=1; i<M; ++i) if(a[i] > a[j]) j = i; return j; }

    // Component-wise standard library math functions on vectors
    template<class A> vec_apply_t<detail::std_abs,   A> abs  (const A & a) { return apply(detail::std_abs{},   a); }
    template<class A> vec_apply_t<detail::std_floor, A> floor(const A & a) { return apply(detail::std_floor{}, a); }
    template<class A> vec_apply_t<detail::std_ceil,  A> ceil (const A & a) { return apply(detail::std_ceil{},  a); }
    template<class A> vec_apply_t<detail::std_exp,   A> exp  (const A & a) { return apply(detail::std_exp{},   a); }
    template<class A> vec_apply_t<detail::std_log,   A> log  (const A & a) { return apply(detail::std_log{},   a); }
    template<class A> vec_apply_t<detail::std_log10, A> log10(const A & a) { return apply(detail::std_log10{}, a); }
    template<class A> vec_apply_t<detail::std_sqrt,  A> sqrt (const A & a) { return apply(detail::std_sqrt{},  a); }
    template<class A> vec_apply_t<detail::std_sin,   A> sin  (const A & a) { return apply(detail::std_sin{},   a); }
    template<class A> vec_apply_t<detail::std_cos,   A> cos  (const A & a) { return apply(detail::std_cos{},   a); }
    template<class A> vec_apply_t<detail::std_tan,   A> tan  (const A & a) { return apply(detail::std_tan{},   a); }
    template<class A> vec_apply_t<detail::std_asin,  A> asin (const A & a) { return apply(detail::std_asin{},  a); }
    template<class A> vec_apply_t<detail::std_acos,  A> acos (const A & a) { return apply(detail::std_acos{},  a); }
    template<class A> vec_apply_t<detail::std_atan,  A> atan (const A & a) { return apply(detail::std_atan{},  a); }
    template<class A> vec_apply_t<detail::std_sinh,  A> sinh (const A & a) { return apply(detail::std_sinh{},  a); }
    template<class A> vec_apply_t<detail::std_cosh,  A> cosh (const A & a) { return apply(detail::std_cosh{},  a); }
    template<class A> vec_apply_t<detail::std_tanh,  A> tanh (const A & a) { return apply(detail::std_tanh{},  a); }
    template<class A> vec_apply_t<detail::std_round, A> round(const A & a) { return apply(detail::std_round{}, a); }

    template<class A, class B> vec_apply_t<detail::std_fmod,     A, B> fmod    (const A & a, const B & b) { return apply(detail::std_fmod{},     a, b); }
    template<class A, class B> vec_apply_t<detail::std_pow,      A, B> pow     (const A & a, const B & b) { return apply(detail::std_pow{},      a, b); }
    template<class A, class B> vec_apply_t<detail::std_atan2,    A, B> atan2   (const A & a, const B & b) { return apply(detail::std_atan2{},    a, b); }
    template<class A, class B> vec_apply_t<detail::std_copysign, A, B> copysign(const A & a, const B & b) { return apply(detail::std_copysign{}, a, b); }

    // Component-wise relational functions on vectors
    template<class A, class B> constexpr vec_apply_t<detail::op_eq, A, B> equal  (const A & a, const B & b) { return apply(detail::op_eq{}, a, b); }
    template<class A, class B> constexpr vec_apply_t<detail::op_ne, A, B> nequal (const A & a, const B & b) { return apply(detail::op_ne{}, a, b); }
    template<class A, class B> constexpr vec_apply_t<detail::op_lt, A, B> less   (const A & a, const B & b) { return apply(detail::op_lt{}, a, b); }
    template<class A, class B> constexpr vec_apply_t<detail::op_gt, A, B> greater(const A & a, const B & b) { return apply(detail::op_gt{}, a, b); }
    template<class A, class B> constexpr vec_apply_t<detail::op_le, A, B> lequal (const A & a, const B & b) { return apply(detail::op_le{}, a, b); }
    template<class A, class B> constexpr vec_apply_t<detail::op_ge, A, B> gequal (const A & a, const B & b) { return apply(detail::op_ge{}, a, b); }

    // Component-wise selection functions on vectors
    template<class A, class B> constexpr vec_apply_t<detail::min, A, B> min(const A & a, const B & b) { return apply(detail::min{}, a, b); }
    template<class A, class B> constexpr vec_apply_t<detail::max, A, B> max(const A & a, const B & b) { return apply(detail::max{}, a, b); }
    template<class X, class L, class H> constexpr vec_apply_t<detail::clamp,  X, L, H> clamp (const X & x, const L & l, const H & h) { return apply(detail::clamp{},  x, l, h); }
    template<class P, class A, class B> constexpr vec_apply_t<detail::select, P, A, B> select(const P & p, const A & a, const B & b) { return apply(detail::select{}, p, a, b); }
    template<class A, class B, class T> constexpr vec_apply_t<detail::lerp,   A, B, T> lerp  (const A & a, const B & b, const T & t) { return apply(detail::lerp{},   a, b, t); }

    // Vector algebra functions
    template<class T> constexpr T               cross    (const vec<T,2> & a, const vec<T,2> & b)      { return a[0]*b[1]-a[1]*b[0]; }
    template<class T> constexpr vec<T,2>        cross    (T a, const vec<T,2> & b)                     { return {-a*b[1], a*b[0]}; }
    template<class T> constexpr vec<T,2>        cross    (const vec<T,2> & a, T b)                     { return {a[1]*b, -a[0]*b}; }
    template<class T> constexpr vec<T,3>        cross    (const vec<T,3> & a, const vec<T,3> & b)      { return {a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]}; }
    template<class T, int M> constexpr T        dot      (const vec<T,M> & a, const vec<T,M> & b)      { return sum(a*b); }
    template<class T, int M> constexpr T        length2  (const vec<T,M> & a)                          { return dot(a,a); }
    template<class T, int M> T                  length   (const vec<T,M> & a)                          { return std::sqrt(length2(a)); }
    template<class T, int M> vec<T,M>           normalize(const vec<T,M> & a)                          { return a / length(a); }
    template<class T, int M> constexpr T        distance2(const vec<T,M> & a, const vec<T,M> & b)      { return length2(b-a); }
    template<class T, int M> T                  distance (const vec<T,M> & a, const vec<T,M> & b)      { return length(b-a); }
    template<class T, int M> T                  uangle   (const vec<T,M> & a, const vec<T,M> & b)      { T d=dot(a,b); return d > 1 ? 0 : std::acos(d < -1 ? -1 : d); }
    template<class T, int M> T                  angle    (const vec<T,M> & a, const vec<T,M> & b)      { return uangle(normalize(a), normalize(b)); }
    template<class T> vec<T,2>                  rot      (T a, const vec<T,2> & v)                     { const T s = std::sin(a), c = std::cos(a); return {v[0]*c - v[1]*s, v[0]*s + v[1]*c}; }
    template<class T, int M> vec<T,M>           nlerp    (const vec<T,M> & a, const vec<T,M> & b, T t) { return normalize(lerp(a,b,t)); }
    template<class T, int M> vec<T,M>           slerp    (const vec<T,M> & a, const vec<T,M> & b, T t) { T th=uangle(a,b); return th == 0 ? a : a*(std::sin(th*(1-t))/std::sin(th)) + b*(std::sin(th*t)/std::sin(th)); }

    ////////////////////////////////////
    // Matrix operators and functions //
    ////////////////////////////////////

    // Matrix product
    template<class T, int M> constexpr vec<T,M> operator * (const mat<T,M,1> & a, const vec<T,1> & b) { return a[0]*b[0]; }
    template<class T, int M> constexpr vec<T,M> operator * (const mat<T,M,2> & a, const vec<T,2> & b) { return a[0]*b[0] + a[1]*b[1]; }
    template<class T, int M> constexpr vec<T,M> operator * (const mat<T,M,3> & a, const vec<T,3> & b) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }
    template<class T, int M> constexpr vec<T,M> operator * (const mat<T,M,4> & a, const vec<T,4> & b) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3]; }
    template<class T, int M, int N> constexpr mat<T,M,1> operator * (const mat<T,M,N> & a, const mat<T,N,1> & b) { return {a*b[0]}; }
    template<class T, int M, int N> constexpr mat<T,M,2> operator * (const mat<T,M,N> & a, const mat<T,N,2> & b) { return {a*b[0], a*b[1]}; }
    template<class T, int M, int N> constexpr mat<T,M,3> operator * (const mat<T,M,N> & a, const mat<T,N,3> & b) { return {a*b[0], a*b[1], a*b[2]}; }
    template<class T, int M, int N> constexpr mat<T,M,4> operator * (const mat<T,M,N> & a, const mat<T,N,4> & b) { return {a*b[0], a*b[1], a*b[2], a*b[3]}; }

    // Matrix algebra functions
    template<class T> constexpr vec<T,1>          diagonal    (const mat<T,1,1> & a) { return {a[0][0]}; }
    template<class T> constexpr vec<T,2>          diagonal    (const mat<T,2,2> & a) { return {a[0][0], a[1][1]}; }
    template<class T> constexpr vec<T,3>          diagonal    (const mat<T,3,3> & a) { return {a[0][0], a[1][1], a[2][2]}; }
    template<class T> constexpr vec<T,4>          diagonal    (const mat<T,4,4> & a) { return {a[0][0], a[1][1], a[2][2], a[3][3]}; }
    template<class T, int M> constexpr mat<T,M,1> outerprod   (const vec<T,M> & a, const vec<T,1> & b) { return {a*b[0]}; }
    template<class T, int M> constexpr mat<T,M,2> outerprod   (const vec<T,M> & a, const vec<T,2> & b) { return {a*b[0], a*b[1]}; }
    template<class T, int M> constexpr mat<T,M,3> outerprod   (const vec<T,M> & a, const vec<T,3> & b) { return {a*b[0], a*b[1], a*b[2]}; }
    template<class T, int M> constexpr mat<T,M,4> outerprod   (const vec<T,M> & a, const vec<T,4> & b) { return {a*b[0], a*b[1], a*b[2], a*b[3]}; }
    template<class T, int M> constexpr mat<T,M,1> transpose   (const mat<T,1,M> & m) { return {m.row(0)}; }
    template<class T, int M> constexpr mat<T,M,2> transpose   (const mat<T,2,M> & m) { return {m.row(0), m.row(1)}; }
    template<class T, int M> constexpr mat<T,M,3> transpose   (const mat<T,3,M> & m) { return {m.row(0), m.row(1), m.row(2)}; }
    template<class T, int M> constexpr mat<T,M,4> transpose   (const mat<T,4,M> & m) { return {m.row(0), m.row(1), m.row(2), m.row(3)}; }
    template<class T> constexpr mat<T,1,1>        adjugate    (const mat<T,1,1> & a) { return {vec<T,1>{1}}; }
    template<class T> constexpr mat<T,2,2>        adjugate    (const mat<T,2,2> & a) { return {{a[1][1], -a[0][1]}, {-a[1][0], a[0][0]}}; }
    template<class T> constexpr mat<T,3,3>        adjugate    (const mat<T,3,3> & a);
    template<class T> constexpr mat<T,4,4>        adjugate    (const mat<T,4,4> & a);
    template<class T> constexpr T                 determinant (const mat<T,1,1> & a) { return a[0][0]; }
    template<class T> constexpr T                 determinant (const mat<T,2,2> & a) { return a[0][0]*a[1][1] - a[0][1]*a[1][0]; }
    template<class T> constexpr T                 determinant (const mat<T,3,3> & a) { return a[0][0]*(a[1][1]*a[2][2] - a[2][1]*a[1][2]) + a[0][1]*(a[1][2]*a[2][0] - a[2][2]*a[1][0]) + a[0][2]*(a[1][0]*a[2][1] - a[2][0]*a[1][1]); }
    template<class T> constexpr T                 determinant (const mat<T,4,4> & a);
    template<class T, int N> constexpr T          trace       (const mat<T,N,N> & a) { return sum(diagonal(a)); }
    template<class T, int N> constexpr mat<T,N,N> inverse     (const mat<T,N,N> & a) { return adjugate(a)/determinant(a); }

    ////////////////////////////////////////
    // Quaternion operators and functions //
    ////////////////////////////////////////

    // Quaternion product
    template<class T> constexpr quat<T> operator * (const quat<T> & a, const quat<T> & b) { return {a.x*b.w+a.w*b.x+a.y*b.z-a.z*b.y, a.y*b.w+a.w*b.y+a.z*b.x-a.x*b.z, a.z*b.w+a.w*b.z+a.x*b.y-a.y*b.x, a.w*b.w-a.x*b.x-a.y*b.y-a.z*b.z}; }

    // Quaternion algebra functions
    template<class T> quat<T>           exp      (const quat<T> & q)                         { const auto v = q.xyz(); const auto vv = length(v); return std::exp(q.w)*quat<T>{v * (vv > 0 ? std::sin(vv)/vv : 0), std::cos(vv)}; }
    template<class T> quat<T>           log      (const quat<T> & q)                         { const auto v = q.xyz(); const auto vv = length(v), qq = length(q); return {v * (vv > 0 ? std::acos(q.w/qq)/vv : 0), std::log(qq)}; }
    template<class T> quat<T>           pow      (const quat<T> & q, const T & p)            { const auto v = q.xyz(); const auto vv = length(v), qq = length(q), th = std::acos(q.w/qq); return std::pow(qq,p)*quat<T>{v * (vv > 0 ? std::sin(p*th)/vv : 0), std::cos(p*th)}; }
    template<class T> constexpr quat<T> conjugate(const quat<T> & a)                         { return {-a.x, -a.y, -a.z, a.w}; }
    template<class T> constexpr T       dot      (const quat<T> & a, const quat<T> & b)      { return a.x*b.x + a.y*b.y + a.z*b.z + a.w*b.w; }
    template<class T> constexpr T       length2  (const quat<T> & a)                         { return dot(a,a); }
    template<class T> T                 length   (const quat<T> & a)                         { return std::sqrt(length2(a)); }
    template<class T> constexpr quat<T> inverse  (const quat<T> & a)                         { return conjugate(a) / length2(a); }
    template<class T> quat<T>           normalize(const quat<T> & a)                         { return a / length(a); }
    template<class T> T                 uangle   (const quat<T> & a, const quat<T> & b)      { T d=dot(a,b); return d > 1 ? 0 : std::acos(d < -1 ? -1 : d); }
    template<class T> constexpr quat<T> lerp     (const quat<T> & a, const quat<T> & b, T t) { return a*(1-t) + b*t; }
    template<class T> quat<T>           nlerp    (const quat<T> & a, const quat<T> & b, T t) { return normalize(lerp(a,b,t)); }
    template<class T> quat<T>           slerp    (const quat<T> & a, const quat<T> & b, T t) { T th=uangle(a,b); return th == 0 ? a : a*(std::sin(th*(1-t))/std::sin(th)) + b*(std::sin(th*t)/std::sin(th)); }

    // Support for 3D spatial rotations using quaternions, via q * v * conjugate(q)
    template<class T> constexpr vec<T,3>   qxdir (const quat<T> & q)                          { return {q.w*q.w+q.x*q.x-q.y*q.y-q.z*q.z, (q.x*q.y+q.z*q.w)*2, (q.z*q.x-q.y*q.w)*2}; }
    template<class T> constexpr vec<T,3>   qydir (const quat<T> & q)                          { return {(q.x*q.y-q.z*q.w)*2, q.w*q.w-q.x*q.x+q.y*q.y-q.z*q.z, (q.y*q.z+q.x*q.w)*2}; }
    template<class T> constexpr vec<T,3>   qzdir (const quat<T> & q)                          { return {(q.z*q.x+q.y*q.w)*2, (q.y*q.z-q.x*q.w)*2, q.w*q.w-q.x*q.x-q.y*q.y+q.z*q.z}; }
    template<class T> constexpr mat<T,3,3> qmat  (const quat<T> & q)                          { return {qxdir(q), qydir(q), qzdir(q)}; }
    template<class T> constexpr vec<T,3>   qrot  (const quat<T> & q, const vec<T,3> & v)      { return qxdir(q)*v.x + qydir(q)*v.y + qzdir(q)*v.z; }
    template<class T> T                    qangle(const quat<T> & q)                          { return std::atan2(length(q.xyz()), q.w)*2; }
    template<class T> vec<T,3>             qaxis (const quat<T> & q)                          { return normalize(q.xyz()); }

    //////////////////////////////////////////////////////////////////////////////////
    // Standard type aliases, can be brought in via using namespace linalg::aliases //
    //////////////////////////////////////////////////////////////////////////////////

    namespace aliases
    {
        using bool1=vec<bool,1>; using byte1=vec<uint8_t,1>; using short1=vec<int16_t,1>; using ushort1=vec<uint16_t,1>;
        using bool2=vec<bool,2>; using byte2=vec<uint8_t,2>; using short2=vec<int16_t,2>; using ushort2=vec<uint16_t,2>;
        using bool3=vec<bool,3>; using byte3=vec<uint8_t,3>; using short3=vec<int16_t,3>; using ushort3=vec<uint16_t,3>; 
        using bool4=vec<bool,4>; using byte4=vec<uint8_t,4>; using short4=vec<int16_t,4>; using ushort4=vec<uint16_t,4>;
        using int1=vec<int,1>; using uint1=vec<unsigned,1>; using float1=vec<float,1>; using double1=vec<double,1>;
        using int2=vec<int,2>; using uint2=vec<unsigned,2>; using float2=vec<float,2>; using double2=vec<double,2>;
        using int3=vec<int,3>; using uint3=vec<unsigned,3>; using float3=vec<float,3>; using double3=vec<double,3>;
        using int4=vec<int,4>; using uint4=vec<unsigned,4>; using float4=vec<float,4>; using double4=vec<double,4>;
        using bool1x1=mat<bool,1,1>; using int1x1=mat<int,1,1>; using float1x1=mat<float,1,1>; using double1x1=mat<double,1,1>;
        using bool1x2=mat<bool,1,2>; using int1x2=mat<int,1,2>; using float1x2=mat<float,1,2>; using double1x2=mat<double,1,2>;
        using bool1x3=mat<bool,1,3>; using int1x3=mat<int,1,3>; using float1x3=mat<float,1,3>; using double1x3=mat<double,1,3>;
        using bool1x4=mat<bool,1,4>; using int1x4=mat<int,1,4>; using float1x4=mat<float,1,4>; using double1x4=mat<double,1,4>;
        using bool2x1=mat<bool,2,1>; using int2x1=mat<int,2,1>; using float2x1=mat<float,2,1>; using double2x1=mat<double,2,1>;
        using bool2x2=mat<bool,2,2>; using int2x2=mat<int,2,2>; using float2x2=mat<float,2,2>; using double2x2=mat<double,2,2>;
        using bool2x3=mat<bool,2,3>; using int2x3=mat<int,2,3>; using float2x3=mat<float,2,3>; using double2x3=mat<double,2,3>;
        using bool2x4=mat<bool,2,4>; using int2x4=mat<int,2,4>; using float2x4=mat<float,2,4>; using double2x4=mat<double,2,4>;
        using bool3x1=mat<bool,3,1>; using int3x1=mat<int,3,1>; using float3x1=mat<float,3,1>; using double3x1=mat<double,3,1>;
        using bool3x2=mat<bool,3,2>; using int3x2=mat<int,3,2>; using float3x2=mat<float,3,2>; using double3x2=mat<double,3,2>;
        using bool3x3=mat<bool,3,3>; using int3x3=mat<int,3,3>; using float3x3=mat<float,3,3>; using double3x3=mat<double,3,3>;
        using bool3x4=mat<bool,3,4>; using int3x4=mat<int,3,4>; using float3x4=mat<float,3,4>; using double3x4=mat<double,3,4>;
        using bool4x1=mat<bool,4,1>; using int4x1=mat<int,4,1>; using float4x1=mat<float,4,1>; using double4x1=mat<double,4,1>;
        using bool4x2=mat<bool,4,2>; using int4x2=mat<int,4,2>; using float4x2=mat<float,4,2>; using double4x2=mat<double,4,2>;
        using bool4x3=mat<bool,4,3>; using int4x3=mat<int,4,3>; using float4x3=mat<float,4,3>; using double4x3=mat<double,4,3>;
        using bool4x4=mat<bool,4,4>; using int4x4=mat<int,4,4>; using float4x4=mat<float,4,4>; using double4x4=mat<double,4,4>;
        using quatf=quat<float>; using quatd=quat<double>;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Reasonable default printing behavior, can be brought in via using namespace linalg::ostream_overloads //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////

    namespace ostream_overloads
    {
        // These overloads stream out something that resembles an aggregate literal that could be used to construct the specified value
        template<class C, class T> std::basic_ostream<C> & operator << (std::basic_ostream<C> & out, const vec<T,1> & v) { return out << '{' << v[0] << '}'; }
        template<class C, class T> std::basic_ostream<C> & operator << (std::basic_ostream<C> & out, const vec<T,2> & v) { return out << '{' << v[0] << ',' << v[1] << '}'; }
        template<class C, class T> std::basic_ostream<C> & operator << (std::basic_ostream<C> & out, const vec<T,3> & v) { return out << '{' << v[0] << ',' << v[1] << ',' << v[2] << '}'; }
        template<class C, class T> std::basic_ostream<C> & operator << (std::basic_ostream<C> & out, const vec<T,4> & v) { return out << '{' << v[0] << ',' << v[1] << ',' << v[2] << ',' << v[3] << '}'; }
        template<class C, class T, int M> std::basic_ostream<C> & operator << (std::basic_ostream<C> & out, const mat<T,M,1> & m) { return out << '{' << m[0] << '}'; }
        template<class C, class T, int M> std::basic_ostream<C> & operator << (std::basic_ostream<C> & out, const mat<T,M,2> & m) { return out << '{' << m[0] << ',' << m[1] << '}'; }
        template<class C, class T, int M> std::basic_ostream<C> & operator << (std::basic_ostream<C> & out, const mat<T,M,3> & m) { return out << '{' << m[0] << ',' << m[1] << ',' << m[2] << '}'; }
        template<class C, class T, int M> std::basic_ostream<C> & operator << (std::basic_ostream<C> & out, const mat<T,M,4> & m) { return out << '{' << m[0] << ',' << m[1] << ',' << m[2] << ',' << m[3] << '}'; }
        template<class C, class T> std::basic_ostream<C> & operator << (std::basic_ostream<C> & out, const quat<T> & q) { return out << '{' << q.x << ',' << q.y << ',' << q.z << ',' << q.w << '}'; }
        template<class C, class T, class A, int... I> std::basic_ostream<C> & operator << (std::basic_ostream<C> & out, const _lswizzle<T,A,I...> & s) { vec<T,sizeof...(I)> v = s; return out << v; }
        template<class C, class T, class A, int... I> std::basic_ostream<C> & operator << (std::basic_ostream<C> & out, const _rswizzle<T,A,I...> & s) { vec<T,sizeof...(I)> v = s; return out << v; }
    }

    ///////////////////////////////////////////////////////////////////////////
    // Chopping block - Functions in this section may be reworked or removed //
    ///////////////////////////////////////////////////////////////////////////

    // Factory functions for 3D spatial transformations (will possibly be removed or changed in a future version)
    enum fwd_axis { neg_z, pos_z };                 // Should projection matrices be generated assuming forward is {0,0,-1} or {0,0,1}
    enum z_range { neg_one_to_one, zero_to_one };   // Should projection matrices map z into the range of [-1,1] or [0,1]?
    template<class T> quat<T>              rotation_quat     (const vec<T,3> & axis, T angle)             { return {axis*std::sin(angle/2), std::cos(angle/2)}; }
    template<class T> quat<T>              rotation_quat     (const vec<T,3> & from, const vec<T,3> & to) { return rotation_quat(normalize(cross(from,to)), angle(from,to)); }
    template<class T> quat<T>              rotation_quat     (const mat<T,3,3> & m);
    template<class T> constexpr mat<T,4,4> translation_matrix(const vec<T,3> & translation)          { return {{1,0,0,0},{0,1,0,0},{0,0,1,0},{translation,1}}; }
    template<class T> constexpr mat<T,4,4> rotation_matrix   (const quat<T> & rotation)              { return {{qxdir(rotation),0}, {qydir(rotation),0}, {qzdir(rotation),0}, {0,0,0,1}}; }
    template<class T> constexpr mat<T,4,4> scaling_matrix    (const vec<T,3> & scaling)              { return {{scaling[0],0,0,0}, {0,scaling[1],0,0}, {0,0,scaling[2],0}, {0,0,0,1}}; }
    template<class T> constexpr mat<T,4,4> pose_matrix       (const quat<T> & q, const vec<T,3> & p) { return {{qxdir(q),0}, {qydir(q),0}, {qzdir(q),0}, {p,1}}; }
    template<class T> constexpr mat<T,4,4> frustum_matrix    (T x0, T x1, T y0, T y1, T n, T f, fwd_axis a = neg_z, z_range z = neg_one_to_one) { return {{2*n/(x1-x0),0,0,0}, {0,2*n/(y1-y0),0,0}, vec<T,4>{-(x0+x1)/(x1-x0), -(y0+y1)/(y1-y0), (z == zero_to_one ? f : f+n)/(f-n), 1} * (a == pos_z ? T(1) : T(-1)), {0,0,(z == zero_to_one ? -1 : -2)*n*f/(f-n),0}}; }
    template<class T> mat<T,4,4>           perspective_matrix(T fovy, T aspect, T n, T f, fwd_axis a = neg_z, z_range z = neg_one_to_one)       { T y = n*std::tan(fovy / 2), x = y*aspect; return frustum_matrix(-x, x, -y, y, n, f, a, z); }
}

////////////////////////////////////////////////////////
// Specializations of std::hash<...> for linalg types //
////////////////////////////////////////////////////////

namespace std 
{ 
    template<class T> struct hash<linalg::vec<T,2>> { std::size_t operator()(const linalg::vec<T,2> & v) const { std::hash<T> h; return h(v[0]) ^ (h(v[1]) << 1); } };
    template<class T> struct hash<linalg::vec<T,3>> { std::size_t operator()(const linalg::vec<T,3> & v) const { std::hash<T> h; return h(v[0]) ^ (h(v[1]) << 1) ^ (h(v[2]) << 2); } };
    template<class T> struct hash<linalg::vec<T,4>> { std::size_t operator()(const linalg::vec<T,4> & v) const { std::hash<T> h; return h(v[0]) ^ (h(v[1]) << 1) ^ (h(v[2]) << 2) ^ (h(v[3]) << 3); } };
    template<class T, int M> struct hash<linalg::mat<T,M,2>> { std::size_t operator()(const linalg::mat<T,M,2> & m) const { std::hash<linalg::vec<T,M>> h; return h(m[0]) ^ (h(m[1]) << M); } };
    template<class T, int M> struct hash<linalg::mat<T,M,3>> { std::size_t operator()(const linalg::mat<T,M,3> & m) const { std::hash<linalg::vec<T,M>> h; return h(m[0]) ^ (h(m[1]) << M) ^ (h(m[2]) << (M*2)); } };
    template<class T, int M> struct hash<linalg::mat<T,M,4>> { std::size_t operator()(const linalg::mat<T,M,4> & m) const { std::hash<linalg::vec<T,M>> h; return h(m[0]) ^ (h(m[1]) << M) ^ (h(m[2]) << (M*2)) ^ (h(m[3]) << (M*3)); } };
    template<class T> struct hash<linalg::quat<T>> { std::size_t operator()(const linalg::quat<T> & q) const { std::hash<T> h; return h(q.x) ^ (h(q.y) << 1) ^ (h(q.z) << 2) ^ (h(q.w) << 3); } };
}

////////////////////////////////////////////////////////////
// Definitions of functions too long to be defined inline //
////////////////////////////////////////////////////////////

template<class T> constexpr linalg::mat<T,3,3> linalg::adjugate(const mat<T,3,3> & a) 
{ 
    return {{a[1][1]*a[2][2] - a[2][1]*a[1][2], a[2][1]*a[0][2] - a[0][1]*a[2][2], a[0][1]*a[1][2] - a[1][1]*a[0][2]},
            {a[1][2]*a[2][0] - a[2][2]*a[1][0], a[2][2]*a[0][0] - a[0][2]*a[2][0], a[0][2]*a[1][0] - a[1][2]*a[0][0]},
            {a[1][0]*a[2][1] - a[2][0]*a[1][1], a[2][0]*a[0][1] - a[0][0]*a[2][1], a[0][0]*a[1][1] - a[1][0]*a[0][1]}}; 
}

template<class T> constexpr linalg::mat<T,4,4> linalg::adjugate(const mat<T,4,4> & a) 
{ 
    return {{a[1][1]*a[2][2]*a[3][3] + a[3][1]*a[1][2]*a[2][3] + a[2][1]*a[3][2]*a[1][3] - a[1][1]*a[3][2]*a[2][3] - a[2][1]*a[1][2]*a[3][3] - a[3][1]*a[2][2]*a[1][3],
             a[0][1]*a[3][2]*a[2][3] + a[2][1]*a[0][2]*a[3][3] + a[3][1]*a[2][2]*a[0][3] - a[3][1]*a[0][2]*a[2][3] - a[2][1]*a[3][2]*a[0][3] - a[0][1]*a[2][2]*a[3][3],
             a[0][1]*a[1][2]*a[3][3] + a[3][1]*a[0][2]*a[1][3] + a[1][1]*a[3][2]*a[0][3] - a[0][1]*a[3][2]*a[1][3] - a[1][1]*a[0][2]*a[3][3] - a[3][1]*a[1][2]*a[0][3],
             a[0][1]*a[2][2]*a[1][3] + a[1][1]*a[0][2]*a[2][3] + a[2][1]*a[1][2]*a[0][3] - a[0][1]*a[1][2]*a[2][3] - a[2][1]*a[0][2]*a[1][3] - a[1][1]*a[2][2]*a[0][3]},
            {a[1][2]*a[3][3]*a[2][0] + a[2][2]*a[1][3]*a[3][0] + a[3][2]*a[2][3]*a[1][0] - a[1][2]*a[2][3]*a[3][0] - a[3][2]*a[1][3]*a[2][0] - a[2][2]*a[3][3]*a[1][0],
             a[0][2]*a[2][3]*a[3][0] + a[3][2]*a[0][3]*a[2][0] + a[2][2]*a[3][3]*a[0][0] - a[0][2]*a[3][3]*a[2][0] - a[2][2]*a[0][3]*a[3][0] - a[3][2]*a[2][3]*a[0][0],
             a[0][2]*a[3][3]*a[1][0] + a[1][2]*a[0][3]*a[3][0] + a[3][2]*a[1][3]*a[0][0] - a[0][2]*a[1][3]*a[3][0] - a[3][2]*a[0][3]*a[1][0] - a[1][2]*a[3][3]*a[0][0],
             a[0][2]*a[1][3]*a[2][0] + a[2][2]*a[0][3]*a[1][0] + a[1][2]*a[2][3]*a[0][0] - a[0][2]*a[2][3]*a[1][0] - a[1][2]*a[0][3]*a[2][0] - a[2][2]*a[1][3]*a[0][0]},
            {a[1][3]*a[2][0]*a[3][1] + a[3][3]*a[1][0]*a[2][1] + a[2][3]*a[3][0]*a[1][1] - a[1][3]*a[3][0]*a[2][1] - a[2][3]*a[1][0]*a[3][1] - a[3][3]*a[2][0]*a[1][1],
             a[0][3]*a[3][0]*a[2][1] + a[2][3]*a[0][0]*a[3][1] + a[3][3]*a[2][0]*a[0][1] - a[0][3]*a[2][0]*a[3][1] - a[3][3]*a[0][0]*a[2][1] - a[2][3]*a[3][0]*a[0][1],
             a[0][3]*a[1][0]*a[3][1] + a[3][3]*a[0][0]*a[1][1] + a[1][3]*a[3][0]*a[0][1] - a[0][3]*a[3][0]*a[1][1] - a[1][3]*a[0][0]*a[3][1] - a[3][3]*a[1][0]*a[0][1],
             a[0][3]*a[2][0]*a[1][1] + a[1][3]*a[0][0]*a[2][1] + a[2][3]*a[1][0]*a[0][1] - a[0][3]*a[1][0]*a[2][1] - a[2][3]*a[0][0]*a[1][1] - a[1][3]*a[2][0]*a[0][1]},
            {a[1][0]*a[3][1]*a[2][2] + a[2][0]*a[1][1]*a[3][2] + a[3][0]*a[2][1]*a[1][2] - a[1][0]*a[2][1]*a[3][2] - a[3][0]*a[1][1]*a[2][2] - a[2][0]*a[3][1]*a[1][2],
             a[0][0]*a[2][1]*a[3][2] + a[3][0]*a[0][1]*a[2][2] + a[2][0]*a[3][1]*a[0][2] - a[0][0]*a[3][1]*a[2][2] - a[2][0]*a[0][1]*a[3][2] - a[3][0]*a[2][1]*a[0][2],
             a[0][0]*a[3][1]*a[1][2] + a[1][0]*a[0][1]*a[3][2] + a[3][0]*a[1][1]*a[0][2] - a[0][0]*a[1][1]*a[3][2] - a[3][0]*a[0][1]*a[1][2] - a[1][0]*a[3][1]*a[0][2],
             a[0][0]*a[1][1]*a[2][2] + a[2][0]*a[0][1]*a[1][2] + a[1][0]*a[2][1]*a[0][2] - a[0][0]*a[2][1]*a[1][2] - a[1][0]*a[0][1]*a[2][2] - a[2][0]*a[1][1]*a[0][2]}}; 
}

template<class T> constexpr T linalg::determinant(const mat<T,4,4> & a) 
{ 
    return a[0][0]*(a[1][1]*a[2][2]*a[3][3] + a[3][1]*a[1][2]*a[2][3] + a[2][1]*a[3][2]*a[1][3] - a[1][1]*a[3][2]*a[2][3] - a[2][1]*a[1][2]*a[3][3] - a[3][1]*a[2][2]*a[1][3])
         + a[0][1]*(a[1][2]*a[3][3]*a[2][0] + a[2][2]*a[1][3]*a[3][0] + a[3][2]*a[2][3]*a[1][0] - a[1][2]*a[2][3]*a[3][0] - a[3][2]*a[1][3]*a[2][0] - a[2][2]*a[3][3]*a[1][0])
         + a[0][2]*(a[1][3]*a[2][0]*a[3][1] + a[3][3]*a[1][0]*a[2][1] + a[2][3]*a[3][0]*a[1][1] - a[1][3]*a[3][0]*a[2][1] - a[2][3]*a[1][0]*a[3][1] - a[3][3]*a[2][0]*a[1][1])
         + a[0][3]*(a[1][0]*a[3][1]*a[2][2] + a[2][0]*a[1][1]*a[3][2] + a[3][0]*a[2][1]*a[1][2] - a[1][0]*a[2][1]*a[3][2] - a[3][0]*a[1][1]*a[2][2] - a[2][0]*a[3][1]*a[1][2]); 
}

template<class T> linalg::quat<T> linalg::rotation_quat(const mat<T,3,3> & m)
{
    const vec<T,4> q {m[0][0]-m[1][1]-m[2][2], m[1][1]-m[0][0]-m[2][2], m[2][2]-m[0][0]-m[1][1], m[0][0]+m[1][1]+m[2][2]}, s[] {
        {1, m[0][1] + m[1][0], m[2][0] + m[0][2], m[1][2] - m[2][1]}, 
        {m[0][1] + m[1][0], 1, m[1][2] + m[2][1], m[2][0] - m[0][2]},
        {m[0][2] + m[2][0], m[1][2] + m[2][1], 1, m[0][1] - m[1][0]},
        {m[1][2] - m[2][1], m[2][0] - m[0][2], m[0][1] - m[1][0], 1}};
    return quat<T>{copysign(normalize(sqrt(max(T(0), T(1)+q))), s[argmax(q)])};
}

#endif
