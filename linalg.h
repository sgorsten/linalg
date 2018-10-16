// linalg.h - v2.0 - Single-header public domain linear algebra library
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
#include <type_traits>  // For std::is_arithmetic, std::is_same, std::enable_if, std::conditional
#include <functional>   // For std::hash

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
    template<class T> struct converter<mat<T,2,2>, identity_t> { mat<T,2,2> operator() (identity_t) const { return {{1,0},{0,1}}; } };
    template<class T> struct converter<mat<T,3,3>, identity_t> { mat<T,3,3> operator() (identity_t) const { return {{1,0,0},{0,1,0},{0,0,1}}; } };
    template<class T> struct converter<mat<T,4,4>, identity_t> { mat<T,4,4> operator() (identity_t) const { return {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}}; } };
    template<class T> struct converter<quat<T>, identity_t> { quat<T> operator() (identity_t) const { return {0,0,0,1}; } };
    constexpr identity_t identity {1};

    // Support for named scalar types, as in https://t0rakka.silvrback.com/simd-scalar-accessor
    template<class T, class A, int I> struct scalar_accessor 
    { 
        A elems;
        scalar_accessor()                          = delete;
        operator const T & () const                { return elems[I]; }
        operator T & ()                            { return elems[I]; }
        const T * operator & () const              { return &elems[I]; }
        T * operator & ()                          { return &elems[I]; }
        T & operator = (const T & value)           { return elems[I] = value; }
        T & operator = (const scalar_accessor & r) { return elems[I] = r; } // Default copy operator has incorrect semantics, force interpretation as scalar
    };
    template<class T, class A, int I0, int I1> struct swizzle2
    {
        A                           e;
        constexpr                   swizzle2(T e0, T e1)                    : e{} { e[I0]=e0; e[I1]=e1; } // Aggregate constructor (braced-initialization) has incorrect semantics
                                    swizzle2(const vec<T,2> & r)            : swizzle2(r[0], r[1]) {}
        template<class B, int... J> swizzle2(const swizzle2<T,B,J...> & r)  : swizzle2(vec<T,2>(r)) {}
                                    operator vec<T,2> () const              { return {e[I0], e[I1]}; }           
        swizzle2 &                  operator = (const swizzle2 & r)         { e[I0]=r.e[I0]; e[I1]=r.e[I1]; return *this; } // Default copy operator has incorrect semantics
    };
    template<class T, class A, int I0, int I1, int I2> struct swizzle3
    {
        A                           e;
        constexpr                   swizzle3(T e0, T e1, T e2)              : e{} { e[I0]=e0; e[I1]=e1; e[I2]=e2; } // Aggregate constructor (braced-initialization) has incorrect semantics
                                    swizzle3(const vec<T,3> & r)            : swizzle3(r[0], r[1], r[2]) {}
        template<class B, int... J> swizzle3(const swizzle3<T,B,J...> & r)  : swizzle3(vec<T,3>(r)) {}
                                    operator vec<T,3> () const              { return {e[I0], e[I1], e[I2]}; }           
        swizzle3 &                  operator = (const swizzle3 & r)         { e[I0]=r.e[I0]; e[I1]=r.e[I1]; e[I2]=r.e[I2]; return *this; } // Default copy operator has incorrect semantics
    };
    template<class T, class A, int I0, int I1, int I2, int I3> struct swizzle4
    {
        A                           e;
        constexpr                   swizzle4(T e0, T e1, T e2, T e3)        : e{} { e[I0]=e0; e[I1]=e1; e[I2]=e2; e[I3]=e3; } // Aggregate constructor (braced-initialization) has incorrect semantics
                                    swizzle4(const vec<T,4> & r)            : swizzle4(r[0], r[1], r[2], r[3]) {}
        template<class B, int... J> swizzle4(const swizzle4<T,B,J...> & r)  : swizzle4(vec<T,4>(r)) {}
                                    operator vec<T,4> () const              { return {e[I0], e[I1], e[I2], e[I3]}; }           
        swizzle4 &                  operator = (const swizzle4 & r)         { e[I0]=r.e[I0]; e[I1]=r.e[I1]; e[I2]=r.e[I2]; e[I3]=r.e[I3]; return *this; } // Default copy operator has incorrect semantics
    };

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Implementation details. Do not make use of the contents of this namespace from outside the library. //
    /////////////////////////////////////////////////////////////////////////////////////////////////////////

    namespace detail
    {
        template<class A> struct element_storage { A elems; };

        // Support for user defined conversions
        template<class T, class U> using conv_t = typename std::enable_if<!std::is_same<T,U>::value, decltype(converter<T,U>()(std::declval<U>()))>::type;
        template<class T, class U> constexpr T convert(const U & u) { return converter<T,U>()(u); }

        // Helper which reduces scalar_accessor and swizzleN types to their intended targets
        template<class T> struct unpack { using type=T; };
        template<class T, class A, int I> struct unpack<scalar_accessor<T,A,I>> { using type=T; };
        template<class T, class A, int... I> struct unpack<swizzle2<T,A,I...>> { using type=vec<T,2>; };
        template<class T, class A, int... I> struct unpack<swizzle3<T,A,I...>> { using type=vec<T,3>; };
        template<class T, class A, int... I> struct unpack<swizzle4<T,A,I...>> { using type=vec<T,4>; };
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
        template<class T> struct any_compare<vec<T,2>,vec<T,2>> { using type=ord<T>; constexpr ord<T> operator() (const vec<T,2> & a, const vec<T,2> & b) const { return !(a[0]==b[0]) ? ord<T>{a[0],b[0]} : ord<T>{a[1],b[1]}; } };
        template<class T> struct any_compare<vec<T,3>,vec<T,3>> { using type=ord<T>; constexpr ord<T> operator() (const vec<T,3> & a, const vec<T,3> & b) const { return !(a[0]==b[0]) ? ord<T>{a[0],b[0]} : !(a[1]==b[1]) ? ord<T>{a[1],b[1]} : ord<T>{a[2],b[2]}; } };
        template<class T> struct any_compare<vec<T,4>,vec<T,4>> { using type=ord<T>; constexpr ord<T> operator() (const vec<T,4> & a, const vec<T,4> & b) const { return !(a[0]==b[0]) ? ord<T>{a[0],b[0]} : !(a[1]==b[1]) ? ord<T>{a[1],b[1]} : !(a[2]==b[2]) ? ord<T>{a[2],b[2]} : ord<T>{a[3],b[3]}; } };
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

        // Define vec/vec and vec/scalar patterns of up to three arguments to apply(...)
        template<class F, class Void, class... T> struct vec_apply {};
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

        // Define mat/mat and mat/scalar patterns of up to two arguments to apply(...)
        template<class F, class Void, class... T> struct mat_apply {};
        template<class F, int M, int N, class A         > struct mat_apply<F, scalars_t< >, mat<A,M,N>            > { using type = mat<ret_t<F,A  >,M,N>; enum {size=N}; template<int... J> static constexpr type impl(seq<J...>, F f, const mat<A,M,N> & a                      ) { return {vec_apply<F, void, vec<A,M>          >::impl(make_seq<M>{}, f, a[J]      )...}; } };
        template<class F, int M, int N, class A, class B> struct mat_apply<F, scalars_t< >, mat<A,M,N>, mat<B,M,N>> { using type = mat<ret_t<F,A,B>,M,N>; enum {size=N}; template<int... J> static constexpr type impl(seq<J...>, F f, const mat<A,M,N> & a, const mat<B,M,N> & b) { return {vec_apply<F, void, vec<A,M>, vec<B,M>>::impl(make_seq<M>{}, f, a[J], b[J])...}; } };
        template<class F, int M, int N, class A, class B> struct mat_apply<F, scalars_t<B>, mat<A,M,N>, B         > { using type = mat<ret_t<F,A,B>,M,N>; enum {size=N}; template<int... J> static constexpr type impl(seq<J...>, F f, const mat<A,M,N> & a, B                  b) { return {vec_apply<F, void, vec<A,M>, B       >::impl(make_seq<M>{}, f, a[J], b   )...}; } };
        template<class F, int M, int N, class A, class B> struct mat_apply<F, scalars_t<A>, A,          mat<B,M,N>> { using type = mat<ret_t<F,A,B>,M,N>; enum {size=N}; template<int... J> static constexpr type impl(seq<J...>, F f, A                  a, const mat<B,M,N> & b) { return {vec_apply<F, void, A,        vec<B,M>>::impl(make_seq<M>{}, f, a,    b[J])...}; } };

        // Define quat/quat and quat/scalar patterns of up to two arguments to apply(...)
        template<class F, class Void, class... T> struct quat_apply {};
        template<class F, class A         > struct quat_apply<F, scalars_t< >, quat<A>         > { using type = quat<ret_t<F,A  >>; enum {size=0}; static constexpr type impl(seq<>, F f, const quat<A> & a                   ) { return {f(a.x     ), f(a.y     ), f(a.z     ), f(a.w     )}; } };
        template<class F, class A, class B> struct quat_apply<F, scalars_t< >, quat<A>, quat<B>> { using type = quat<ret_t<F,A,B>>; enum {size=0}; static constexpr type impl(seq<>, F f, const quat<A> & a, const quat<B> & b) { return {f(a.x, b.x), f(a.y, b.y), f(a.z, b.z), f(a.w, b.w)}; } };
        template<class F, class A, class B> struct quat_apply<F, scalars_t<B>, quat<A>, B      > { using type = quat<ret_t<F,A,B>>; enum {size=0}; static constexpr type impl(seq<>, F f, const quat<A> & a, B               b) { return {f(a.x, b  ), f(a.y, b  ), f(a.z, b  ), f(a.w, b  )}; } };
        template<class F, class A, class B> struct quat_apply<F, scalars_t<A>, A,       quat<B>> { using type = quat<ret_t<F,A,B>>; enum {size=0}; static constexpr type impl(seq<>, F f, A               a, const quat<B> & b) { return {f(a,   b.x), f(a,   b.y), f(a,   b.z), f(a,   b.w)}; } };

        // Define scalar/scalar patterns for arbitrary numbers of arguments to apply(...)
        template<class F, class Void, class... T> struct scalar_apply {};
        template<class F, class... A> struct scalar_apply<F, scalars_t<A...>, A...> { using type = ret_t<F,A...>; enum {size=0}; static constexpr type impl(seq<>, F f, A... a) { return f(a...); } };

        // Aggregate all valid patterns for apply(...), apply_t<...> can be used for return-type SFINAE to control overload set
        template<class F, class... T> struct any_apply : vec_apply<F,void,T...>, mat_apply<F,void,T...>, quat_apply<F,void,T...>, scalar_apply<F,void,T...> {};

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

    //////////////////////////////////////////////////////////////
    // vec<T,M> specializations for 2, 3, and 4 element vectors //
    //////////////////////////////////////////////////////////////

    template<class T> struct vec<T,2>
    {
        union
        {
            detail::element_storage<T[2]> elems;
            scalar_accessor<T,T[2],0> x, r, s;
            scalar_accessor<T,T[2],1> y, g, t;
            swizzle2<T,T[2],0,1> xy, rg, st;
            swizzle2<T,T[2],1,0> yx, gr, ts;
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
            scalar_accessor<T,T[3],0> x, r, s;
            scalar_accessor<T,T[3],1> y, g, t;
            scalar_accessor<T,T[3],2> z, b, p;
            swizzle2<T,T[3],0,1> xy, rg, st;
            swizzle2<T,T[3],0,2> xz, rb, sp;
            swizzle2<T,T[3],1,0> yx, gr, ts;
            swizzle2<T,T[3],1,2> yz, gb, tp;
            swizzle2<T,T[3],2,0> zx, br, ps;
            swizzle2<T,T[3],2,1> zy, bg, pt;
            swizzle3<T,T[3],0,1,2> xyz, rgb, stp;
            swizzle3<T,T[3],0,2,1> xzy, rbg, spt;
            swizzle3<T,T[3],1,0,2> yxz, grb, tsp;
            swizzle3<T,T[3],1,2,0> yzx, gbr, tps;
            swizzle3<T,T[3],2,0,1> zxy, brg, pst;
            swizzle3<T,T[3],2,1,0> zyx, bgr, pts;
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
            scalar_accessor<T,T[4],0> x, r, s;
            scalar_accessor<T,T[4],1> y, g, t;
            scalar_accessor<T,T[4],2> z, b, p;
            scalar_accessor<T,T[4],3> w, a, q;
            swizzle2<T,T[4],0,1> xy, rg, st;
            swizzle2<T,T[4],0,2> xz, rb, sp;
            swizzle2<T,T[4],0,3> xw, ra, sq;
            swizzle2<T,T[4],1,0> yx, gr, ts;
            swizzle2<T,T[4],1,2> yz, gb, tp;
            swizzle2<T,T[4],1,3> yw, ga, tq;
            swizzle2<T,T[4],2,0> zx, br, ps;
            swizzle2<T,T[4],2,1> zy, bg, pt;
            swizzle2<T,T[4],2,3> zw, ba, pq;
            swizzle2<T,T[4],3,0> wx, ar, qs;
            swizzle2<T,T[4],3,1> wy, ag, qt;
            swizzle2<T,T[4],3,2> wz, ab, qp;
            swizzle3<T,T[4],0,1,2> xyz, rgb, stp;
            swizzle3<T,T[4],0,1,3> xyw, rga, stq;
            swizzle3<T,T[4],0,2,1> xzy, rbg, spt;
            swizzle3<T,T[4],0,2,3> xzw, rba, spq;
            swizzle3<T,T[4],0,3,1> xwy, rag, sqt;
            swizzle3<T,T[4],0,3,2> xwz, rab, sqp;
            swizzle3<T,T[4],1,0,2> yxz, grb, tsp;
            swizzle3<T,T[4],1,0,3> yxw, gra, tsq;
            swizzle3<T,T[4],1,2,0> yzx, gbr, tps;
            swizzle3<T,T[4],1,2,3> yzw, gba, tpq;
            swizzle3<T,T[4],1,3,0> ywx, gar, tqs;
            swizzle3<T,T[4],1,3,2> ywz, gab, tqp;
            swizzle3<T,T[4],2,0,1> zxy, brg, pst;
            swizzle3<T,T[4],2,0,3> zxw, bra, psq;
            swizzle3<T,T[4],2,1,0> zyx, bgr, pts;
            swizzle3<T,T[4],2,1,3> zyw, bga, ptq;
            swizzle3<T,T[4],2,3,0> zwx, bar, pqs;
            swizzle3<T,T[4],2,3,1> zwy, bag, pqt;
            swizzle3<T,T[4],3,0,1> wxy, arg, qst;
            swizzle3<T,T[4],3,0,2> wxz, arb, qsp;
            swizzle3<T,T[4],3,1,0> wyx, agr, qts;
            swizzle3<T,T[4],3,1,2> wyz, agb, qtp;
            swizzle3<T,T[4],3,2,0> wzx, abr, qps;
            swizzle3<T,T[4],3,2,1> wzy, abg, qpt;
            swizzle4<T,T[4],0,1,2,3> xyzw, rgba, stpq;
            swizzle4<T,T[4],0,1,3,2> xywz, rgab, stqp;
            swizzle4<T,T[4],0,2,1,3> xzyw, rbga, sptq;
            swizzle4<T,T[4],0,2,3,1> xzwy, rbag, spqt;
            swizzle4<T,T[4],0,3,1,2> xwyz, ragb, sqtp;
            swizzle4<T,T[4],0,3,2,1> xwzy, rabg, sqpt;
            swizzle4<T,T[4],1,0,2,3> yxzw, grba, tspq;
            swizzle4<T,T[4],1,0,3,2> yxwz, grab, tsqp;
            swizzle4<T,T[4],1,2,0,3> yzxw, gbra, tpsq;
            swizzle4<T,T[4],1,2,3,0> yzwx, gbar, tpqs;
            swizzle4<T,T[4],1,3,0,2> ywxz, garb, tqsp;
            swizzle4<T,T[4],1,3,2,0> ywzx, gabr, tqps;
            swizzle4<T,T[4],2,0,1,3> zxyw, brga, pstq;
            swizzle4<T,T[4],2,0,3,1> zxwy, brag, psqt;
            swizzle4<T,T[4],2,1,0,3> zyxw, bgra, ptsq;
            swizzle4<T,T[4],2,1,3,0> zywx, bgar, ptqs;
            swizzle4<T,T[4],2,3,0,1> zwxy, barg, pqst;
            swizzle4<T,T[4],2,3,1,0> zwyx, bagr, pqts;
            swizzle4<T,T[4],3,0,1,2> wxyz, argb, qstp;
            swizzle4<T,T[4],3,0,2,1> wxzy, arbg, qspt;
            swizzle4<T,T[4],3,1,0,2> wyxz, agrb, qtsp;
            swizzle4<T,T[4],3,1,2,0> wyzx, agbr, qtps;
            swizzle4<T,T[4],3,2,0,1> wzxy, abrg, qpst;
            swizzle4<T,T[4],3,2,1,0> wzyx, abgr, qpts;
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

    ////////////////////////////////////////////////////////////////
    // mat<T,M,N> specializations for 2, 3, and 4 column matrices //
    ////////////////////////////////////////////////////////////////
    
    template<class T, int M> struct mat<T,M,2>
    {
        using V=vec<T,M>;
        V                                       cols[2];
        constexpr                               mat()                                                    : cols{} {}
        constexpr                               mat(const V & x_, const V & y_)                          : cols{x_, y_} {}
        constexpr explicit                      mat(const T & s)                                         : cols{V(s), V(s)} {}
        template<class U> constexpr explicit    mat(const mat<U,M,2> & m)                                : cols{V(m[0]), V(m[1])} {}
        constexpr const V &                     operator[] (int j) const                                 { return cols[j]; }
        LINALG_CONSTEXPR14 V &                  operator[] (int j)                                       { return cols[j]; }
        constexpr vec<T,2>                      row(int i) const                                         { return {cols[0][i], cols[1][i]}; }

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
    template<class T, class F> constexpr T fold(const vec<T,2> & a, F f) { return f(a[0],a[1]); }
    template<class T, class F> constexpr T fold(const vec<T,3> & a, F f) { return f(f(a[0],a[1]),a[2]); }
    template<class T, class F> constexpr T fold(const vec<T,4> & a, F f) { return f(f(f(a[0],a[1]),a[2]),a[3]); }

    template<class T, int M, class F> constexpr T fold(const mat<T,M,2> & a, F f) { return f(fold(a[0],f),fold(a[1],f)); }
    template<class T, int M, class F> constexpr T fold(const mat<T,M,3> & a, F f) { return f(f(fold(a[0],f),fold(a[1],f)),fold(a[2],f)); }
    template<class T, int M, class F> constexpr T fold(const mat<T,M,4> & a, F f) { return f(f(f(fold(a[0],f),fold(a[1],f)),fold(a[2],f)),fold(a[3],f)); }   

    template<class T, class F> constexpr T fold(const quat<T> & a, F f) { return f(f(f(a.x,a.y),a.z),a.w); }

    // Type aliases for the result of calling apply(...) with various arguments, can be used with return type SFINAE to constrian overload sets
    template<class F, class... A> using apply_t = typename detail::any_apply<F,detail::unpack_t<A>...>::type;
    template<class F, class... A> using vec_apply_t = typename detail::vec_apply<F,void,detail::unpack_t<A>...>::type;
    template<class F, class... A> using mat_apply_t = typename detail::mat_apply<F,void,detail::unpack_t<A>...>::type;
    template<class F, class... A> using quat_apply_t = typename detail::quat_apply<F,void,detail::unpack_t<A>...>::type;

    // apply(f,...) applies the provided function in an elementwise fashion to its arguments, producing an object of the same dimensions
    template<class F, class... A> constexpr apply_t<F,A...> apply(F func, const A & ... args) { return detail::any_apply<F,detail::unpack_t<A>...>::impl(detail::make_seq<detail::any_apply<F,detail::unpack_t<A>...>::size>{}, func, args...); }

    // map(a,f) is equivalent to apply(f,a)
    template<class A, class F> constexpr apply_t<F,A> map(const A & a, F func) { return apply(func, a); }

    // zip(a,b,f) is equivalent to apply(f,a,b)
    template<class A, class B, class F> constexpr apply_t<F,A,B> zip(const A & a, const B & b, F func) { return apply(func, a, b); }

    ////////////////////////////////////
    // Vector operators and functions //
    ////////////////////////////////////

    // Component-wise unary operators: $vector
    template<class A> constexpr vec_apply_t<detail::op_pos, A> operator + (const A & a) { return apply(detail::op_pos{}, a); }
    template<class A> constexpr vec_apply_t<detail::op_neg, A> operator - (const A & a) { return apply(detail::op_neg{}, a); }
    template<class A> constexpr vec_apply_t<detail::op_cmp, A> operator ~ (const A & a) { return apply(detail::op_cmp{}, a); }
    template<class A> constexpr vec_apply_t<detail::op_not, A> operator ! (const A & a) { return apply(detail::op_not{}, a); }

    // Component-wise binary operators: vector $ vector; vector $ scalar; scalar $ vector
    template<class A, class B> constexpr vec_apply_t<detail::op_add, A, B> operator +  (const A & a, const B & b) { return apply(detail::op_add{}, a, b); }
    template<class A, class B> constexpr vec_apply_t<detail::op_sub, A, B> operator -  (const A & a, const B & b) { return apply(detail::op_sub{}, a, b); }
    template<class A, class B> constexpr vec_apply_t<detail::op_mul, A, B> operator *  (const A & a, const B & b) { return apply(detail::op_mul{}, a, b); }
    template<class A, class B> vec_apply_t<detail::op_div, A, B> operator /  (const A & a, const B & b) { return apply(detail::op_div{}, a, b); }
    template<class A, class B> constexpr vec_apply_t<detail::op_mod, A, B> operator %  (const A & a, const B & b) { return apply(detail::op_mod{}, a, b); }
    template<class A, class B> constexpr vec_apply_t<detail::op_un,  A, B> operator |  (const A & a, const B & b) { return apply(detail::op_un{},  a, b); }
    template<class A, class B> constexpr vec_apply_t<detail::op_xor, A, B> operator ^  (const A & a, const B & b) { return apply(detail::op_xor{}, a, b); }
    template<class A, class B> constexpr vec_apply_t<detail::op_int, A, B> operator &  (const A & a, const B & b) { return apply(detail::op_int{}, a, b); }
    template<class A, class B> constexpr vec_apply_t<detail::op_lsh, A, B> operator << (const A & a, const B & b) { return apply(detail::op_lsh{}, a, b); }
    template<class A, class B> constexpr vec_apply_t<detail::op_rsh, A, B> operator >> (const A & a, const B & b) { return apply(detail::op_rsh{}, a, b); }

    // Binary assignment operators: vector $= vector, vector $= scalar
    template<class T, int M, class B> constexpr vec<T,M> & operator +=  (vec<T,M> & a, const B & b) { return a = a + b; }
    template<class T, int M, class B> constexpr vec<T,M> & operator -=  (vec<T,M> & a, const B & b) { return a = a - b; }
    template<class T, int M, class B> constexpr vec<T,M> & operator *=  (vec<T,M> & a, const B & b) { return a = a * b; }
    template<class T, int M, class B> constexpr vec<T,M> & operator /=  (vec<T,M> & a, const B & b) { return a = a / b; }
    template<class T, int M, class B> constexpr vec<T,M> & operator %=  (vec<T,M> & a, const B & b) { return a = a % b; }
    template<class T, int M, class B> constexpr vec<T,M> & operator |=  (vec<T,M> & a, const B & b) { return a = a | b; }
    template<class T, int M, class B> constexpr vec<T,M> & operator ^=  (vec<T,M> & a, const B & b) { return a = a ^ b; }
    template<class T, int M, class B> constexpr vec<T,M> & operator &=  (vec<T,M> & a, const B & b) { return a = a & b; }
    template<class T, int M, class B> constexpr vec<T,M> & operator <<= (vec<T,M> & a, const B & b) { return a = a << b; }
    template<class T, int M, class B> constexpr vec<T,M> & operator >>= (vec<T,M> & a, const B & b) { return a = a >> b; }

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

    // Unary operators
    template<class T, int M, int N> constexpr mat<T,M,N> operator + (const mat<T,M,N> & a) { return apply(detail::op_pos{}, a); }
    template<class T, int M, int N> constexpr mat<T,M,N> operator - (const mat<T,M,N> & a) { return apply(detail::op_neg{}, a); }

    // Binary operators: matrix $ vector
    template<class T, int M> constexpr vec<T,M> operator * (const mat<T,M,2> & a, const vec<T,2> & b) { return a[0]*b[0] + a[1]*b[1]; }
    template<class T, int M> constexpr vec<T,M> operator * (const mat<T,M,3> & a, const vec<T,3> & b) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }
    template<class T, int M> constexpr vec<T,M> operator * (const mat<T,M,4> & a, const vec<T,4> & b) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3]; }

    // Binary operators: matrix $ matrix, matrix $ scalar, scalar $ matrix
    template<class T, int M, int N> constexpr mat<T,M,N> operator + (const mat<T,M,N> & a, const mat<T,M,N> & b) { return apply(detail::op_add{}, a, b); }
    template<class T, int M, int N> constexpr mat<T,M,N> operator - (const mat<T,M,N> & a, const mat<T,M,N> & b) { return apply(detail::op_sub{}, a, b); }
    template<class T, int M, int N> constexpr mat<T,M,2> operator * (const mat<T,M,N> & a, const mat<T,N,2> & b) { return {a*b[0], a*b[1]}; }
    template<class T, int M, int N> constexpr mat<T,M,3> operator * (const mat<T,M,N> & a, const mat<T,N,3> & b) { return {a*b[0], a*b[1], a*b[2]}; }
    template<class T, int M, int N> constexpr mat<T,M,4> operator * (const mat<T,M,N> & a, const mat<T,N,4> & b) { return {a*b[0], a*b[1], a*b[2], a*b[3]}; }
    template<class T, int M, int N> constexpr mat<T,M,N> operator * (const mat<T,M,N> & a, T b) { return apply(detail::op_mul{}, a, b); }
    template<class T, int M, int N> constexpr mat<T,M,N> operator * (T a, const mat<T,M,N> & b) { return apply(detail::op_mul{}, a, b); }
    template<class T, int M, int N> constexpr mat<T,M,N> operator / (const mat<T,M,N> & a, T b) { return apply(detail::op_div{}, a, b); }

    // Binary assignment operators
    template<class T, int M, int N> constexpr mat<T,M,N> & operator += (mat<T,M,N> & a, const mat<T,M,N> & b) { return a = a + b; }
    template<class T, int M, int N> constexpr mat<T,M,N> & operator -= (mat<T,M,N> & a, const mat<T,M,N> & b) { return a = a - b; }
    template<class T, int M, int N> constexpr mat<T,M,N> & operator *= (mat<T,M,N> & a, const mat<T,N,N> & b) { return a = a * b; }
    template<class T, int M, int N> constexpr mat<T,M,N> & operator *= (mat<T,M,N> & a, const T & b) { return a = a * b; }
    template<class T, int M, int N> constexpr mat<T,M,N> & operator /= (mat<T,M,N> & a, const T & b) { return a = a / b; }

    // Matrix algebra functions
    template<class T> constexpr vec<T,2>          diagonal    (const mat<T,2,2> & a) { return {a[0][0], a[1][1]}; }
    template<class T> constexpr vec<T,3>          diagonal    (const mat<T,3,3> & a) { return {a[0][0], a[1][1], a[2][2]}; }
    template<class T> constexpr vec<T,4>          diagonal    (const mat<T,4,4> & a) { return {a[0][0], a[1][1], a[2][2], a[3][3]}; }
    template<class T, int M> constexpr mat<T,M,2> outerprod   (const vec<T,M> & a, const vec<T,2> & b) { return {a*b[0], a*b[1]}; }
    template<class T, int M> constexpr mat<T,M,3> outerprod   (const vec<T,M> & a, const vec<T,3> & b) { return {a*b[0], a*b[1], a*b[2]}; }
    template<class T, int M> constexpr mat<T,M,4> outerprod   (const vec<T,M> & a, const vec<T,4> & b) { return {a*b[0], a*b[1], a*b[2], a*b[3]}; }
    template<class T, int M> constexpr mat<T,M,2> transpose   (const mat<T,2,M> & m) { return {m.row(0), m.row(1)}; }
    template<class T, int M> constexpr mat<T,M,3> transpose   (const mat<T,3,M> & m) { return {m.row(0), m.row(1), m.row(2)}; }
    template<class T, int M> constexpr mat<T,M,4> transpose   (const mat<T,4,M> & m) { return {m.row(0), m.row(1), m.row(2), m.row(3)}; }
    template<class T> constexpr mat<T,2,2>        adjugate    (const mat<T,2,2> & a) { return {{a[1][1], -a[0][1]}, {-a[1][0], a[0][0]}}; }
    template<class T> constexpr mat<T,3,3>        adjugate    (const mat<T,3,3> & a);
    template<class T> constexpr mat<T,4,4>        adjugate    (const mat<T,4,4> & a);
    template<class T> constexpr T                 determinant (const mat<T,2,2> & a) { return a[0][0]*a[1][1] - a[0][1]*a[1][0]; }
    template<class T> constexpr T                 determinant (const mat<T,3,3> & a) { return a[0][0]*(a[1][1]*a[2][2] - a[2][1]*a[1][2]) + a[0][1]*(a[1][2]*a[2][0] - a[2][2]*a[1][0]) + a[0][2]*(a[1][0]*a[2][1] - a[2][0]*a[1][1]); }
    template<class T> constexpr T                 determinant (const mat<T,4,4> & a);
    template<class T, int N> constexpr mat<T,N,N> inverse     (const mat<T,N,N> & a) { return adjugate(a)/determinant(a); }

    ////////////////////////////////////////
    // Quaternion operators and functions //
    ////////////////////////////////////////

    // Unary operators
    template<class T> constexpr quat<T> operator + (const quat<T> & a) { return apply(detail::op_pos{}, a); }
    template<class T> constexpr quat<T> operator - (const quat<T> & a) { return apply(detail::op_neg{}, a); }

    // Binary operators
    template<class T> constexpr quat<T> operator + (const quat<T> & a, const quat<T> & b) { return apply(detail::op_add{}, a, b); }
    template<class T> constexpr quat<T> operator - (const quat<T> & a, const quat<T> & b) { return apply(detail::op_sub{}, a, b); }
    template<class T> constexpr quat<T> operator * (const quat<T> & a, const quat<T> & b) { return {a.x*b.w+a.w*b.x+a.y*b.z-a.z*b.y, a.y*b.w+a.w*b.y+a.z*b.x-a.x*b.z, a.z*b.w+a.w*b.z+a.x*b.y-a.y*b.x, a.w*b.w-a.x*b.x-a.y*b.y-a.z*b.z}; }
    template<class T> constexpr quat<T> operator * (const quat<T> & a, T b) { return apply(detail::op_mul{}, a, b); }
    template<class T> constexpr quat<T> operator * (T a, const quat<T> & b) { return apply(detail::op_mul{}, a, b); }
    template<class T> constexpr quat<T> operator / (const quat<T> & a, T b) { return apply(detail::op_div{}, a, b); }

    // Binary assignment operators
    template<class T> constexpr quat<T> & operator += (quat<T> & a, const quat<T> & b) { return a = a + b; }
    template<class T> constexpr quat<T> & operator -= (quat<T> & a, const quat<T> & b) { return a = a - b; }
    template<class T> constexpr quat<T> & operator *= (quat<T> & a, const quat<T> & b) { return a = a * b; }
    template<class T> constexpr quat<T> & operator *= (quat<T> & a, const T & b) { return a = a * b; }
    template<class T> constexpr quat<T> & operator /= (quat<T> & a, const T & b) { return a = a / b; }

    // Quaternion algebra functions
    template<class T> quat<T>           qexp     (const quat<T> & q)                         { const auto v = q.xyz(); const auto vv = length(v); return std::exp(q.w) * vec<T,4>{v * (vv > 0 ? std::sin(vv)/vv : 0), std::cos(vv)}; }
    template<class T> quat<T>           qlog     (const quat<T> & q)                         { const auto v = q.xyz(); const auto vv = length(v), qq = length(q); return {v * (vv > 0 ? std::acos(q.w/qq)/vv : 0), std::log(qq)}; }
    template<class T> quat<T>           qpow     (const quat<T> & q, const T & p)            { const auto v = q.xyz(); const auto vv = length(v), qq = length(q), th = std::acos(q.w/qq); return std::pow(qq,p)*vec<T,4>{v * (vv > 0 ? std::sin(p*th)/vv : 0), std::cos(p*th)}; }
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
    template<class T> quat<T>              qnlerp(const quat<T> & a, const vec<T,4> & b, T t) { return nlerp(a, dot(a,b) < 0 ? -b : b, t); }
    template<class T> quat<T>              qslerp(const quat<T> & a, const vec<T,4> & b, T t) { return slerp(a, dot(a,b) < 0 ? -b : b, t); }

    //////////////////////////////////////////////////////////////////////////////////
    // Standard type aliases, can be brought in via using namespace linalg::aliases //
    //////////////////////////////////////////////////////////////////////////////////

    namespace aliases
    {
        using bool2=vec<bool,2>; using byte2=vec<uint8_t,2>; using short2=vec<int16_t,2>; using ushort2=vec<uint16_t,2>; 
        using bool3=vec<bool,3>; using byte3=vec<uint8_t,3>; using short3=vec<int16_t,3>; using ushort3=vec<uint16_t,3>; 
        using bool4=vec<bool,4>; using byte4=vec<uint8_t,4>; using short4=vec<int16_t,4>; using ushort4=vec<uint16_t,4>;
        using int2=vec<int,2>; using uint2=vec<unsigned,2>; using float2=vec<float,2>; using double2=vec<double,2>;
        using int3=vec<int,3>; using uint3=vec<unsigned,3>; using float3=vec<float,3>; using double3=vec<double,3>;
        using int4=vec<int,4>; using uint4=vec<unsigned,4>; using float4=vec<float,4>; using double4=vec<double,4>;
        using bool2x2=mat<bool,2,2>; using int2x2=mat<int,2,2>; using float2x2=mat<float,2,2>; using double2x2=mat<double,2,2>;
        using bool2x3=mat<bool,2,3>; using int2x3=mat<int,2,3>; using float2x3=mat<float,2,3>; using double2x3=mat<double,2,3>;
        using bool2x4=mat<bool,2,4>; using int2x4=mat<int,2,4>; using float2x4=mat<float,2,4>; using double2x4=mat<double,2,4>;
        using bool3x2=mat<bool,3,2>; using int3x2=mat<int,3,2>; using float3x2=mat<float,3,2>; using double3x2=mat<double,3,2>;
        using bool3x3=mat<bool,3,3>; using int3x3=mat<int,3,3>; using float3x3=mat<float,3,3>; using double3x3=mat<double,3,3>;
        using bool3x4=mat<bool,3,4>; using int3x4=mat<int,3,4>; using float3x4=mat<float,3,4>; using double3x4=mat<double,3,4>;
        using bool4x2=mat<bool,4,2>; using int4x2=mat<int,4,2>; using float4x2=mat<float,4,2>; using double4x2=mat<double,4,2>;
        using bool4x3=mat<bool,4,3>; using int4x3=mat<int,4,3>; using float4x3=mat<float,4,3>; using double4x3=mat<double,4,3>;
        using bool4x4=mat<bool,4,4>; using int4x4=mat<int,4,4>; using float4x4=mat<float,4,4>; using double4x4=mat<double,4,4>;
        using quatf=quat<float>; using quatd=quat<double>;
    }

    ///////////////////////////////////////////////////////////////////////////
    // Chopping block - Functions in this section may be reworked or removed //
    ///////////////////////////////////////////////////////////////////////////

    // Factory functions for 3D spatial transformations (will possibly be removed or changed in a future version)
    enum fwd_axis { neg_z, pos_z };                 // Should projection matrices be generated assuming forward is {0,0,-1} or {0,0,1}
    enum z_range { neg_one_to_one, zero_to_one };   // Should projection matrices map z into the range of [-1,1] or [0,1]?
    template<class T> vec<T,4>             rotation_quat     (const vec<T,3> & axis, T angle)         { return {axis*std::sin(angle/2), std::cos(angle/2)}; }
    template<class T> vec<T,4>             rotation_quat     (const mat<T,3,3> & m);
    template<class T> constexpr mat<T,4,4> translation_matrix(const vec<T,3> & translation)           { return {{1,0,0,0},{0,1,0,0},{0,0,1,0},{translation,1}}; }
    template<class T> constexpr mat<T,4,4> rotation_matrix   (const vec<T,4> & rotation)              { return {{qxdir(rotation),0}, {qydir(rotation),0}, {qzdir(rotation),0}, {0,0,0,1}}; }
    template<class T> constexpr mat<T,4,4> scaling_matrix    (const vec<T,3> & scaling)               { return {{scaling[0],0,0,0}, {0,scaling[1],0,0}, {0,0,scaling[2],0}, {0,0,0,1}}; }
    template<class T> constexpr mat<T,4,4> pose_matrix       (const vec<T,4> & q, const vec<T,3> & p) { return {{qxdir(q),0}, {qydir(q),0}, {qzdir(q),0}, {p,1}}; }
    template<class T> constexpr mat<T,4,4> frustum_matrix    (T x0, T x1, T y0, T y1, T n, T f, fwd_axis a = neg_z, z_range z = neg_one_to_one) { return {{2*n/(x1-x0),0,0,0}, {0,2*n/(y1-y0),0,0}, vec<T,4>{-(x0+x1)/(x1-x0), -(y0+y1)/(y1-y0), (z == zero_to_one ? f : f+n)/(f-n), 1} * (a == pos_z ? T(1) : T(-1)), {0,0,(z == zero_to_one ? -1 : -2)*n*f/(f-n),0}}; }
    template<class T> mat<T,4,4>           perspective_matrix(T fovy, T aspect, T n, T f, fwd_axis a = neg_z, z_range z = neg_one_to_one)       { T y = n*std::tan(fovy / 2), x = y*aspect; return frustum_matrix(-x, x, -y, y, n, f, a, z); }

    ////////////////////
    // Legacy support //
    ////////////////////

    // Support for quaternion algebra using 4D vectors, representing xi + yj + zk + w
    template<class T> vec<T,4>           qexp (const vec<T,4> & q)                     { const vec<T,3> v = q.xyz; const auto vv = length(v); return std::exp(q.w) * vec<T,4>{v * (vv > 0 ? std::sin(vv)/vv : 0), std::cos(vv)}; }
    template<class T> vec<T,4>           qlog (const vec<T,4> & q)                     { const vec<T,3> v = q.xyz; const auto vv = length(v), qq = length(q); return {v * (vv > 0 ? std::acos(q.w/qq)/vv : 0), std::log(qq)}; }
    template<class T> vec<T,4>           qpow (const vec<T,4> & q, const T & p)        { const vec<T,3> v = q.xyz; const auto vv = length(v), qq = length(q), th = std::acos(q.w/qq); return std::pow(qq,p)*vec<T,4>{v * (vv > 0 ? std::sin(p*th)/vv : 0), std::cos(p*th)}; }

    // Support for 3D spatial rotations using quaternions, via qmul(qmul(q, v), qconj(q))
    template<class T> constexpr vec<T,3>   qxdir (const vec<T,4> & q)                          { return {q.w*q.w+q.x*q.x-q.y*q.y-q.z*q.z, (q.x*q.y+q.z*q.w)*2, (q.z*q.x-q.y*q.w)*2}; }
    template<class T> constexpr vec<T,3>   qydir (const vec<T,4> & q)                          { return {(q.x*q.y-q.z*q.w)*2, q.w*q.w-q.x*q.x+q.y*q.y-q.z*q.z, (q.y*q.z+q.x*q.w)*2}; }
    template<class T> constexpr vec<T,3>   qzdir (const vec<T,4> & q)                          { return {(q.z*q.x+q.y*q.w)*2, (q.y*q.z-q.x*q.w)*2, q.w*q.w-q.x*q.x-q.y*q.y+q.z*q.z}; }
    template<class T> constexpr mat<T,3,3> qmat  (const vec<T,4> & q)                          { return {qxdir(q), qydir(q), qzdir(q)}; }
    template<class T> constexpr vec<T,3>   qrot  (const vec<T,4> & q, const vec<T,3> & v)      { return qxdir(q)*v.x + qydir(q)*v.y + qzdir(q)*v.z; }
    template<class T> T                    qangle(const vec<T,4> & q)                          { return std::atan2(length(vec<T,3>{q.xyz}), q.w)*2; }
    template<class T> vec<T,3>             qaxis (const vec<T,4> & q)                          { return normalize(vec<T,3>{q.xyz}); }
    template<class T> vec<T,4>             qnlerp(const vec<T,4> & a, const vec<T,4> & b, T t) { return nlerp(a, dot(a,b) < 0 ? -b : b, t); }
    template<class T> vec<T,4>             qslerp(const vec<T,4> & a, const vec<T,4> & b, T t) { return slerp(a, dot(a,b) < 0 ? -b : b, t); }

    // These functions exist to ease the difficulty of porting from older versions of linalg
    template<class T> constexpr vec<T,4> qconj(const vec<T,4> & q) { return vec<T,4>(conjugate(quat<T>(q))); }
    template<class T> constexpr vec<T,4> qinv(const vec<T,4> & q) { return vec<T,4>(inverse(quat<T>(q))); }
    template<class T> constexpr vec<T,4> qmul (const vec<T,4> & a, const vec<T,4> & b) { return vec<T,4>(quat<T>(a) * quat<T>(b)); }
    template<class T, class... R> constexpr vec<T,4> qmul(const vec<T,4> & a, R... r)  { return qmul(a, qmul(r...)); }
    template<class T, int M, int N, class B> constexpr auto mul(const mat<T,M,N> & a, const B & b) -> decltype(a*b) { return a*b; }
    template<class T, int M, int N, class... R> constexpr auto mul(const mat<T,M,N> & a, R... r) -> decltype(mul(a, mul(r...))) { return mul(a, mul(r...)); }
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

template<class T> linalg::vec<T,4> linalg::rotation_quat(const mat<T,3,3> & m)
{
    const vec<T,4> q {m[0][0]-m[1][1]-m[2][2], m[1][1]-m[0][0]-m[2][2], m[2][2]-m[0][0]-m[1][1], m[0][0]+m[1][1]+m[2][2]}, s[] {
        {1, m[0][1] + m[1][0], m[2][0] + m[0][2], m[1][2] - m[2][1]}, 
        {m[0][1] + m[1][0], 1, m[1][2] + m[2][1], m[2][0] - m[0][2]},
        {m[0][2] + m[2][0], m[1][2] + m[2][1], 1, m[0][1] - m[1][0]},
        {m[1][2] - m[2][1], m[2][0] - m[0][2], m[0][1] - m[1][0], 1}};
    return copysign(normalize(sqrt(max(T(0), T(1)+q))), s[argmax(q)]);
}

#endif
