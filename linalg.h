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
#include <utility>      // For std::integer_sequence
#include <array>        // For std::array

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

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Implementation details. Do not make use of the contents of this namespace from outside the library. //
    /////////////////////////////////////////////////////////////////////////////////////////////////////////

    namespace detail
    {
        // Type with an implicit conversion to the multiplicative identity of any given algebraic object
        struct identity_t { constexpr identity_t(int) {}; };

        // Type returned by the compare(...) function which supports all six comparison operators against 0
        template<class T> struct ord { T a,b; };
        template<class T> constexpr bool operator == (const ord<T> & o, std::nullptr_t) { return o.a == o.b; }
        template<class T> constexpr bool operator != (const ord<T> & o, std::nullptr_t) { return !(o.a == o.b); }
        template<class T> constexpr bool operator < (const ord<T> & o, std::nullptr_t) { return o.a < o.b; }
        template<class T> constexpr bool operator > (const ord<T> & o, std::nullptr_t) { return o.b < o.a; }
        template<class T> constexpr bool operator <= (const ord<T> & o, std::nullptr_t) { return !(o.b < o.a); }
        template<class T> constexpr bool operator >= (const ord<T> & o, std::nullptr_t) { return !(o.a < o.b); }

        // Value type wrapper around a C array, for use in single-expression constexpr functions
        template<class T, int N> struct arr { T e[N]; };

        // SFINAE helpers to determine result of function application
        template<class F, class... T> using ret_t = decltype(std::declval<F>()(std::declval<T>()...));

        // SFINAE helper which is defined if all provided types are scalars
        struct empty {};
        template<class... T> struct scalars;
        template<> struct scalars<> { using type=void; };
        template<class T, class... U> struct scalars<T,U...> : std::conditional_t<std::is_arithmetic_v<T>, scalars<U...>, empty> {};
        template<class... T> using scalars_t = typename scalars<T...>::type;

        // Define vec/vec and vec/scalar patterns for apply(...), vec_apply_t can be used for return-type SFINAE to control overload set
        template<class V, class F, class... T> struct vec_apply {};
        template<class F, int M, class A                  > struct vec_apply<scalars_t<   >, F, vec<A,M>                    > { using type=vec<ret_t<F,A    >,M>; enum {size=M}; template<int... I> static constexpr type impl(std::integer_sequence<int,I...>, F f, const vec<A,M> & a                                        ) { return {f(a[I]            )...}; } };
        template<class F, int M, class A, class B         > struct vec_apply<scalars_t<B  >, F, vec<A,M>, B                 > { using type=vec<ret_t<F,A,B  >,M>; enum {size=M}; template<int... I> static constexpr type impl(std::integer_sequence<int,I...>, F f, const vec<A,M> & a, B                b                    ) { return {f(a[I], b         )...}; } };
        template<class F, int M, class A, class B         > struct vec_apply<scalars_t<A  >, F, A,        vec<B,M>          > { using type=vec<ret_t<F,A,B  >,M>; enum {size=M}; template<int... I> static constexpr type impl(std::integer_sequence<int,I...>, F f, A                a, const vec<B,M> & b                    ) { return {f(a,    b[I]      )...}; } };
        template<class F, int M, class A, class B         > struct vec_apply<scalars_t<   >, F, vec<A,M>, vec<B,M>          > { using type=vec<ret_t<F,A,B  >,M>; enum {size=M}; template<int... I> static constexpr type impl(std::integer_sequence<int,I...>, F f, const vec<A,M> & a, const vec<B,M> & b                    ) { return {f(a[I], b[I]      )...}; } };
        template<class F, int M, class A, class B, class C> struct vec_apply<scalars_t<B,C>, F, vec<A,M>, B,        C       > { using type=vec<ret_t<F,A,B,C>,M>; enum {size=M}; template<int... I> static constexpr type impl(std::integer_sequence<int,I...>, F f, const vec<A,M> & a, B                b, C                c) { return {f(a[I], b,    c   )...}; } };
        template<class F, int M, class A, class B, class C> struct vec_apply<scalars_t<A,C>, F, A,        vec<B,M>, C       > { using type=vec<ret_t<F,A,B,C>,M>; enum {size=M}; template<int... I> static constexpr type impl(std::integer_sequence<int,I...>, F f, A                a, const vec<B,M> & b, C                c) { return {f(a,    b[I], c   )...}; } };
        template<class F, int M, class A, class B, class C> struct vec_apply<scalars_t<A,B>, F, A,        B,        vec<C,M>> { using type=vec<ret_t<F,A,B,C>,M>; enum {size=M}; template<int... I> static constexpr type impl(std::integer_sequence<int,I...>, F f, A                a, B                b, const vec<C,M> & c) { return {f(a,    b,    c[I])...}; } };
        template<class F, int M, class A, class B, class C> struct vec_apply<scalars_t<C  >, F, vec<A,M>, vec<B,M>, C       > { using type=vec<ret_t<F,A,B,C>,M>; enum {size=M}; template<int... I> static constexpr type impl(std::integer_sequence<int,I...>, F f, const vec<A,M> & a, const vec<B,M> & b, C                c) { return {f(a[I], b[I], c   )...}; } };
        template<class F, int M, class A, class B, class C> struct vec_apply<scalars_t<B  >, F, vec<A,M>, B,        vec<C,M>> { using type=vec<ret_t<F,A,B,C>,M>; enum {size=M}; template<int... I> static constexpr type impl(std::integer_sequence<int,I...>, F f, const vec<A,M> & a, B                b, const vec<C,M> & c) { return {f(a[I], b,    c[I])...}; } };
        template<class F, int M, class A, class B, class C> struct vec_apply<scalars_t<A  >, F, A,        vec<B,M>, vec<C,M>> { using type=vec<ret_t<F,A,B,C>,M>; enum {size=M}; template<int... I> static constexpr type impl(std::integer_sequence<int,I...>, F f, A                a, const vec<B,M> & b, const vec<C,M> & c) { return {f(a,    b[I], c[I])...}; } };
        template<class F, int M, class A, class B, class C> struct vec_apply<scalars_t<   >, F, vec<A,M>, vec<B,M>, vec<C,M>> { using type=vec<ret_t<F,A,B,C>,M>; enum {size=M}; template<int... I> static constexpr type impl(std::integer_sequence<int,I...>, F f, const vec<A,M> & a, const vec<B,M> & b, const vec<C,M> & c) { return {f(a[I], b[I], c[I])...}; } };

        // Define mat/mat and mat/scalar patterns for apply(...), mat_apply_t<...> can be used for return-type SFINAE to control overload set
        template<class F, class... T> struct mat_apply {};
        template<class F, class T, int M, int N> struct mat_apply<F, mat<T,M,N>>             { using type = mat<ret_t<F,T>,M,N>;   enum {size=N}; template<int... J> static constexpr type impl(std::integer_sequence<int,J...>, F f, const mat<T,M,N> & a)                       { return {vec_apply<void,F,vec<T,M>         >::impl(std::make_integer_sequence<int,M>{}, f, a[J])...}; } };
        template<class F, class T, int M, int N> struct mat_apply<F, mat<T,M,N>, T>          { using type = mat<ret_t<F,T,T>,M,N>; enum {size=N}; template<int... J> static constexpr type impl(std::integer_sequence<int,J...>, F f, const mat<T,M,N> & a, T b)                  { return {vec_apply<void,F,vec<T,M>,T       >::impl(std::make_integer_sequence<int,M>{}, f, a[J], b)...}; }};
        template<class F, class T, int M, int N> struct mat_apply<F, T, mat<T,M,N>>          { using type = mat<ret_t<F,T,T>,M,N>; enum {size=N}; template<int... J> static constexpr type impl(std::integer_sequence<int,J...>, F f, T a, const mat<T,M,N> & b)                  { return {vec_apply<void,F,T,vec<T,M>       >::impl(std::make_integer_sequence<int,M>{}, f, a, b[J])...}; }};
        template<class F, class T, int M, int N> struct mat_apply<F, mat<T,M,N>, mat<T,M,N>> { using type = mat<ret_t<F,T,T>,M,N>; enum {size=N}; template<int... J> static constexpr type impl(std::integer_sequence<int,J...>, F f, const mat<T,M,N> & a, const mat<T,M,N> & b) { return {vec_apply<void,F,vec<T,M>,vec<T,M>>::impl(std::make_integer_sequence<int,M>{}, f, a[J], b[J])...}; }};

        // Define quat/quat and quat/scalar patterns for apply(...), quat_apply_t<...> can be used for return-type SFINAE to control overload set
        template<class F, class... T> struct quat_apply {};
        template<class F, class T> struct quat_apply<F, quat<T>>          { using type = quat<ret_t<F,T>>;   enum {size=0}; static constexpr type impl(std::integer_sequence<int>, F f, const quat<T> & a)                    { return {f(a.x), f(a.y), f(a.z), f(a.w)}; }                 };
        template<class F, class T> struct quat_apply<F, quat<T>, T>       { using type = quat<ret_t<F,T,T>>; enum {size=0}; static constexpr type impl(std::integer_sequence<int>, F f, const quat<T> & a, T b)               { return {f(a.x,b), f(a.y,b), f(a.z,b), f(a.w,b)}; }         };
        template<class F, class T> struct quat_apply<F, T, quat<T>>       { using type = quat<ret_t<F,T,T>>; enum {size=0}; static constexpr type impl(std::integer_sequence<int>, F f, T a, const quat<T> & b)               { return {f(a,b.x), f(a,b.y), f(a,b.z), f(a,b.w)}; }         };
        template<class F, class T> struct quat_apply<F, quat<T>, quat<T>> { using type = quat<ret_t<F,T,T>>; enum {size=0}; static constexpr type impl(std::integer_sequence<int>, F f, const quat<T> & a, const quat<T> & b) { return {f(a.x,b.x), f(a.y,b.y), f(a.z,b.z), f(a.w,b.w)}; } };

        // Aggregate all valid patterns for apply(...), apply_t<...> can be used for return-type SFINAE to control overload set
        template<class F, class... T> struct any_apply : vec_apply<void,F,T...>, mat_apply<F,T...>, quat_apply<F,T...> {};

        // Function objects for selecting between alternatives
        struct min { template<class T> constexpr T operator() (T a, T b) const { return a < b ? a : b; } };
        struct max { template<class T> constexpr T operator() (T a, T b) const { return a < b ? b : a; } };
        struct clamp { template<class T> constexpr T operator() (T a, T b, T c) const { return a < b ? b : a < c ? a : c; } };
        struct select { template<class T> constexpr T operator() (bool a, T b, T c) const { return a ? b : c; } };

        // Function objects for applying operators
        struct op_pos { template<class T> constexpr auto operator() (T a) const { return +a; } };
        struct op_neg { template<class T> constexpr auto operator() (T a) const { return -a; } };
        struct op_not { template<class T> constexpr auto operator() (T a) const { return !a; } };
        struct op_cmp { template<class T> constexpr auto operator() (T a) const { return ~a; } };
        struct op_mul { template<class T> constexpr auto operator() (T a, T b) const { return a * b; } };
        struct op_div { template<class T> constexpr auto operator() (T a, T b) const { return a / b; } };
        struct op_mod { template<class T> constexpr auto operator() (T a, T b) const { return a % b; } };
        struct op_add { template<class T> constexpr auto operator() (T a, T b) const { return a + b; } };
        struct op_sub { template<class T> constexpr auto operator() (T a, T b) const { return a - b; } };
        struct op_lsh { template<class T> constexpr auto operator() (T a, T b) const { return a << b; } };
        struct op_rsh { template<class T> constexpr auto operator() (T a, T b) const { return a >> b; } };
        struct op_lt  { template<class T> constexpr auto operator() (T a, T b) const { return a < b; } };
        struct op_gt  { template<class T> constexpr auto operator() (T a, T b) const { return a > b; } };
        struct op_le  { template<class T> constexpr auto operator() (T a, T b) const { return a <= b; } };
        struct op_ge  { template<class T> constexpr auto operator() (T a, T b) const { return a >= b; } };
        struct op_eq  { template<class T> constexpr auto operator() (T a, T b) const { return a == b; } };
        struct op_ne  { template<class T> constexpr auto operator() (T a, T b) const { return a != b; } };
        struct op_int { template<class T> constexpr auto operator() (T a, T b) const { return a & b; } };        
        struct op_xor { template<class T> constexpr auto operator() (T a, T b) const { return a ^ b; } };
        struct op_un  { template<class T> constexpr auto operator() (T a, T b) const { return a | b; } };
        struct op_and { template<class T> constexpr auto operator() (T a, T b) const { return a && b; } };
        struct op_or  { template<class T> constexpr auto operator() (T a, T b) const { return a || b; } };

        // Function objects for applying standard library math functions
        struct std_abs      { template<class T> auto operator() (T a) const { return std::abs  (a); } };
        struct std_floor    { template<class T> auto operator() (T a) const { return std::floor(a); } };
        struct std_ceil     { template<class T> auto operator() (T a) const { return std::ceil (a); } };
        struct std_exp      { template<class T> auto operator() (T a) const { return std::exp  (a); } };
        struct std_log      { template<class T> auto operator() (T a) const { return std::log  (a); } };
        struct std_log10    { template<class T> auto operator() (T a) const { return std::log10(a); } };
        struct std_sqrt     { template<class T> auto operator() (T a) const { return std::sqrt (a); } };
        struct std_sin      { template<class T> auto operator() (T a) const { return std::sin  (a); } };
        struct std_cos      { template<class T> auto operator() (T a) const { return std::cos  (a); } };
        struct std_tan      { template<class T> auto operator() (T a) const { return std::tan  (a); } };
        struct std_asin     { template<class T> auto operator() (T a) const { return std::asin (a); } };
        struct std_acos     { template<class T> auto operator() (T a) const { return std::acos (a); } };
        struct std_atan     { template<class T> auto operator() (T a) const { return std::atan (a); } };
        struct std_sinh     { template<class T> auto operator() (T a) const { return std::sinh (a); } };
        struct std_cosh     { template<class T> auto operator() (T a) const { return std::cosh (a); } };
        struct std_tanh     { template<class T> auto operator() (T a) const { return std::tanh (a); } };
        struct std_round    { template<class T> auto operator() (T a) const { return std::round(a); } };
        struct std_fmod     { template<class T> auto operator() (T a, T b) const { return std::fmod    (a, b); } };
        struct std_pow      { template<class T> auto operator() (T a, T b) const { return std::pow     (a, b); } };
        struct std_atan2    { template<class T> auto operator() (T a, T b) const { return std::atan2   (a, b); } };
        struct std_copysign { template<class T> auto operator() (T a, T b) const { return std::copysign(a, b); } };
    }

    //////////////////////////////////////////////////////////////
    // vec<T,M> specializations for 2, 3, and 4 element vectors //
    //////////////////////////////////////////////////////////////

    template<class T> struct vec<T,2>
    {
        T                           x,y;
        constexpr                   vec()                               : x(), y() {}
        constexpr                   vec(const T & x_, const T & y_)     : x(x_), y(y_) {}
        constexpr                   vec(const std::array<T,2> & a)      : vec(a[0], a[1]) {}
        constexpr explicit          vec(const T & s)                    : vec(s, s) {}
        constexpr explicit          vec(const T * p)                    : vec(p[0], p[1]) {}
        template<class U>
        constexpr explicit          vec(const vec<U,2> & v)             : vec(static_cast<T>(v.x), static_cast<T>(v.y)) {}
        constexpr const T &         operator[] (int i) const            { return this->*(detail::arr<T vec::*,2>{&vec::x, &vec::y}.e[i]); }
        LINALG_CONSTEXPR14 T &      operator[] (int i)                  { return this->*(detail::arr<T vec::*,2>{&vec::x, &vec::y}.e[i]); }
    };
    template<class T> struct vec<T,3>
    {
        T                           x,y,z;
        constexpr                   vec()                               : x(), y(), z() {}
        constexpr                   vec(const T & x_, const T & y_, 
                                        const T & z_)                   : x(x_), y(y_), z(z_) {}
        constexpr                   vec(const vec<T,2> & xy,
                                        const T & z_)                   : vec(xy.x, xy.y, z_) {}
        constexpr                   vec(const std::array<T,3> & a)      : vec(a[0], a[1], a[2]) {}
        constexpr explicit          vec(const T & s)                    : vec(s, s, s) {}
        constexpr explicit          vec(const T * p)                    : vec(p[0], p[1], p[2]) {}
        template<class U>
        constexpr explicit          vec(const vec<U,3> & v)             : vec(static_cast<T>(v.x), static_cast<T>(v.y), static_cast<T>(v.z)) {}
        constexpr const T &         operator[] (int i) const            { return this->*(detail::arr<T vec::*,3>{&vec::x, &vec::y, &vec::z}.e[i]); }
        LINALG_CONSTEXPR14 T &      operator[] (int i)                  { return this->*(detail::arr<T vec::*,3>{&vec::x, &vec::y, &vec::z}.e[i]); }
        constexpr vec<T,2>          xy() const                          { return {x,y}; }
    };
    template<class T> struct vec<T,4>
    {
        T                           x,y,z,w;
        constexpr                   vec()                               : x(), y(), z(), w() {}
        constexpr                   vec(const T & x_, const T & y_,
                                        const T & z_, const T & w_)     : x(x_), y(y_), z(z_), w(w_) {}
        constexpr                   vec(const vec<T,2> & xy, 
                                        const T & z_, const T & w_)     : vec(xy.x, xy.y, z_, w_) {}
        constexpr                   vec(const vec<T,3> & xyz,
                                        const T & w_)                   : vec(xyz.x, xyz.y, xyz.z, w_) {}
        constexpr                   vec(const std::array<T,4> & a)      : vec(a[0], a[1], a[2], a[3]) {}
        constexpr explicit          vec(const T & s)                    : vec(s, s, s, s) {}
        constexpr explicit          vec(const quat<T> & q)              : vec(q.x, q.y, q.z, q.w) {}
        constexpr explicit          vec(const T * p)                    : vec(p[0], p[1], p[2], p[3]) {}
        template<class U> 
        constexpr explicit          vec(const vec<U,4> & v)             : vec(static_cast<T>(v.x), static_cast<T>(v.y), static_cast<T>(v.z), static_cast<T>(v.w)) {}
        constexpr const T &         operator[] (int i) const            { return this->*(detail::arr<T vec::*,4>{&vec::x, &vec::y, &vec::z, &vec::w}.e[i]); }
        LINALG_CONSTEXPR14 T &      operator[] (int i)                  { return this->*(detail::arr<T vec::*,4>{&vec::x, &vec::y, &vec::z, &vec::w}.e[i]); }
        constexpr vec<T,3>          xyz() const                         { return {x,y,z}; }
        constexpr vec<T,2>          xy() const                          { return {x,y}; }
    };

    ////////////////////////////////////////////////////////////////
    // mat<T,M,N> specializations for 2, 3, and 4 column matrices //
    ////////////////////////////////////////////////////////////////
    
    template<class T, int M> struct mat<T,M,2>
    {
        using V=vec<T,M>;
        V                           cols[2];
        constexpr                   mat()                               : cols{} {}
        constexpr                   mat(detail::identity_t)             : cols{{1,0},{0,1}} {}
        constexpr                   mat(const V & x_, const V & y_)     : cols{x_, y_} {}
        constexpr explicit          mat(const T & s)                    : cols{V(s), V(s)} {}
        constexpr explicit          mat(const T * p)                    : cols{V(p+M*0), V(p+M*1)} {}
        template<class U> 
        constexpr explicit          mat(const mat<U,M,2> & m)           : cols{V(m[0]), V(m[1])} {}
        constexpr const V &         operator[] (int j) const            { return cols[j]; }
        LINALG_CONSTEXPR14 V &      operator[] (int j)                  { return cols[j]; }
        constexpr vec<T,2>          row(int i) const                    { return {cols[0][i], cols[1][i]}; }
    };
    template<class T, int M> struct mat<T,M,3>
    {
        using V=vec<T,M>;
        V                           cols[3];
        constexpr                   mat()                               : cols{} {}
        constexpr                   mat(detail::identity_t)             : cols{{1,0,0},{0,1,0},{0,0,1}} {}
        constexpr                   mat(const V & x_, const V & y_, 
                                        const V & z_)                   : cols{x_, y_, z_} {}
        constexpr explicit          mat(const T & s)                    : cols{V(s), V(s), V(s)} {}
        constexpr explicit          mat(const T * p)                    : cols{V(p+M*0), V(p+M*1), V(p+M*2)} {}
        template<class U> 
        constexpr explicit          mat(const mat<U,M,3> & m)           : cols{V(m[0]), V(m[1]), V(m[2])} {}
        constexpr const V &         operator[] (int j) const            { return cols[j]; }
        LINALG_CONSTEXPR14 V &      operator[] (int j)                  { return cols[j]; }
        constexpr vec<T,3>          row(int i) const                    { return {cols[0][i], cols[1][i], cols[2][i]}; }
    };
    template<class T, int M> struct mat<T,M,4>
    {
        using V=vec<T,M>;
        V                           cols[4];
        constexpr                   mat()                               : cols{} {}
        constexpr                   mat(detail::identity_t)             : cols{{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}} {}
        constexpr                   mat(const V & x_, const V & y_,
                                        const V & z_, const V & w_)     : cols{x_, y_, z_, w_} {}
        constexpr explicit          mat(const T & s)                    : cols{V(s), V(s), V(s), V(s)} {}
        constexpr explicit          mat(const T * p)                    : cols{V(p+M*0), V(p+M*1), V(p+M*2), V(p+M*3)} {}
        template<class U> 
        constexpr explicit          mat(const mat<U,M,4> & m)           : cols{V(m[0]), V(m[1]), V(m[2]), V(m[3])} {}
        constexpr const V &         operator[] (int j) const            { return cols[j]; }
        LINALG_CONSTEXPR14 V &      operator[] (int j)                  { return cols[j]; }
        constexpr vec<T,4>          row(int i) const                    { return {cols[0][i], cols[1][i], cols[2][i], cols[3][i]}; }
    };

    ////////////////////////////
    // quat<T> implementation //
    ////////////////////////////

    template<class T> struct quat
    {
        T x,y,z,w;
        constexpr                   quat()                              : x(), y(), z(), w() {}
        constexpr                   quat(detail::identity_t)            : quat(0,0,0,1) {}
        constexpr                   quat(const T & x_, const T & y_,
                                         const T & z_, const T & w_)    : x(x_), y(y_), z(z_), w(w_) {}
        constexpr                   quat(const vec<T,3> & xyz,
                                         const T & w_)                  : quat(xyz.x, xyz.y, xyz.z, w_) {}
        constexpr explicit          quat(const std::array<T,4> & a)     : quat(a[0], a[1], a[2], a[3]) {}
        constexpr explicit          quat(const vec<T,4> & xyzw)         : quat(xyzw.x, xyzw.y, xyzw.z, xyzw.w) {}
        constexpr explicit          quat(const T * p)                   : quat(p[0], p[1], p[2], p[3]) {}
        template<class U> 
        constexpr explicit          quat(const quat<U> & q)             : quat(static_cast<T>(q.x), static_cast<T>(q.y), static_cast<T>(q.z), static_cast<T>(q.w)) {}
        constexpr vec<T,3>          xyz() const                         { return {x,y,z}; }
    };

    ///////////////
    // Constants //
    ///////////////

    // Converts implicity to the multiplicative identity of any given algebraic object 
    constexpr detail::identity_t identity {1};

    //////////////////////////
    // Relational operators //
    //////////////////////////

    template<class A> constexpr auto operator == (const A & a, const A & b) -> decltype(compare(a,b) == 0) { return compare(a,b) == 0; }
    template<class A> constexpr auto operator != (const A & a, const A & b) -> decltype(compare(a,b) != 0) { return compare(a,b) != 0; }
    template<class A> constexpr auto operator <  (const A & a, const A & b) -> decltype(compare(a,b) <  0) { return compare(a,b) <  0; }
    template<class A> constexpr auto operator >  (const A & a, const A & b) -> decltype(compare(a,b) >  0) { return compare(a,b) >  0; }
    template<class A> constexpr auto operator <= (const A & a, const A & b) -> decltype(compare(a,b) <= 0) { return compare(a,b) <= 0; }
    template<class A> constexpr auto operator >= (const A & a, const A & b) -> decltype(compare(a,b) >= 0) { return compare(a,b) >= 0; }

    template<class T> constexpr detail::ord<T> compare(const vec<T,2> & a, const vec<T,2> & b) { return !(a.x==b.x) ? detail::ord<T>{a.x,b.x} : detail::ord<T>{a.y,b.y}; }
    template<class T> constexpr detail::ord<T> compare(const vec<T,3> & a, const vec<T,3> & b) { return !(a.x==b.x) ? detail::ord<T>{a.x,b.x} : !(a.y==b.y) ? detail::ord<T>{a.y,b.y} : detail::ord<T>{a.z,b.z}; }
    template<class T> constexpr detail::ord<T> compare(const vec<T,4> & a, const vec<T,4> & b) { return !(a.x==b.x) ? detail::ord<T>{a.x,b.x} : !(a.y==b.y) ? detail::ord<T>{a.y,b.y} : !(a.z==b.z) ? detail::ord<T>{a.z,b.z} : detail::ord<T>{a.w,b.w}; }

    template<class T, int M> constexpr detail::ord<T> compare(const mat<T,M,2> & a, const mat<T,M,2> & b) { return a[0]!=b[0] ? compare(a[0],b[0]) : compare(a[1],b[1]); }
    template<class T, int M> constexpr detail::ord<T> compare(const mat<T,M,3> & a, const mat<T,M,3> & b) { return a[0]!=b[0] ? compare(a[0],b[0]) : a[1]!=b[1] ? compare(a[1],b[1]) : compare(a[2],b[2]); }
    template<class T, int M> constexpr detail::ord<T> compare(const mat<T,M,4> & a, const mat<T,M,4> & b) { return a[0]!=b[0] ? compare(a[0],b[0]) : a[1]!=b[1] ? compare(a[1],b[1]) : a[2]!=b[2] ? compare(a[2],b[2]) : compare(a[3],b[3]); }

    template<class T> constexpr detail::ord<T> compare(const quat<T> & a, const quat<T> & b) { return !(a.x==b.x) ? detail::ord<T>{a.x,b.x} : !(a.y==b.y) ? detail::ord<T>{a.y,b.y} : !(a.z==b.z) ? detail::ord<T>{a.z,b.z} : detail::ord<T>{a.w,b.w}; }

    ////////////////////////////
    // Higher-order functions //
    ////////////////////////////

    // Produce a scalar by applying f(T,T) -> T to adjacent pairs of elements from vector/matrix a in left-to-right order (matching the associativity of arithmetic and logical operators)
    template<class T, class F> constexpr T fold(const vec<T,2> & a, F f) { return f(a.x,a.y); }
    template<class T, class F> constexpr T fold(const vec<T,3> & a, F f) { return f(f(a.x,a.y),a.z); }
    template<class T, class F> constexpr T fold(const vec<T,4> & a, F f) { return f(f(f(a.x,a.y),a.z),a.w); }

    template<class T, int M, class F> constexpr T fold(const mat<T,M,2> & a, F f) { return f(fold(a[0],f),fold(a[1],f)); }
    template<class T, int M, class F> constexpr T fold(const mat<T,M,3> & a, F f) { return f(f(fold(a[0],f),fold(a[1],f)),fold(a[2],f)); }
    template<class T, int M, class F> constexpr T fold(const mat<T,M,4> & a, F f) { return f(f(f(fold(a[0],f),fold(a[1],f)),fold(a[2],f)),fold(a[3],f)); }   

    template<class T, class F> constexpr T fold(const quat<T> & a, F f) { return f(f(f(a.x,a.y),a.z),a.w); }

    // apply(f,...) applies the provided function in an elementwise fashion to its arguments, producing an object of the same dimensions
    template<class F, class... A> constexpr auto apply(F func, const A & ... args) { return detail::any_apply<F,A...>::impl(std::make_integer_sequence<int,detail::any_apply<F,A...>::size>{}, func, args...); }

    // map(a,f) is equivalent to apply(f,a)
    template<class A, class F> constexpr auto map(const A & a, F func) { return apply(func, a); }

    // zip(a,b,f) is equivalent to apply(f,a,b)
    template<class A, class B, class F> constexpr auto zip(const A & a, const B & b, F func) { return apply(func, a, b); }

    // Type aliases for the result of calling apply(...) with various arguments, can be used with return type SFINAE to constrian overload sets
    template<class F, class... A> using apply_t = typename detail::any_apply<F,A...>::type;
    template<class F, class... A> using vec_apply_t = typename detail::vec_apply<void,F,A...>::type;
    template<class F, class... A> using mat_apply_t = typename detail::mat_apply<F,A...>::type;
    template<class F, class... A> using quat_apply_t = typename detail::quat_apply<F,A...>::type;

    ////////////////////////////////////
    // Vector operators and functions //
    ////////////////////////////////////

    // Component-wise unary operators: $vector
    template<class T, int M> constexpr auto operator + (const vec<T,M> & a) { return apply(detail::op_pos{}, a); }
    template<class T, int M> constexpr auto operator - (const vec<T,M> & a) { return apply(detail::op_neg{}, a); }
    template<class T, int M> constexpr auto operator ~ (const vec<T,M> & a) { return apply(detail::op_cmp{}, a); }
    template<class T, int M> constexpr auto operator ! (const vec<T,M> & a) { return apply(detail::op_not{}, a); }

    // Component-wise binary operators: vector $ vector; vector $ scalar; scalar $ vector
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
    template<class A, class B, class C> constexpr vec_apply_t<detail::clamp,  A, B, C> clamp (const A & a, const B & b, const C & c) { return apply(detail::clamp{},  a, b, c); }
    template<class A, class B, class C> constexpr vec_apply_t<detail::select, A, B, C> select(const A & a, const B & b, const C & c) { return apply(detail::select{}, a, b, c); }

    // Vector algebra functions
    template<class T> constexpr T               cross    (const vec<T,2> & a, const vec<T,2> & b)      { return a.x*b.y-a.y*b.x; }
    template<class T> constexpr vec<T,2>        cross    (T a, const vec<T,2> & b)                     { return {-a*b.y, a*b.x}; }
    template<class T> constexpr vec<T,2>        cross    (const vec<T,2> & a, T b)                     { return {a.y*b, -a.x*b}; }
    template<class T> constexpr vec<T,3>        cross    (const vec<T,3> & a, const vec<T,3> & b)      { return {a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x}; }
    template<class T, int M> constexpr T        dot      (const vec<T,M> & a, const vec<T,M> & b)      { return sum(a*b); }
    template<class T, int M> constexpr T        length2  (const vec<T,M> & a)                          { return dot(a,a); }
    template<class T, int M> T                  length   (const vec<T,M> & a)                          { return std::sqrt(length2(a)); }
    template<class T, int M> vec<T,M>           normalize(const vec<T,M> & a)                          { return a / length(a); }
    template<class T, int M> constexpr T        distance2(const vec<T,M> & a, const vec<T,M> & b)      { return length2(b-a); }
    template<class T, int M> T                  distance (const vec<T,M> & a, const vec<T,M> & b)      { return length(b-a); }
    template<class T, int M> T                  uangle   (const vec<T,M> & a, const vec<T,M> & b)      { T d=dot(a,b); return d > 1 ? 0 : std::acos(d < -1 ? -1 : d); }
    template<class T, int M> T                  angle    (const vec<T,M> & a, const vec<T,M> & b)      { return uangle(normalize(a), normalize(b)); }
    template<class T> vec<T,2>                  rot      (T a, const vec<T,2> & v)                     { const T s = std::sin(a), c = std::cos(a); return {v.x*c - v.y*s, v.x*s + v.y*c}; }
    template<class T, int M> constexpr vec<T,M> lerp     (const vec<T,M> & a, const vec<T,M> & b, T t) { return a*(1-t) + b*t; }
    template<class T, int M> vec<T,M>           nlerp    (const vec<T,M> & a, const vec<T,M> & b, T t) { return normalize(lerp(a,b,t)); }
    template<class T, int M> vec<T,M>           slerp    (const vec<T,M> & a, const vec<T,M> & b, T t) { T th=uangle(a,b); return th == 0 ? a : a*(std::sin(th*(1-t))/std::sin(th)) + b*(std::sin(th*t)/std::sin(th)); }

    ////////////////////////////////////
    // Matrix operators and functions //
    ////////////////////////////////////

    // Unary operators
    template<class T, int M, int N> constexpr mat<T,M,N> operator + (const mat<T,M,N> & a) { return apply(detail::op_pos{}, a); }
    template<class T, int M, int N> constexpr mat<T,M,N> operator - (const mat<T,M,N> & a) { return apply(detail::op_neg{}, a); }

    // Binary operators: matrix $ vector
    template<class T, int M> constexpr auto operator * (const mat<T,M,2> & a, const vec<T,2> & b) { return a[0]*b.x + a[1]*b.y; }
    template<class T, int M> constexpr auto operator * (const mat<T,M,3> & a, const vec<T,3> & b) { return a[0]*b.x + a[1]*b.y + a[2]*b.z; }
    template<class T, int M> constexpr auto operator * (const mat<T,M,4> & a, const vec<T,4> & b) { return a[0]*b.x + a[1]*b.y + a[2]*b.z + a[3]*b.w; }

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
    template<class T, int M> constexpr mat<T,M,2> outerprod   (const vec<T,M> & a, const vec<T,2> & b) { return {a*b.x, a*b.y}; }
    template<class T, int M> constexpr mat<T,M,3> outerprod   (const vec<T,M> & a, const vec<T,3> & b) { return {a*b.x, a*b.y, a*b.z}; }
    template<class T, int M> constexpr mat<T,M,4> outerprod   (const vec<T,M> & a, const vec<T,4> & b) { return {a*b.x, a*b.y, a*b.z, a*b.w}; }
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

    // Vectors and matrices can be used as ranges
    template<class T, int M>       T * begin(      vec<T,M> & a) { return &a[0]; } // Undefined behavior when used as iterator
    template<class T, int M> const T * begin(const vec<T,M> & a) { return &a[0]; } // Undefined behavior when used as iterator
    template<class T, int M>       T * end  (      vec<T,M> & a) { return begin(a) + M; } // Undefined behavior
    template<class T, int M> const T * end  (const vec<T,M> & a) { return begin(a) + M; } // Undefined behavior
    template<class T, int M, int N>       vec<T,M> * begin(      mat<T,M,N> & a) { return a.cols; }
    template<class T, int M, int N> const vec<T,M> * begin(const mat<T,M,N> & a) { return a.cols; }
    template<class T, int M, int N>       vec<T,M> * end  (      mat<T,M,N> & a) { return a.cols + N; }
    template<class T, int M, int N> const vec<T,M> * end  (const mat<T,M,N> & a) { return a.cols + N; }

    // Factory functions for 3D spatial transformations (will possibly be removed or changed in a future version)
    enum fwd_axis { neg_z, pos_z };                 // Should projection matrices be generated assuming forward is {0,0,-1} or {0,0,1}
    enum z_range { neg_one_to_one, zero_to_one };   // Should projection matrices map z into the range of [-1,1] or [0,1]?
    template<class T> vec<T,4>   rotation_quat     (const vec<T,3> & axis, T angle)         { return {axis*std::sin(angle/2), std::cos(angle/2)}; }
    template<class T> vec<T,4>   rotation_quat     (const mat<T,3,3> & m);
    template<class T> mat<T,4,4> translation_matrix(const vec<T,3> & translation)           { return {{1,0,0,0},{0,1,0,0},{0,0,1,0},{translation,1}}; }
    template<class T> mat<T,4,4> rotation_matrix   (const vec<T,4> & rotation)              { return {{qxdir(rotation),0}, {qydir(rotation),0}, {qzdir(rotation),0}, {0,0,0,1}}; }
    template<class T> mat<T,4,4> scaling_matrix    (const vec<T,3> & scaling)               { return {{scaling.x,0,0,0}, {0,scaling.y,0,0}, {0,0,scaling.z,0}, {0,0,0,1}}; }
    template<class T> mat<T,4,4> pose_matrix       (const vec<T,4> & q, const vec<T,3> & p) { return {{qxdir(q),0}, {qydir(q),0}, {qzdir(q),0}, {p,1}}; }
    template<class T> mat<T,4,4> frustum_matrix    (T x0, T x1, T y0, T y1, T n, T f, fwd_axis a = neg_z, z_range z = neg_one_to_one);
    template<class T> mat<T,4,4> perspective_matrix(T fovy, T aspect, T n, T f, fwd_axis a = neg_z, z_range z = neg_one_to_one) { T y = n*std::tan(fovy / 2), x = y*aspect; return frustum_matrix(-x, x, -y, y, n, f, a, z); }

    ////////////////////
    // Legacy support //
    ////////////////////

    // Support for quaternion algebra using 4D vectors, representing xi + yj + zk + w
    template<class T> vec<T,4>           qexp (const vec<T,4> & q)                     { const auto v = q.xyz(); const auto vv = length(v); return std::exp(q.w) * vec<T,4>{v * (vv > 0 ? std::sin(vv)/vv : 0), std::cos(vv)}; }
    template<class T> vec<T,4>           qlog (const vec<T,4> & q)                     { const auto v = q.xyz(); const auto vv = length(v), qq = length(q); return {v * (vv > 0 ? std::acos(q.w/qq)/vv : 0), std::log(qq)}; }
    template<class T> vec<T,4>           qpow (const vec<T,4> & q, const T & p)        { const auto v = q.xyz(); const auto vv = length(v), qq = length(q), th = std::acos(q.w/qq); return std::pow(qq,p)*vec<T,4>{v * (vv > 0 ? std::sin(p*th)/vv : 0), std::cos(p*th)}; }

    // Support for 3D spatial rotations using quaternions, via qmul(qmul(q, v), qconj(q))
    template<class T> constexpr vec<T,3>   qxdir (const vec<T,4> & q)                          { return {q.w*q.w+q.x*q.x-q.y*q.y-q.z*q.z, (q.x*q.y+q.z*q.w)*2, (q.z*q.x-q.y*q.w)*2}; }
    template<class T> constexpr vec<T,3>   qydir (const vec<T,4> & q)                          { return {(q.x*q.y-q.z*q.w)*2, q.w*q.w-q.x*q.x+q.y*q.y-q.z*q.z, (q.y*q.z+q.x*q.w)*2}; }
    template<class T> constexpr vec<T,3>   qzdir (const vec<T,4> & q)                          { return {(q.z*q.x+q.y*q.w)*2, (q.y*q.z-q.x*q.w)*2, q.w*q.w-q.x*q.x-q.y*q.y+q.z*q.z}; }
    template<class T> constexpr mat<T,3,3> qmat  (const vec<T,4> & q)                          { return {qxdir(q), qydir(q), qzdir(q)}; }
    template<class T> constexpr vec<T,3>   qrot  (const vec<T,4> & q, const vec<T,3> & v)      { return qxdir(q)*v.x + qydir(q)*v.y + qzdir(q)*v.z; }
    template<class T> T                    qangle(const vec<T,4> & q)                          { return std::atan2(length(q.xyz()), q.w)*2; }
    template<class T> vec<T,3>             qaxis (const vec<T,4> & q)                          { return normalize(q.xyz()); }
    template<class T> vec<T,4>             qnlerp(const vec<T,4> & a, const vec<T,4> & b, T t) { return nlerp(a, dot(a,b) < 0 ? -b : b, t); }
    template<class T> vec<T,4>             qslerp(const vec<T,4> & a, const vec<T,4> & b, T t) { return slerp(a, dot(a,b) < 0 ? -b : b, t); }

    // These functions exist to ease the difficulty of porting from older versions of linalg
    template<class T> constexpr vec<T,4> qconj(const vec<T,4> & q) { return vec<T,4>(conjugate(quat<T>(q))); }
    template<class T> constexpr vec<T,4> qinv(const vec<T,4> & q) { return vec<T,4>(inverse(quat<T>(q))); }
    template<class T> constexpr vec<T,4> qmul (const vec<T,4> & a, const vec<T,4> & b) { return vec<T,4>(quat<T>(a) * quat<T>(b)); }
    template<class T, class... R> constexpr vec<T,4> qmul(const vec<T,4> & a, R... r)  { return qmul(a, qmul(r...)); }
    template<class T, int M, int N, class B> constexpr auto mul(const mat<T,M,N> & a, const B & b) { return a*b; }
    template<class T, int M, int N, class... R> constexpr auto mul(const mat<T,M,N> & a, R... r) { return mul(a, mul(r...)); }
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

template<class T> linalg::mat<T,4,4> linalg::frustum_matrix(T x0, T x1, T y0, T y1, T n, T f, fwd_axis a, z_range z) 
{ 
    const T s = a == pos_z ? T(1) : T(-1), o = z == neg_one_to_one ? n : 0;
    return {{2*n/(x1-x0),0,0,0}, {0,2*n/(y1-y0),0,0}, {-s*(x0+x1)/(x1-x0),-s*(y0+y1)/(y1-y0),s*(f+o)/(f-n),s}, {0,0,-(n+o)*f/(f-n),0}};
}

#endif
