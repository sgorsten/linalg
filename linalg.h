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



// Option: Define LINALG_LEGACY_MATRIX_PRODUCT to define the matrix product as mul(a,b) instead of a*b.
// Note that either way, the old element-wise definition of a*b has been removed. This can be used to audit
// codebases transitioning from previous versions of linalg.

#pragma once
#ifndef LINALG_H
#define LINALG_H

#include <cmath>        // For various unary math functions, such as std::sqrt
#include <cstdlib>      // To resolve std::abs ambiguity on clang
#include <cstddef>      // For std::nullptr_t
#include <cstdint>      // For implementing namespace linalg::aliases
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
        struct identity_t
        {
            constexpr identity_t(int) {};
            template<class T> constexpr operator mat<T,2,2>() const { return {{1,0},{0,1}}; }
            template<class T> constexpr operator mat<T,3,3>() const { return {{1,0,0},{0,1,0},{0,0,1}}; }
            template<class T> constexpr operator mat<T,4,4>() const { return {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}}; }
            template<class T> constexpr operator quat<T>() const { return {0,0,0,1}; }
        };

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

        // SFINAE helper which determines type of zip(a,b,f) where at least one of a and b is a vec
        template<class A, class B, class F> struct vec_result {};
        template<class T, int M, class F> struct vec_result<vec<T,M>, vec<T,M>, F> { using type = vec<decltype(std::declval<F>()(std::declval<T>(), std::declval<T>())), M>; };
        template<class T, int M, class F> struct vec_result<vec<T,M>, T, F> { using type = vec<decltype(std::declval<F>()(std::declval<T>(), std::declval<T>())), M>; };
        template<class T, int M, class F> struct vec_result<T, vec<T,M>, F> { using type = vec<decltype(std::declval<F>()(std::declval<T>(), std::declval<T>())), M>; };
        template<class A, class B, class F> using vec_result_t = typename vec_result<A,B,F>::type;

        // SFINAE helper which determines type of zip(a,b,f) where at least one of a and b is a vec or mat
        template<class A, class B, class F> struct zip_result : vec_result<A,B,F> {};
        template<class T, int M, int N, class F> struct zip_result<mat<T,M,N>, mat<T,M,N>, F> { using type = vec<decltype(std::declval<F>()(std::declval<T>(), std::declval<T>())), M>; };
        template<class T, int M, int N, class F> struct zip_result<mat<T,M,N>, T, F> { using type = vec<decltype(std::declval<F>()(std::declval<T>(), std::declval<T>())), M>; };
        template<class T, int M, int N, class F> struct zip_result<T, mat<T,M,N>, F> { using type = vec<decltype(std::declval<F>()(std::declval<T>(), std::declval<T>())), M>; };
        template<class A, class B, class F> using zip_result_t = typename zip_result<A,B,F>::type;

        // Lambdas are not constexpr in C++14, so we provide explicit function objects instead
        struct pos { template<class T> constexpr auto operator() (T r) const { return +r; } };
        struct neg { template<class T> constexpr auto operator() (T r) const { return -r; } };
        struct add { template<class T> constexpr auto operator() (T l, T r) const { return l + r; } };
        struct sub { template<class T> constexpr auto operator() (T l, T r) const { return l - r; } };
        struct mul { template<class T> constexpr auto operator() (T l, T r) const { return l * r; } };
        struct div { template<class T> constexpr auto operator() (T l, T r) const { return l / r; } };
        struct mod { template<class T> constexpr auto operator() (T l, T r) const { return l % r; } };
        struct lshift { template<class T> constexpr auto operator() (T l, T r) const { return l << r; } };
        struct rshift { template<class T> constexpr auto operator() (T l, T r) const { return l >> r; } };

        struct bit_not { template<class T> constexpr auto operator() (T r) const { return ~r; } };
        struct bit_or { template<class T> constexpr auto operator() (T l, T r) const { return l | r; } };
        struct bit_xor { template<class T> constexpr auto operator() (T l, T r) const { return l ^ r; } };
        struct bit_and { template<class T> constexpr auto operator() (T l, T r) const { return l & r; } };

        struct bool_not { template<class T> constexpr bool operator() (T r) const { return !r; } };
        struct bool_or { template<class T> constexpr bool operator() (T l, T r) const { return l || r; } };
        struct bool_and { template<class T> constexpr bool operator() (T l, T r) const { return l && r; } };

        // TODO: Decide if we should implement these only in terms of == and <, or keep them passing through to the underlying operator
        struct equal { template<class T> constexpr bool operator() (T l, T r) const { return l == r; } };
        struct nequal { template<class T> constexpr bool operator() (T l, T r) const { return l != r; } };
        struct less { template<class T> constexpr bool operator() (T l, T r) const { return l < r; } };
        struct greater { template<class T> constexpr bool operator() (T l, T r) const { return l > r; } };
        struct lequal { template<class T> constexpr bool operator() (T l, T r) const { return l <= r; } };
        struct gequal { template<class T> constexpr bool operator() (T l, T r) const { return l >= r; } };

        struct min { template<class T> constexpr T operator() (T l, T r) const { return l < r ? l : r; } };
        struct max { template<class T> constexpr T operator() (T l, T r) const { return l < r ? r : l; } };
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
        typedef vec<T,M>            V;
        V                           cols[2];
        constexpr                   mat()                               : cols{} {}
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
        typedef vec<T,M>            V;
        V                           cols[3];
        constexpr                   mat()                               : cols{} {}
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
        typedef vec<T,M>            V;
        V                           cols[4];
        constexpr                   mat()                               : cols{} {}
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

    // Produce a vector/matrix/quaternion by applying f(T,T) to corresponding pairs of elements from vectors/matrices/quaternions a and b
    template<class T,               class F> constexpr auto zip(const vec<T,2  > & a, const vec<T,2  > & b, F f) { return vec<decltype(f(T(),T())),2>{f(a.x,b.x), f(a.y,b.y)}; }
    template<class T,               class F> constexpr auto zip(const vec<T,3  > & a, const vec<T,3  > & b, F f) { return vec<decltype(f(T(),T())),3>{f(a.x,b.x), f(a.y,b.y), f(a.z,b.z)}; }
    template<class T,               class F> constexpr auto zip(const vec<T,4  > & a, const vec<T,4  > & b, F f) { return vec<decltype(f(T(),T())),4>{f(a.x,b.x), f(a.y,b.y), f(a.z,b.z), f(a.w,b.w)}; }
    template<class T, int M,        class F> constexpr auto zip(const vec<T,M  > & a,                  T b, F f) { return zip(a, vec<T,M>(b), f); }
    template<class T, int M,        class F> constexpr auto zip(                 T a, const vec<T,M  > & b, F f) { return zip(vec<T,M>(a), b, f); }

    template<class T, int M,        class F> constexpr auto zip(const mat<T,M,2> & a, const mat<T,M,2> & b, F f) { return mat<decltype(f(T(),T())),M,2>{zip(a[0],b[0],f), zip(a[1],b[1],f)}; }
    template<class T, int M,        class F> constexpr auto zip(const mat<T,M,3> & a, const mat<T,M,3> & b, F f) { return mat<decltype(f(T(),T())),M,3>{zip(a[0],b[0],f), zip(a[1],b[1],f), zip(a[2],b[2],f)}; }
    template<class T, int M,        class F> constexpr auto zip(const mat<T,M,4> & a, const mat<T,M,4> & b, F f) { return mat<decltype(f(T(),T())),M,4>{zip(a[0],b[0],f), zip(a[1],b[1],f), zip(a[2],b[2],f), zip(a[3],b[3],f)}; }
    template<class T, int M, int N, class F> constexpr auto zip(const mat<T,M,N> & a,                  T b, F f) { return zip(a, mat<T,M,N>(b), f); }
    template<class T, int M, int N, class F> constexpr auto zip(                 T a, const mat<T,M,N> & b, F f) { return zip(mat<T,M,N>(a), b, f); }

    // Produce a vector/matrix/quaternion by applying f(T) to elements from vector/matrix/quaternion a
    template<class T, int M,        class F> constexpr auto map(const vec<T,M  > & a, F f) { return zip(a, a, [f](T l, T) { return f(l); }); }
    template<class T, int M, int N, class F> constexpr auto map(const mat<T,M,N> & a, F f) { return zip(a, a, [f](T l, T) { return f(l); }); }

    //////////////////////////////////////
    // Vector-valued operator overloads //
    //////////////////////////////////////

    // Component-wise unary operators: $vector
    template<class T, int M> constexpr auto operator + (const vec<T,M> & a) { return map(a, detail::pos{}); }
    template<class T, int M> constexpr auto operator - (const vec<T,M> & a) { return map(a, detail::neg{}); }
    template<class T, int M> constexpr auto operator ~ (const vec<T,M> & a) { return map(a, detail::bit_not{}); }
    template<class T, int M> constexpr auto operator ! (const vec<T,M> & a) { return map(a, detail::bool_not{}); }

    // Component-wise binary operators: vector $ vector; vector $ scalar; scalar $ vector
    template<class A, class B> constexpr detail::vec_result_t<A,B,detail::add>     operator +  (const A & a, const B & b) { return zip(a, b, detail::add{}); }
    template<class A, class B> constexpr detail::vec_result_t<A,B,detail::sub>     operator -  (const A & a, const B & b) { return zip(a, b, detail::sub{}); }
    template<class A, class B> constexpr detail::vec_result_t<A,B,detail::mul>     operator *  (const A & a, const B & b) { return zip(a, b, detail::mul{}); }
    template<class A, class B> constexpr detail::vec_result_t<A,B,detail::div>     operator /  (const A & a, const B & b) { return zip(a, b, detail::div{}); }
    template<class A, class B> constexpr detail::vec_result_t<A,B,detail::mod>     operator %  (const A & a, const B & b) { return zip(a, b, detail::mod{}); }
    template<class A, class B> constexpr detail::vec_result_t<A,B,detail::bit_or>  operator |  (const A & a, const B & b) { return zip(a, b, detail::bit_or{}); }
    template<class A, class B> constexpr detail::vec_result_t<A,B,detail::bit_xor> operator ^  (const A & a, const B & b) { return zip(a, b, detail::bit_xor{}); }
    template<class A, class B> constexpr detail::vec_result_t<A,B,detail::bit_and> operator &  (const A & a, const B & b) { return zip(a, b, detail::bit_and{}); }
    template<class A, class B> constexpr detail::vec_result_t<A,B,detail::lshift>  operator << (const A & a, const B & b) { return zip(a, b, detail::lshift{}); }
    template<class A, class B> constexpr detail::vec_result_t<A,B,detail::rshift>  operator >> (const A & a, const B & b) { return zip(a, b, detail::rshift{}); }

    // Algebraic binary operators: matrix * vector
    template<class T, int M> constexpr auto operator * (const mat<T,M,2> & a, const vec<T,2> & b) { return a[0]*b.x + a[1]*b.y; }
    template<class T, int M> constexpr auto operator * (const mat<T,M,3> & a, const vec<T,3> & b) { return a[0]*b.x + a[1]*b.y + a[2]*b.z; }
    template<class T, int M> constexpr auto operator * (const mat<T,M,4> & a, const vec<T,4> & b) { return a[0]*b.x + a[1]*b.y + a[2]*b.z + a[3]*b.w; }

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


    //////////////////////////////////////
    // Matrix-valued operator overloads //
    //////////////////////////////////////

    // Unary operators
    template<class T, int M, int N> constexpr auto operator + (const mat<T,M,N> & a) { return map(a, detail::pos{}); }
    template<class T, int M, int N> constexpr auto operator - (const mat<T,M,N> & a) { return map(a, detail::neg{}); }

    // Binary operators
    template<class T, int M, int N> constexpr auto operator + (const mat<T,M,N> & a, const mat<T,M,N> & b) { return zip(a, b, detail::add{}); }
    template<class T, int M, int N> constexpr auto operator - (const mat<T,M,N> & a, const mat<T,M,N> & b) { return zip(a, b, detail::sub{}); }
    template<class T, int M, int N> constexpr auto operator * (const mat<T,M,N> & a, const mat<T,N,2> & b) { return mat<T,M,2>{a*b[0], a*b[1]}; }
    template<class T, int M, int N> constexpr auto operator * (const mat<T,M,N> & a, const mat<T,N,3> & b) { return mat<T,M,3>{a*b[0], a*b[1], a*b[2]}; }
    template<class T, int M, int N> constexpr auto operator * (const mat<T,M,N> & a, const mat<T,N,4> & b) { return mat<T,M,4>{a*b[0], a*b[1], a*b[2], a*b[3]}; }
    template<class T, int M, int N> constexpr auto operator * (const mat<T,M,N> & a, T b) { return zip(a, b, detail::mul{}); }
    template<class T, int M, int N> constexpr auto operator * (T a, const mat<T,M,N> & b) { return zip(a, b, detail::mul{}); }
    template<class T, int M, int N> constexpr auto operator / (const mat<T,M,N> & a, T b) { return zip(a, b, detail::div{}); }

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

    //////////////////////////////////////////
    // Quaternion-valued operator overloads //
    //////////////////////////////////////////

    // Unary operators
    template<class T> constexpr quat<T> operator + (const quat<T> & a) { return {+a.x, +a.y, +a.z, +a.w}; }
    template<class T> constexpr quat<T> operator - (const quat<T> & a) { return {-a.x, -a.y, -a.z, -a.w}; }

    // Binary operators
    template<class T> constexpr quat<T> operator + (const quat<T> & a, const quat<T> & b) { return {a.x+b.x, a.y+b.y, a.z+b.z, a.w+b.w}; }
    template<class T> constexpr quat<T> operator - (const quat<T> & a, const quat<T> & b) { return {a.x-b.x, a.y-b.y, a.z-b.z, a.w-b.w}; }
    template<class T> constexpr quat<T> operator * (const quat<T> & a, const quat<T> & b) { return {a.x*b.w+a.w*b.x+a.y*b.z-a.z*b.y, a.y*b.w+a.w*b.y+a.z*b.x-a.x*b.z, a.z*b.w+a.w*b.z+a.x*b.y-a.y*b.x, a.w*b.w-a.x*b.x-a.y*b.y-a.z*b.z}; }
    template<class T> constexpr quat<T> operator * (const quat<T> & a, T b) { return {a.x*b, a.y*b, a.z*b, a.w*b}; }
    template<class T> constexpr quat<T> operator * (T a, const quat<T> & b) { return {a*b.x, a*b.y, a*b.z, a*b.w}; }
    template<class T> constexpr quat<T> operator / (const quat<T> & a, T b) { return {a.x/b, a.y/b, a.z/b, a.w/b}; }

    // Binary assignment operators
    template<class T> constexpr quat<T> & operator += (quat<T> & a, const quat<T> & b) { return a = a + b; }
    template<class T> constexpr quat<T> & operator -= (quat<T> & a, const quat<T> & b) { return a = a - b; }
    template<class T> constexpr quat<T> & operator *= (quat<T> & a, const quat<T> & b) { return a = a * b; }
    template<class T> constexpr quat<T> & operator *= (quat<T> & a, const T & b) { return a = a * b; }
    template<class T> constexpr quat<T> & operator /= (quat<T> & a, const T & b) { return a = a / b; }

    // Quaternion algebra functions
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

    ///////////////////////
    // Various functions //
    ///////////////////////

    // Reduction functions
    template<class A> constexpr auto any    (const A & a) { return fold(a, detail::bool_or{}); }
    template<class A> constexpr auto all    (const A & a) { return fold(a, detail::bool_and{}); }
    template<class A> constexpr auto sum    (const A & a) { return fold(a, detail::add{}); }
    template<class A> constexpr auto product(const A & a) { return fold(a, detail::mul{}); }
    template<class A> constexpr auto minelem(const A & a) { return fold(a, detail::min{}); }
    template<class A> constexpr auto maxelem(const A & a) { return fold(a, detail::max{}); }

    // Component-wise standard library math functions
    template<class A> auto abs  (const A & a) { return map(a, [](auto l) { return std::abs  (l); }); }
    template<class A> auto floor(const A & a) { return map(a, [](auto l) { return std::floor(l); }); }
    template<class A> auto ceil (const A & a) { return map(a, [](auto l) { return std::ceil (l); }); }
    template<class A> auto exp  (const A & a) { return map(a, [](auto l) { return std::exp  (l); }); }
    template<class A> auto log  (const A & a) { return map(a, [](auto l) { return std::log  (l); }); }
    template<class A> auto log10(const A & a) { return map(a, [](auto l) { return std::log10(l); }); }
    template<class A> auto sqrt (const A & a) { return map(a, [](auto l) { return std::sqrt (l); }); }
    template<class A> auto sin  (const A & a) { return map(a, [](auto l) { return std::sin  (l); }); }
    template<class A> auto cos  (const A & a) { return map(a, [](auto l) { return std::cos  (l); }); }
    template<class A> auto tan  (const A & a) { return map(a, [](auto l) { return std::tan  (l); }); }
    template<class A> auto asin (const A & a) { return map(a, [](auto l) { return std::asin (l); }); }
    template<class A> auto acos (const A & a) { return map(a, [](auto l) { return std::acos (l); }); }
    template<class A> auto atan (const A & a) { return map(a, [](auto l) { return std::atan (l); }); }
    template<class A> auto sinh (const A & a) { return map(a, [](auto l) { return std::sinh (l); }); }
    template<class A> auto cosh (const A & a) { return map(a, [](auto l) { return std::cosh (l); }); }
    template<class A> auto tanh (const A & a) { return map(a, [](auto l) { return std::tanh (l); }); }
    template<class A> auto round(const A & a) { return map(a, [](auto l) { return std::round(l); }); }

    template<class A, class B> auto fmod    (const A & a, const B & b) { return zip(a, b, [](auto l, auto r) { return std::fmod    (l, r); }); }
    template<class A, class B> auto pow     (const A & a, const B & b) { return zip(a, b, [](auto l, auto r) { return std::pow     (l, r); }); }
    template<class A, class B> auto atan2   (const A & a, const B & b) { return zip(a, b, [](auto l, auto r) { return std::atan2   (l, r); }); }
    template<class A, class B> auto copysign(const A & a, const B & b) { return zip(a, b, [](auto l, auto r) { return std::copysign(l, r); }); }

    // Component-wise relational functions
    template<class A, class B> constexpr detail::zip_result_t<A,B,detail::equal  > equal  (const A & a, const B & b) { return zip(a, b, detail::equal  {}); }
    template<class A, class B> constexpr detail::zip_result_t<A,B,detail::nequal > nequal (const A & a, const B & b) { return zip(a, b, detail::nequal {}); }
    template<class A, class B> constexpr detail::zip_result_t<A,B,detail::less   > less   (const A & a, const B & b) { return zip(a, b, detail::less   {}); }
    template<class A, class B> constexpr detail::zip_result_t<A,B,detail::greater> greater(const A & a, const B & b) { return zip(a, b, detail::greater{}); }
    template<class A, class B> constexpr detail::zip_result_t<A,B,detail::lequal > lequal (const A & a, const B & b) { return zip(a, b, detail::lequal {}); }
    template<class A, class B> constexpr detail::zip_result_t<A,B,detail::gequal > gequal (const A & a, const B & b) { return zip(a, b, detail::gequal {}); }

    // Component-wise min, max, and clamp to a range
    template<class A, class B> constexpr auto min  (const A & a, const B & b) { return zip(a, b, detail::min{}); }
    template<class A, class B> constexpr auto max  (const A & a, const B & b) { return zip(a, b, detail::max{}); }
    template<class A, class B> constexpr auto clamp(const A & a, const B & b, const B & c) { return min(max(a,b),c); } // TODO: Revisit

    // Search functions
    template<class T, int M> int argmin(const vec<T,M> & a) { int j=0; for(int i=1; i<M; ++i) if(a[i] < a[j]) j = i; return j; }
    template<class T, int M> int argmax(const vec<T,M> & a) { int j=0; for(int i=1; i<M; ++i) if(a[i] > a[j]) j = i; return j; }

    // Vectors and matrices can be used as ranges
    template<class T, int M>       T * begin(      vec<T,M> & a) { return &a[0]; }
    template<class T, int M> const T * begin(const vec<T,M> & a) { return &a[0]; }
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

    // Provide typedefs for common element types and vector/matrix sizes
    namespace aliases
    {
        typedef vec<bool,2> bool2; typedef vec<uint8_t,2> byte2; typedef vec<int16_t,2> short2; typedef vec<uint16_t,2> ushort2; 
        typedef vec<bool,3> bool3; typedef vec<uint8_t,3> byte3; typedef vec<int16_t,3> short3; typedef vec<uint16_t,3> ushort3; 
        typedef vec<bool,4> bool4; typedef vec<uint8_t,4> byte4; typedef vec<int16_t,4> short4; typedef vec<uint16_t,4> ushort4;
        typedef vec<int,2> int2; typedef vec<unsigned,2> uint2; typedef vec<float,2> float2; typedef vec<double,2> double2;
        typedef vec<int,3> int3; typedef vec<unsigned,3> uint3; typedef vec<float,3> float3; typedef vec<double,3> double3;
        typedef vec<int,4> int4; typedef vec<unsigned,4> uint4; typedef vec<float,4> float4; typedef vec<double,4> double4;
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
    template<class T> [[deprecated("prefer conjugate(quat<T>)")]] constexpr vec<T,4> qconj(const vec<T,4> & q) { return vec<T,4>(conjugate(quat<T>(q))); }
    template<class T> [[deprecated("prefer inverse(quat<T>)")]] constexpr vec<T,4> qinv(const vec<T,4> & q) { return vec<T,4>(inverse(quat<T>(q))); }
    template<class T> [[deprecated("prefer quat<T> * ...")]] constexpr vec<T,4> qmul (const vec<T,4> & a, const vec<T,4> & b) { return vec<T,4>(quat<T>(a) * quat<T>(b)); }
    template<class T, class... R> [[deprecated("prefer quat<T> * ...")]] constexpr vec<T,4> qmul(const vec<T,4> & a, R... r)  { return qmul(a, qmul(r...)); }
    template<class T, int M, int N, class B> [[deprecated("prefer mat<T,M,N> * ...")]] constexpr auto mul(const mat<T,M,N> & a, const B & b) { return a*b; }
    template<class T, int M, int N, class... R> [[deprecated("prefer mat<T,M,N> * ...")]] constexpr auto mul(const mat<T,M,N> & a, R... r) { return mul(a, mul(r...)); }
}

// Provide specializations for std::hash<...> with linalg types
namespace std 
{ 
    template<class T> struct hash<linalg::vec<T,2>> { std::size_t operator()(const linalg::vec<T,2> & v) const { std::hash<T> h; return h(v[0]) ^ (h(v[1]) << 1); } };
    template<class T> struct hash<linalg::vec<T,3>> { std::size_t operator()(const linalg::vec<T,3> & v) const { std::hash<T> h; return h(v[0]) ^ (h(v[1]) << 1) ^ (h(v[2]) << 2); } };
    template<class T> struct hash<linalg::vec<T,4>> { std::size_t operator()(const linalg::vec<T,4> & v) const { std::hash<T> h; return h(v[0]) ^ (h(v[1]) << 1) ^ (h(v[2]) << 2) ^ (h(v[3]) << 3); } };
    template<class T, int M> struct hash<linalg::mat<T,M,2>> { std::size_t operator()(const linalg::mat<T,M,2> & v) const { std::hash<linalg::vec<T,M>> h; return h(v[0]) ^ (h(v[1]) << M); } };
    template<class T, int M> struct hash<linalg::mat<T,M,3>> { std::size_t operator()(const linalg::mat<T,M,3> & v) const { std::hash<linalg::vec<T,M>> h; return h(v[0]) ^ (h(v[1]) << M) ^ (h(v[2]) << (M*2)); } };
    template<class T, int M> struct hash<linalg::mat<T,M,4>> { std::size_t operator()(const linalg::mat<T,M,4> & v) const { std::hash<linalg::vec<T,M>> h; return h(v[0]) ^ (h(v[1]) << M) ^ (h(v[2]) << (M*2)) ^ (h(v[3]) << (M*3)); } };
}

// Definitions of functions too long to be defined inline
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
