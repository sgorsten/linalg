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

#include <cmath>        // For various unary math functions, such as std::sqrt
#include <cstdlib>      // To resolve std::abs ambiguity on clang
#include <cstdint>      // For implementing namespace linalg::aliases
#include <type_traits>  // For std::enable_if, std::is_same, std::declval

// Visual Studio versions prior to 2015 lack constexpr support
#if defined(_MSC_VER) && _MSC_VER < 1900 && !defined(constexpr)
#define constexpr
#endif

namespace linalg
{
    // Specialize converter<T,U> with a function application operator that converts type U to type T to enable implicit conversions
    template<class T, class U> struct converter {};
    namespace detail { template<class T, class U> using conv_t = typename std::enable_if<!std::is_same<T,U>::value, decltype(converter<T,U>{}(std::declval<U>()))>::type; }

    // Small, fixed-length vector type, consisting of exactly M elements of type T, and presumed to be a column-vector unless otherwise noted
    template<class T, int M> struct vec;
    template<class T> struct vec<T,2>
    {
        T                           x,y;
        constexpr                   vec()                               : x(), y() {}
        constexpr                   vec(const T & x_, const T & y_)     : x(x_), y(y_) {}
        constexpr explicit          vec(const T & s)                    : vec(s, s) {}
        constexpr explicit          vec(const T * p)                    : vec(p[0], p[1]) {}
        template<class U>
        constexpr explicit          vec(const vec<U,2> & v)             : vec(static_cast<T>(v.x), static_cast<T>(v.y)) {}
        constexpr const T &         operator[] (int i) const            { return (&x)[i]; }
        T &                         operator[] (int i)                  { return (&x)[i]; }

        template<class U, class=detail::conv_t<vec,U>> constexpr vec(const U & u) : vec(converter<vec,U>{}(u)) {}
        template<class U, class=detail::conv_t<U,vec>> constexpr operator U () const { return converter<U,vec>{}(*this); }
    };
    template<class T> struct vec<T,3>
    {
        T                           x,y,z;
        constexpr                   vec()                               : x(), y(), z() {}
        constexpr                   vec(const T & x_, const T & y_, 
                                        const T & z_)                   : x(x_), y(y_), z(z_) {}
        constexpr                   vec(const vec<T,2> & xy,
                                        const T & z_)                   : vec(xy.x, xy.y, z_) {}
        constexpr explicit          vec(const T & s)                    : vec(s, s, s) {}
        constexpr explicit          vec(const T * p)                    : vec(p[0], p[1], p[2]) {}
        template<class U>
        constexpr explicit          vec(const vec<U,3> & v)             : vec(static_cast<T>(v.x), static_cast<T>(v.y), static_cast<T>(v.z)) {}
        constexpr const T &         operator[] (int i) const            { return (&x)[i]; }
        T &                         operator[] (int i)                  { return (&x)[i]; }
        constexpr const vec<T,2> &  xy() const                          { return *reinterpret_cast<const vec<T,2> *>(this); }
        vec<T,2> &                  xy()                                { return *reinterpret_cast<vec<T,2> *>(this); }

        template<class U, class=detail::conv_t<vec,U>> constexpr vec(const U & u) : vec(converter<vec,U>{}(u)) {}
        template<class U, class=detail::conv_t<U,vec>> constexpr operator U () const { return converter<U,vec>{}(*this); }
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
        constexpr explicit          vec(const T & s)                    : vec(s, s, s, s) {}
        constexpr explicit          vec(const T * p)                    : vec(p[0], p[1], p[2], p[3]) {}
        template<class U> 
        constexpr explicit          vec(const vec<U,4> & v)             : vec(static_cast<T>(v.x), static_cast<T>(v.y), static_cast<T>(v.z), static_cast<T>(v.w)) {}
        constexpr const T &         operator[] (int i) const            { return (&x)[i]; }
        T &                         operator[] (int i)                  { return (&x)[i]; }
        constexpr const vec<T,2> &  xy() const                          { return *reinterpret_cast<const vec<T,2> *>(this); }
        constexpr const vec<T,3> &  xyz() const                         { return *reinterpret_cast<const vec<T,3> *>(this); }
        vec<T,2> &                  xy()                                { return *reinterpret_cast<vec<T,2> *>(this); }                
        vec<T,3> &                  xyz()                               { return *reinterpret_cast<vec<T,3> *>(this); }

        template<class U, class=detail::conv_t<vec,U>> constexpr vec(const U & u) : vec(converter<vec,U>{}(u)) {}
        template<class U, class=detail::conv_t<U,vec>> constexpr operator U () const { return converter<U,vec>{}(*this); }
    };

    // Small, fixed-size matrix type, consisting of exactly M rows and N columns of type T, stored in column-major order.
    template<class T, int M, int N> struct mat;
    template<class T, int M> struct mat<T,M,2>
    {
        typedef vec<T,M>            V;
        V                           x,y;
        constexpr                   mat()                               : x(), y() {}
        constexpr                   mat(const V & x_, const V & y_)     : x(x_), y(y_) {}
        constexpr explicit          mat(const T & s)                    : x(s), y(s) {}
        constexpr explicit          mat(const T * p)                    : x(p+M*0), y(p+M*1) {}
        template<class U> 
        constexpr explicit          mat(const mat<U,M,2> & m)           : mat(V(m.x), V(m.y)) {}
        constexpr vec<T,2>          row(int i) const                    { return {x[i], y[i]}; }
        constexpr const V &         operator[] (int j) const            { return (&x)[j]; }
        V &                         operator[] (int j)                  { return (&x)[j]; }

        template<class U, class=detail::conv_t<mat,U>> constexpr mat(const U & u) : mat(converter<mat,U>{}(u)) {}
        template<class U, class=detail::conv_t<U,mat>> constexpr operator U () const { return converter<U,mat>{}(*this); }
    };
    template<class T, int M> struct mat<T,M,3>
    {
        typedef vec<T,M>            V;
        V                           x,y,z;
        constexpr                   mat()                               : x(), y(), z() {}
        constexpr                   mat(const V & x_, const V & y_, 
                                        const V & z_)                   : x(x_), y(y_), z(z_) {}
        constexpr explicit          mat(const T & s)                    : x(s), y(s), z(s) {}
        constexpr explicit          mat(const T * p)                    : x(p+M*0), y(p+M*1), z(p+M*2) {}
        template<class U> 
        constexpr explicit          mat(const mat<U,M,3> & m)           : mat(V(m.x), V(m.y), V(m.z)) {}
        constexpr vec<T,3>          row(int i) const                    { return {x[i], y[i], z[i]}; }
        constexpr const V &         operator[] (int j) const            { return (&x)[j]; }
        V &                         operator[] (int j)                  { return (&x)[j]; }

        template<class U, class=detail::conv_t<mat,U>> constexpr mat(const U & u) : mat(converter<mat,U>{}(u)) {}
        template<class U, class=detail::conv_t<U,mat>> constexpr operator U () const { return converter<U,mat>{}(*this); }
    };
    template<class T, int M> struct mat<T,M,4>
    {
        typedef vec<T,M>            V;
        V                           x,y,z,w;
        constexpr                   mat()                               : x(), y(), z(), w() {}
        constexpr                   mat(const V & x_, const V & y_,
                                        const V & z_, const V & w_)     : x(x_), y(y_), z(z_), w(w_) {}
        constexpr explicit          mat(const T & s)                    : x(s), y(s), z(s), w(s) {}
        constexpr explicit          mat(const T * p)                    : x(p+M*0), y(p+M*1), z(p+M*2), w(p+M*3) {}
        template<class U> 
        constexpr explicit          mat(const mat<U,M,4> & m)           : mat(V(m.x), V(m.y), V(m.z), V(m.w)) {}
        constexpr vec<T,4>          row(int i) const                    { return {x[i], y[i], z[i], w[i]}; }
        constexpr const V &         operator[] (int j) const            { return (&x)[j]; }
        V &                         operator[] (int j)                  { return (&x)[j]; }

        template<class U, class=detail::conv_t<mat,U>> constexpr mat(const U & u) : mat(converter<mat,U>{}(u)) {}
        template<class U, class=detail::conv_t<U,mat>> constexpr operator U () const { return converter<U,mat>{}(*this); }
    };

    // Define a type which will convert to the multiplicative identity of any square matrix
    struct identity_t { constexpr explicit identity_t(int) {} };
    template<class T> struct converter<mat<T,2,2>, identity_t> { mat<T,2,2> operator() (identity_t) const { return {{1,0},{0,1}}; } };
    template<class T> struct converter<mat<T,3,3>, identity_t> { mat<T,3,3> operator() (identity_t) const { return {{1,0,0},{0,1,0},{0,0,1}}; } };
    template<class T> struct converter<mat<T,4,4>, identity_t> { mat<T,4,4> operator() (identity_t) const { return {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}}; } };
    constexpr identity_t identity {1};

    // Type traits for a binary operation involving linear algebra types, used for SFINAE on templated functions and operator overloads
    template<class A, class B> struct traits {};
    template<class T, int M       > struct traits<vec<T,M  >, vec<T,M  >> { typedef T scalar; typedef vec<T,M  > result; typedef vec<bool,M  > bool_result; typedef vec<decltype(+T()),M  > arith_result; };
    template<class T, int M       > struct traits<vec<T,M  >, T         > { typedef T scalar; typedef vec<T,M  > result; typedef vec<bool,M  > bool_result; typedef vec<decltype(+T()),M  > arith_result; };
    template<class T, int M       > struct traits<T,          vec<T,M  >> { typedef T scalar; typedef vec<T,M  > result; typedef vec<bool,M  > bool_result; typedef vec<decltype(+T()),M  > arith_result; };
    template<class T, int M, int N> struct traits<mat<T,M,N>, mat<T,M,N>> { typedef T scalar; typedef mat<T,M,N> result; typedef mat<bool,M,N> bool_result; typedef mat<decltype(+T()),M,N> arith_result; };
    template<class T, int M, int N> struct traits<mat<T,M,N>, T         > { typedef T scalar; typedef mat<T,M,N> result; typedef mat<bool,M,N> bool_result; typedef mat<decltype(+T()),M,N> arith_result; };
    template<class T, int M, int N> struct traits<T,          mat<T,M,N>> { typedef T scalar; typedef mat<T,M,N> result; typedef mat<bool,M,N> bool_result; typedef mat<decltype(+T()),M,N> arith_result; };
    template<class A, class B=A> using scalar_t = typename traits<A,B>::scalar; // Underlying scalar type when performing elementwise operations
    template<class A, class B=A> using result_t = typename traits<A,B>::result; // Result of calling a function on linear algebra types
    template<class A, class B=A> using bool_result_t = typename traits<A,B>::bool_result; // Result of a comparison or unary not operation on linear algebra types
    template<class A, class B=A> using arith_result_t = typename traits<A,B>::arith_result; // Result of an arithmetic operation on linear algebra types (accounts for integer promotion)

    // Produce a scalar by applying f(T,T) -> T to adjacent pairs of elements from vector/matrix a in left-to-right order (matching the associativity of arithmetic and logical operators)
    template<class T, class F> constexpr T fold(const vec<T,2> & a, F f) { return f(a.x,a.y); }
    template<class T, class F> constexpr T fold(const vec<T,3> & a, F f) { return f(f(a.x,a.y),a.z); }
    template<class T, class F> constexpr T fold(const vec<T,4> & a, F f) { return f(f(f(a.x,a.y),a.z),a.w); }
    template<class T, int M, class F> constexpr T fold(const mat<T,M,2> & a, F f) { return f(fold(a.x,f),fold(a.y,f)); }
    template<class T, int M, class F> constexpr T fold(const mat<T,M,3> & a, F f) { return f(f(fold(a.x,f),fold(a.y,f)),fold(a.z,f)); }
    template<class T, int M, class F> constexpr T fold(const mat<T,M,4> & a, F f) { return f(f(f(fold(a.x,f),fold(a.y,f)),fold(a.z,f)),fold(a.w,f)); }

    // Produce a vector/matrix by applying f(T,T) to corresponding pairs of elements from vectors/matrix a and b
    template<class T,               class F> constexpr auto zip(const vec<T,2  > & a, const vec<T,2  > & b, F f) -> vec<decltype(f(T(),T())),2  > { return {f(a.x,b.x), f(a.y,b.y)}; }
    template<class T,               class F> constexpr auto zip(const vec<T,3  > & a, const vec<T,3  > & b, F f) -> vec<decltype(f(T(),T())),3  > { return {f(a.x,b.x), f(a.y,b.y), f(a.z,b.z)}; }
    template<class T,               class F> constexpr auto zip(const vec<T,4  > & a, const vec<T,4  > & b, F f) -> vec<decltype(f(T(),T())),4  > { return {f(a.x,b.x), f(a.y,b.y), f(a.z,b.z), f(a.w,b.w)}; }
    template<class T, int M,        class F> constexpr auto zip(const vec<T,M  > & a,                  T b, F f) -> vec<decltype(f(T(),T())),M  > { return zip(a, vec<T,M>(b), f); }
    template<class T, int M,        class F> constexpr auto zip(                 T a, const vec<T,M  > & b, F f) -> vec<decltype(f(T(),T())),M  > { return zip(vec<T,M>(a), b, f); }
    template<class T, int M,        class F> constexpr auto zip(const mat<T,M,2> & a, const mat<T,M,2> & b, F f) -> mat<decltype(f(T(),T())),M,2> { return {zip(a.x,b.x,f), zip(a.y,b.y,f)}; }
    template<class T, int M,        class F> constexpr auto zip(const mat<T,M,3> & a, const mat<T,M,3> & b, F f) -> mat<decltype(f(T(),T())),M,3> { return {zip(a.x,b.x,f), zip(a.y,b.y,f), zip(a.z,b.z,f)}; }
    template<class T, int M,        class F> constexpr auto zip(const mat<T,M,4> & a, const mat<T,M,4> & b, F f) -> mat<decltype(f(T(),T())),M,4> { return {zip(a.x,b.x,f), zip(a.y,b.y,f), zip(a.z,b.z,f), zip(a.w,b.w,f)}; }
    template<class T, int M, int N, class F> constexpr auto zip(const mat<T,M,N> & a,                  T b, F f) -> mat<decltype(f(T(),T())),M,N> { return zip(a, mat<T,M,N>(b), f); }
    template<class T, int M, int N, class F> constexpr auto zip(                 T a, const mat<T,M,N> & b, F f) -> mat<decltype(f(T(),T())),M,N> { return zip(mat<T,M,N>(a), b, f); }

    // Produce a vector/matrix by applying f(T) to elements from vector/matrix a
    template<class T, int M,        class F> constexpr auto map(const vec<T,M  > & a, F f) -> vec<decltype(f(T())),M  > { return zip(a, a, [f](T l, T) { return f(l); }); }
    template<class T, int M, int N, class F> constexpr auto map(const mat<T,M,N> & a, F f) -> mat<decltype(f(T())),M,N> { return zip(a, a, [f](T l, T) { return f(l); }); }

    // Relational operators are defined to compare the elements of two vectors or matrices lexicographically, in column-major order
    namespace detail
    {
        // Type returned by the compare(...) function which supports all six comparison operators against 0
        template<class T> struct ord { T a,b; };
        template<class T> constexpr bool operator == (const ord<T> & o, std::nullptr_t) { return o.a == o.b; }
        template<class T> constexpr bool operator != (const ord<T> & o, std::nullptr_t) { return !(o.a == o.b); }
        template<class T> constexpr bool operator < (const ord<T> & o, std::nullptr_t) { return o.a < o.b; }
        template<class T> constexpr bool operator > (const ord<T> & o, std::nullptr_t) { return o.b < o.a; }
        template<class T> constexpr bool operator <= (const ord<T> & o, std::nullptr_t) { return !(o.b < o.a); }
        template<class T> constexpr bool operator >= (const ord<T> & o, std::nullptr_t) { return !(o.a < o.b); }

        template<class A, class B> struct any_compare {};
        template<class T> struct any_compare<vec<T,2>,vec<T,2>> { using type=ord<T>; constexpr ord<T> operator() (const vec<T,2> & a, const vec<T,2> & b) const { return !(a[0]==b[0]) ? ord<T>{a[0],b[0]} : ord<T>{a[1],b[1]}; } };
        template<class T> struct any_compare<vec<T,3>,vec<T,3>> { using type=ord<T>; constexpr ord<T> operator() (const vec<T,3> & a, const vec<T,3> & b) const { return !(a[0]==b[0]) ? ord<T>{a[0],b[0]} : !(a[1]==b[1]) ? ord<T>{a[1],b[1]} : ord<T>{a[2],b[2]}; } };
        template<class T> struct any_compare<vec<T,4>,vec<T,4>> { using type=ord<T>; constexpr ord<T> operator() (const vec<T,4> & a, const vec<T,4> & b) const { return !(a[0]==b[0]) ? ord<T>{a[0],b[0]} : !(a[1]==b[1]) ? ord<T>{a[1],b[1]} : !(a[2]==b[2]) ? ord<T>{a[2],b[2]} : ord<T>{a[3],b[3]}; } };
        template<class T, int M> struct any_compare<mat<T,M,2>,mat<T,M,2>> { using type=ord<T>; constexpr ord<T> operator() (const mat<T,M,2> & a, const mat<T,M,2> & b) const { return a[0]!=b[0] ? compare(a[0],b[0]) : compare(a[1],b[1]); } };
        template<class T, int M> struct any_compare<mat<T,M,3>,mat<T,M,3>> { using type=ord<T>; constexpr ord<T> operator() (const mat<T,M,3> & a, const mat<T,M,3> & b) const { return a[0]!=b[0] ? compare(a[0],b[0]) : a[1]!=b[1] ? compare(a[1],b[1]) : compare(a[2],b[2]); } };
        template<class T, int M> struct any_compare<mat<T,M,4>,mat<T,M,4>> { using type=ord<T>; constexpr ord<T> operator() (const mat<T,M,4> & a, const mat<T,M,4> & b) const { return a[0]!=b[0] ? compare(a[0],b[0]) : a[1]!=b[1] ? compare(a[1],b[1]) : a[2]!=b[2] ? compare(a[2],b[2]) : compare(a[3],b[3]); } };
    }
    template<class A, class B> constexpr typename detail::any_compare<A,B>::type compare(const A & a, const B & b) { return detail::any_compare<A,B>()(a,b); }
    template<class A, class B> constexpr auto operator == (const A & a, const B & b) -> decltype(compare(a,b) == 0) { return compare(a,b) == 0; }
    template<class A, class B> constexpr auto operator != (const A & a, const B & b) -> decltype(compare(a,b) != 0) { return compare(a,b) != 0; }
    template<class A, class B> constexpr auto operator <  (const A & a, const B & b) -> decltype(compare(a,b) <  0) { return compare(a,b) <  0; }
    template<class A, class B> constexpr auto operator >  (const A & a, const B & b) -> decltype(compare(a,b) >  0) { return compare(a,b) >  0; }
    template<class A, class B> constexpr auto operator <= (const A & a, const B & b) -> decltype(compare(a,b) <= 0) { return compare(a,b) <= 0; }
    template<class A, class B> constexpr auto operator >= (const A & a, const B & b) -> decltype(compare(a,b) >= 0) { return compare(a,b) >= 0; }

    // Lambdas are not permitted inside constexpr functions, so we provide explicit function objects instead
    namespace op
    {
        template<class T> struct pos { constexpr auto operator() (T r) const -> decltype(+r) { return +r; } };
        template<class T> struct neg { constexpr auto operator() (T r) const -> decltype(-r) { return -r; } };
        template<class T> struct add { constexpr auto operator() (T l, T r) const -> decltype(l + r) { return l + r; } };
        template<class T> struct sub { constexpr auto operator() (T l, T r) const -> decltype(l - r) { return l - r; } };
        template<class T> struct mul { constexpr auto operator() (T l, T r) const -> decltype(l * r) { return l * r; } };
        template<class T> struct div { constexpr auto operator() (T l, T r) const -> decltype(l / r) { return l / r; } };
        template<class T> struct mod { constexpr auto operator() (T l, T r) const -> decltype(l % r) { return l % r; } };
        template<class T> struct lshift { constexpr auto operator() (T l, T r) const -> decltype(l << r) { return l << r; } };
        template<class T> struct rshift { constexpr auto operator() (T l, T r) const -> decltype(l >> r) { return l >> r; } };        

        template<class T> struct binary_not { constexpr auto operator() (T r) const -> decltype(+r) { return ~r; } };
        template<class T> struct binary_or  { constexpr auto operator() (T l, T r) const -> decltype(l | r) { return l | r; } };
        template<class T> struct binary_xor { constexpr auto operator() (T l, T r) const -> decltype(l ^ r) { return l ^ r; } };
        template<class T> struct binary_and { constexpr auto operator() (T l, T r) const -> decltype(l & r) { return l & r; } };

        template<class T> struct logical_not { constexpr bool operator() (T r) const { return !r; } };
        template<class T> struct logical_or  { constexpr bool operator() (T l, T r) const { return l || r; } };
        template<class T> struct logical_and { constexpr bool operator() (T l, T r) const { return l && r; } };        

        template<class T> struct equal   { constexpr bool operator() (T l, T r) const { return l == r; } };
        template<class T> struct nequal  { constexpr bool operator() (T l, T r) const { return l != r; } };
        template<class T> struct less    { constexpr bool operator() (T l, T r) const { return l <  r; } };
        template<class T> struct greater { constexpr bool operator() (T l, T r) const { return l >  r; } };
        template<class T> struct lequal  { constexpr bool operator() (T l, T r) const { return l <= r; } };
        template<class T> struct gequal  { constexpr bool operator() (T l, T r) const { return l >= r; } };

        template<class T> struct min { constexpr T operator() (T l, T r) const { return l < r ? l : r; } };
        template<class T> struct max { constexpr T operator() (T l, T r) const { return l > r ? l : r; } };
    }

    // Functions for coalescing scalar values
    template<class A> constexpr scalar_t<A> any    (const A & a) { return fold(a, op::logical_or<scalar_t<A>>{}); }
    template<class A> constexpr scalar_t<A> all    (const A & a) { return fold(a, op::logical_and<scalar_t<A>>{}); }
    template<class A> constexpr scalar_t<A> sum    (const A & a) { return fold(a, op::add<scalar_t<A>>{}); }
    template<class A> constexpr scalar_t<A> product(const A & a) { return fold(a, op::mul<scalar_t<A>>{}); }
    template<class T, int M> int argmin(const vec<T,M> & a) { int j=0; for(int i=1; i<M; ++i) if(a[i] < a[j]) j = i; return j; }
    template<class T, int M> int argmax(const vec<T,M> & a) { int j=0; for(int i=1; i<M; ++i) if(a[i] > a[j]) j = i; return j; }
    template<class T, int M> T minelem(const vec<T,M> & a) { return a[argmin(a)]; }
    template<class T, int M> T maxelem(const vec<T,M> & a) { return a[argmax(a)]; }

    // Overloads for unary operators on vectors are implemented in terms of elementwise application of the operator
    template<class A> constexpr arith_result_t<A> operator + (const A & a) { return map(a, op::pos<scalar_t<A>>{}); }
    template<class A> constexpr arith_result_t<A> operator - (const A & a) { return map(a, op::neg<scalar_t<A>>{}); }
    template<class A> constexpr arith_result_t<A> operator ~ (const A & a) { return map(a, op::binary_not<scalar_t<A>>{}); }
    template<class A> constexpr bool_result_t<A> operator ! (const A & a) { return map(a, op::logical_not<scalar_t<A>>{}); }

    // Mirror the set of unary scalar math functions to apply elementwise to vectors
    template<class A> result_t<A> abs  (const A & a) { return map(a, [](scalar_t<A> l) { return std::abs  (l); }); }
    template<class A> result_t<A> floor(const A & a) { return map(a, [](scalar_t<A> l) { return std::floor(l); }); }
    template<class A> result_t<A> ceil (const A & a) { return map(a, [](scalar_t<A> l) { return std::ceil (l); }); }
    template<class A> result_t<A> exp  (const A & a) { return map(a, [](scalar_t<A> l) { return std::exp  (l); }); }
    template<class A> result_t<A> log  (const A & a) { return map(a, [](scalar_t<A> l) { return std::log  (l); }); }
    template<class A> result_t<A> log10(const A & a) { return map(a, [](scalar_t<A> l) { return std::log10(l); }); }
    template<class A> result_t<A> sqrt (const A & a) { return map(a, [](scalar_t<A> l) { return std::sqrt (l); }); }
    template<class A> result_t<A> sin  (const A & a) { return map(a, [](scalar_t<A> l) { return std::sin  (l); }); }
    template<class A> result_t<A> cos  (const A & a) { return map(a, [](scalar_t<A> l) { return std::cos  (l); }); }
    template<class A> result_t<A> tan  (const A & a) { return map(a, [](scalar_t<A> l) { return std::tan  (l); }); }
    template<class A> result_t<A> asin (const A & a) { return map(a, [](scalar_t<A> l) { return std::asin (l); }); }
    template<class A> result_t<A> acos (const A & a) { return map(a, [](scalar_t<A> l) { return std::acos (l); }); }
    template<class A> result_t<A> atan (const A & a) { return map(a, [](scalar_t<A> l) { return std::atan (l); }); }
    template<class A> result_t<A> sinh (const A & a) { return map(a, [](scalar_t<A> l) { return std::sinh (l); }); }
    template<class A> result_t<A> cosh (const A & a) { return map(a, [](scalar_t<A> l) { return std::cosh (l); }); }
    template<class A> result_t<A> tanh (const A & a) { return map(a, [](scalar_t<A> l) { return std::tanh (l); }); }
    template<class A> result_t<A> round(const A & a) { return map(a, [](scalar_t<A> l) { return std::round(l); }); }
    template<class A> result_t<A> fract(const A & a) { return map(a, [](scalar_t<A> l) { return l - std::floor(l); }); }

    // Overloads for vector op vector are implemented in terms of elementwise application of the operator, followed by casting back to the original type (integer promotion is suppressed)
    template<class A, class B> constexpr arith_result_t<A,B> operator +  (const A & a, const B & b) { return zip(a, b, op::add<scalar_t<A,B>>{}); }
    template<class A, class B> constexpr arith_result_t<A,B> operator -  (const A & a, const B & b) { return zip(a, b, op::sub<scalar_t<A,B>>{}); }
    template<class A, class B> constexpr arith_result_t<A,B> operator *  (const A & a, const B & b) { return zip(a, b, op::mul<scalar_t<A,B>>{}); }
    template<class A, class B> constexpr arith_result_t<A,B> operator /  (const A & a, const B & b) { return zip(a, b, op::div<scalar_t<A,B>>{}); }
    template<class A, class B> constexpr arith_result_t<A,B> operator %  (const A & a, const B & b) { return zip(a, b, op::mod<scalar_t<A,B>>{}); }
    template<class A, class B> constexpr arith_result_t<A,B> operator |  (const A & a, const B & b) { return zip(a, b, op::binary_or<scalar_t<A,B>>{}); }
    template<class A, class B> constexpr arith_result_t<A,B> operator ^  (const A & a, const B & b) { return zip(a, b, op::binary_xor<scalar_t<A,B>>{}); }
    template<class A, class B> constexpr arith_result_t<A,B> operator &  (const A & a, const B & b) { return zip(a, b, op::binary_and<scalar_t<A,B>>{}); }
    template<class A, class B> constexpr arith_result_t<A,B> operator << (const A & a, const B & b) { return zip(a, b, op::lshift<scalar_t<A,B>>{}); }
    template<class A, class B> constexpr arith_result_t<A,B> operator >> (const A & a, const B & b) { return zip(a, b, op::rshift<scalar_t<A,B>>{}); }
    
    // Overloads for assignment operators are implemented trivially
    template<class A, class B> result_t<A,A> & operator +=  (A & a, const B & b) { return a = a + b; }
    template<class A, class B> result_t<A,A> & operator -=  (A & a, const B & b) { return a = a - b; }
    template<class A, class B> result_t<A,A> & operator *=  (A & a, const B & b) { return a = a * b; }
    template<class A, class B> result_t<A,A> & operator /=  (A & a, const B & b) { return a = a / b; }
    template<class A, class B> result_t<A,A> & operator %=  (A & a, const B & b) { return a = a % b; }
    template<class A, class B> result_t<A,A> & operator |=  (A & a, const B & b) { return a = a | b; }
    template<class A, class B> result_t<A,A> & operator ^=  (A & a, const B & b) { return a = a ^ b; }
    template<class A, class B> result_t<A,A> & operator &=  (A & a, const B & b) { return a = a & b; }
    template<class A, class B> result_t<A,A> & operator <<= (A & a, const B & b) { return a = a << b; }
    template<class A, class B> result_t<A,A> & operator >>= (A & a, const B & b) { return a = a >> b; }

    // Mirror the set of binary scalar math functions to apply elementwise to vectors
    template<class A, class B> constexpr result_t<A,B> min  (const A & a, const B & b) { return zip(a, b, op::min<scalar_t<A,B>>{}); }
    template<class A, class B> constexpr result_t<A,B> max  (const A & a, const B & b) { return zip(a, b, op::max<scalar_t<A,B>>{}); }
    template<class A, class B> constexpr result_t<A,B> clamp(const A & a, const B & b, const B & c) { return min(max(a,b),c); } // TODO: Revisit
    template<class A, class B> result_t<A,B> fmod    (const A & a, const B & b) { return zip(a, b, [](scalar_t<A,B> l, scalar_t<A,B> r) { return std::fmod    (l, r); }); }
    template<class A, class B> result_t<A,B> pow     (const A & a, const B & b) { return zip(a, b, [](scalar_t<A,B> l, scalar_t<A,B> r) { return std::pow     (l, r); }); }
    template<class A, class B> result_t<A,B> atan2   (const A & a, const B & b) { return zip(a, b, [](scalar_t<A,B> l, scalar_t<A,B> r) { return std::atan2   (l, r); }); }
    template<class A, class B> result_t<A,B> copysign(const A & a, const B & b) { return zip(a, b, [](scalar_t<A,B> l, scalar_t<A,B> r) { return std::copysign(l, r); }); }

    // Functions for componentwise application of equivalence and relational operators
    template<class A, class B> bool_result_t<A,B> equal  (const A & a, const B & b) { return zip(a, b, op::equal  <scalar_t<A,B>>{}); }
    template<class A, class B> bool_result_t<A,B> nequal (const A & a, const B & b) { return zip(a, b, op::nequal <scalar_t<A,B>>{}); }
    template<class A, class B> bool_result_t<A,B> less   (const A & a, const B & b) { return zip(a, b, op::less   <scalar_t<A,B>>{}); }
    template<class A, class B> bool_result_t<A,B> greater(const A & a, const B & b) { return zip(a, b, op::greater<scalar_t<A,B>>{}); }
    template<class A, class B> bool_result_t<A,B> lequal (const A & a, const B & b) { return zip(a, b, op::lequal <scalar_t<A,B>>{}); }
    template<class A, class B> bool_result_t<A,B> gequal (const A & a, const B & b) { return zip(a, b, op::gequal <scalar_t<A,B>>{}); }

    // Support for vector algebra
    template<class T> constexpr T                 cross    (const vec<T,2> & a, const vec<T,2> & b)      { return a.x*b.y-a.y*b.x; }
    template<class T> constexpr vec<T,2>          cross    (T a, const vec<T,2> & b)                     { return {-a*b.y, a*b.x}; }
    template<class T> constexpr vec<T,2>          cross    (const vec<T,2> & a, T b)                     { return {a.y*b, -a.x*b}; }
    template<class T> constexpr vec<T,3>          cross    (const vec<T,3> & a, const vec<T,3> & b)      { return {a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x}; }
    template<class T, int M> constexpr T          dot      (const vec<T,M> & a, const vec<T,M> & b)      { return sum(a*b); }
    template<class T, int M> constexpr T          length2  (const vec<T,M> & a)                          { return dot(a,a); }
    template<class T, int M> T                    length   (const vec<T,M> & a)                          { return std::sqrt(length2(a)); }
    template<class T, int M> vec<T,M>             normalize(const vec<T,M> & a)                          { return a / length(a); }
    template<class T, int M> constexpr T          distance2(const vec<T,M> & a, const vec<T,M> & b)      { return length2(b-a); }
    template<class T, int M> T                    distance (const vec<T,M> & a, const vec<T,M> & b)      { return length(b-a); }
    template<class T, int M> T                    uangle   (const vec<T,M> & a, const vec<T,M> & b)      { T d=dot(a,b); return d > 1 ? 0 : std::acos(d < -1 ? -1 : d); }
    template<class T, int M> T                    angle    (const vec<T,M> & a, const vec<T,M> & b)      { return uangle(normalize(a), normalize(b)); }
    template<class T> vec<T,2>                    rot      (T a, const vec<T,2> & v)                     { const T s = std::sin(a), c = std::cos(a); return {v.x*c - v.y*s, v.x*s + v.y*c}; }
    template<class T, int M> constexpr vec<T,M>   lerp     (const vec<T,M> & a, const vec<T,M> & b, T t) { return a*(1-t) + b*t; }
    template<class T, int M> vec<T,M>             nlerp    (const vec<T,M> & a, const vec<T,M> & b, T t) { return normalize(lerp(a,b,t)); }
    template<class T, int M> vec<T,M>             slerp    (const vec<T,M> & a, const vec<T,M> & b, T t) { T th=uangle(a,b); return th == 0 ? a : a*(std::sin(th*(1-t))/std::sin(th)) + b*(std::sin(th*t)/std::sin(th)); }
    template<class T, int M> constexpr mat<T,M,2> outerprod(const vec<T,M> & a, const vec<T,2> & b)      { return {a*b.x, a*b.y}; }
    template<class T, int M> constexpr mat<T,M,3> outerprod(const vec<T,M> & a, const vec<T,3> & b)      { return {a*b.x, a*b.y, a*b.z}; }
    template<class T, int M> constexpr mat<T,M,4> outerprod(const vec<T,M> & a, const vec<T,4> & b)      { return {a*b.x, a*b.y, a*b.z, a*b.w}; }

    // Support for quaternion algebra using 4D vectors, representing xi + yj + zk + w
    template<class T> constexpr vec<T,4> qconj(const vec<T,4> & q)                     { return {-q.x,-q.y,-q.z,q.w}; }
    template<class T> vec<T,4>           qinv (const vec<T,4> & q)                     { return qconj(q)/length2(q); }
    template<class T> vec<T,4>           qexp (const vec<T,4> & q)                     { const auto v = q.xyz(); const auto vv = length(v); return std::exp(q.w) * vec<T,4>{v * (vv > 0 ? std::sin(vv)/vv : 0), std::cos(vv)}; }
    template<class T> vec<T,4>           qlog (const vec<T,4> & q)                     { const auto v = q.xyz(); const auto vv = length(v), qq = length(q); return {v * (vv > 0 ? std::acos(q.w/qq)/vv : 0), std::log(qq)}; }
    template<class T> vec<T,4>           qpow (const vec<T,4> & q, const T & p)        { const auto v = q.xyz(); const auto vv = length(v), qq = length(q), th = std::acos(q.w/qq); return std::pow(qq,p)*vec<T,4>{v * (vv > 0 ? std::sin(p*th)/vv : 0), std::cos(p*th)}; }
    template<class T> constexpr vec<T,4> qmul (const vec<T,4> & a, const vec<T,4> & b) { return {a.x*b.w+a.w*b.x+a.y*b.z-a.z*b.y, a.y*b.w+a.w*b.y+a.z*b.x-a.x*b.z, a.z*b.w+a.w*b.z+a.x*b.y-a.y*b.x, a.w*b.w-a.x*b.x-a.y*b.y-a.z*b.z}; }
    template<class T, class... R> constexpr vec<T,4> qmul(const vec<T,4> & a, R... r)  { return qmul(a, qmul(r...)); }

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

    // Support for matrix algebra
    template<class T, int M> constexpr vec<T,M> mul(const mat<T,M,2> & a, const vec<T,2> & b) { return a.x*b.x + a.y*b.y; }
    template<class T, int M> constexpr vec<T,M> mul(const mat<T,M,3> & a, const vec<T,3> & b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
    template<class T, int M> constexpr vec<T,M> mul(const mat<T,M,4> & a, const vec<T,4> & b) { return a.x*b.x + a.y*b.y + a.z*b.z + a.w*b.w; }
    template<class T, int M, int N> constexpr mat<T,M,2> mul(const mat<T,M,N> & a, const mat<T,N,2> & b) { return {mul(a,b.x), mul(a,b.y)}; }
    template<class T, int M, int N> constexpr mat<T,M,3> mul(const mat<T,M,N> & a, const mat<T,N,3> & b) { return {mul(a,b.x), mul(a,b.y), mul(a,b.z)}; }
    template<class T, int M, int N> constexpr mat<T,M,4> mul(const mat<T,M,N> & a, const mat<T,N,4> & b) { return {mul(a,b.x), mul(a,b.y), mul(a,b.z), mul(a,b.w)}; }
    #if _MSC_VER >= 1910
    template<class T, int M, int N, class... R> constexpr auto mul(const mat<T,M,N> & a, R... r) { return mul(a, mul(r...)); }
    #else
    template<class T, int M, int N, class... R> constexpr auto mul(const mat<T,M,N> & a, R... r) -> decltype(mul(a, mul(r...))) { return mul(a, mul(r...)); }
    #endif    
    template<class T> constexpr vec<T,2> diagonal(const mat<T,2,2> & a) { return {a.x.x, a.y.y}; }
    template<class T> constexpr vec<T,3> diagonal(const mat<T,3,3> & a) { return {a.x.x, a.y.y, a.z.z}; }
    template<class T> constexpr vec<T,4> diagonal(const mat<T,4,4> & a) { return {a.x.x, a.y.y, a.z.z, a.w.w}; }
    template<class T, int M> constexpr mat<T,M,2> transpose(const mat<T,2,M> & m) { return {m.row(0), m.row(1)}; }
    template<class T, int M> constexpr mat<T,M,3> transpose(const mat<T,3,M> & m) { return {m.row(0), m.row(1), m.row(2)}; }
    template<class T, int M> constexpr mat<T,M,4> transpose(const mat<T,4,M> & m) { return {m.row(0), m.row(1), m.row(2), m.row(3)}; }
    template<class T> constexpr mat<T,2,2> adjugate(const mat<T,2,2> & a) { return {{a.y.y, -a.x.y}, {-a.y.x, a.x.x}}; }
    template<class T> constexpr mat<T,3,3> adjugate(const mat<T,3,3> & a);
    template<class T> constexpr mat<T,4,4> adjugate(const mat<T,4,4> & a);
    template<class T> constexpr T determinant(const mat<T,2,2> & a) { return a.x.x*a.y.y - a.x.y*a.y.x; }
    template<class T> constexpr T determinant(const mat<T,3,3> & a) { return a.x.x*(a.y.y*a.z.z - a.z.y*a.y.z) + a.x.y*(a.y.z*a.z.x - a.z.z*a.y.x) + a.x.z*(a.y.x*a.z.y - a.z.x*a.y.y); }
    template<class T> constexpr T determinant(const mat<T,4,4> & a);
    template<class T, int N> constexpr mat<T,N,N> inverse(const mat<T,N,N> & a) { return adjugate(a)/determinant(a); }

    // Vectors and matrices can be used as ranges
    template<class T, int M>       T * begin(      vec<T,M> & a) { return &a[0]; }
    template<class T, int M> const T * begin(const vec<T,M> & a) { return &a[0]; }
    template<class T, int M>       T * end  (      vec<T,M> & a) { return begin(a) + M; }
    template<class T, int M> const T * end  (const vec<T,M> & a) { return begin(a) + M; }
    template<class T, int M, int N>       vec<T,M> * begin(      mat<T,M,N> & a) { return &a[0]; }
    template<class T, int M, int N> const vec<T,M> * begin(const mat<T,M,N> & a) { return &a[0]; }
    template<class T, int M, int N>       vec<T,M> * end  (      mat<T,M,N> & a) { return begin(a) + N; }
    template<class T, int M, int N> const vec<T,M> * end  (const mat<T,M,N> & a) { return begin(a) + N; }

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
        typedef mat<bool,2,2> bool2x2; typedef mat<int,2,2> int2x2; typedef mat<float,2,2> float2x2; typedef mat<double,2,2> double2x2;
        typedef mat<bool,2,3> bool2x3; typedef mat<int,2,3> int2x3; typedef mat<float,2,3> float2x3; typedef mat<double,2,3> double2x3;
        typedef mat<bool,2,4> bool2x4; typedef mat<int,2,4> int2x4; typedef mat<float,2,4> float2x4; typedef mat<double,2,4> double2x4;
        typedef mat<bool,3,2> bool3x2; typedef mat<int,3,2> int3x2; typedef mat<float,3,2> float3x2; typedef mat<double,3,2> double3x2;
        typedef mat<bool,3,3> bool3x3; typedef mat<int,3,3> int3x3; typedef mat<float,3,3> float3x3; typedef mat<double,3,3> double3x3;
        typedef mat<bool,3,4> bool3x4; typedef mat<int,3,4> int3x4; typedef mat<float,3,4> float3x4; typedef mat<double,3,4> double3x4;
        typedef mat<bool,4,2> bool4x2; typedef mat<int,4,2> int4x2; typedef mat<float,4,2> float4x2; typedef mat<double,4,2> double4x2;
        typedef mat<bool,4,3> bool4x3; typedef mat<int,4,3> int4x3; typedef mat<float,4,3> float4x3; typedef mat<double,4,3> double4x3;
        typedef mat<bool,4,4> bool4x4; typedef mat<int,4,4> int4x4; typedef mat<float,4,4> float4x4; typedef mat<double,4,4> double4x4;
    }
}

// Definitions of functions too long to be defined inline
template<class T> constexpr linalg::mat<T,3,3> linalg::adjugate(const mat<T,3,3> & a) 
{ 
    return {{a.y.y*a.z.z - a.z.y*a.y.z, a.z.y*a.x.z - a.x.y*a.z.z, a.x.y*a.y.z - a.y.y*a.x.z},
            {a.y.z*a.z.x - a.z.z*a.y.x, a.z.z*a.x.x - a.x.z*a.z.x, a.x.z*a.y.x - a.y.z*a.x.x},
            {a.y.x*a.z.y - a.z.x*a.y.y, a.z.x*a.x.y - a.x.x*a.z.y, a.x.x*a.y.y - a.y.x*a.x.y}}; 
}

template<class T> constexpr linalg::mat<T,4,4> linalg::adjugate(const mat<T,4,4> & a) 
{ 
    return {{a.y.y*a.z.z*a.w.w + a.w.y*a.y.z*a.z.w + a.z.y*a.w.z*a.y.w - a.y.y*a.w.z*a.z.w - a.z.y*a.y.z*a.w.w - a.w.y*a.z.z*a.y.w,
             a.x.y*a.w.z*a.z.w + a.z.y*a.x.z*a.w.w + a.w.y*a.z.z*a.x.w - a.w.y*a.x.z*a.z.w - a.z.y*a.w.z*a.x.w - a.x.y*a.z.z*a.w.w,
             a.x.y*a.y.z*a.w.w + a.w.y*a.x.z*a.y.w + a.y.y*a.w.z*a.x.w - a.x.y*a.w.z*a.y.w - a.y.y*a.x.z*a.w.w - a.w.y*a.y.z*a.x.w,
             a.x.y*a.z.z*a.y.w + a.y.y*a.x.z*a.z.w + a.z.y*a.y.z*a.x.w - a.x.y*a.y.z*a.z.w - a.z.y*a.x.z*a.y.w - a.y.y*a.z.z*a.x.w},
            {a.y.z*a.w.w*a.z.x + a.z.z*a.y.w*a.w.x + a.w.z*a.z.w*a.y.x - a.y.z*a.z.w*a.w.x - a.w.z*a.y.w*a.z.x - a.z.z*a.w.w*a.y.x,
             a.x.z*a.z.w*a.w.x + a.w.z*a.x.w*a.z.x + a.z.z*a.w.w*a.x.x - a.x.z*a.w.w*a.z.x - a.z.z*a.x.w*a.w.x - a.w.z*a.z.w*a.x.x,
             a.x.z*a.w.w*a.y.x + a.y.z*a.x.w*a.w.x + a.w.z*a.y.w*a.x.x - a.x.z*a.y.w*a.w.x - a.w.z*a.x.w*a.y.x - a.y.z*a.w.w*a.x.x,
             a.x.z*a.y.w*a.z.x + a.z.z*a.x.w*a.y.x + a.y.z*a.z.w*a.x.x - a.x.z*a.z.w*a.y.x - a.y.z*a.x.w*a.z.x - a.z.z*a.y.w*a.x.x},
            {a.y.w*a.z.x*a.w.y + a.w.w*a.y.x*a.z.y + a.z.w*a.w.x*a.y.y - a.y.w*a.w.x*a.z.y - a.z.w*a.y.x*a.w.y - a.w.w*a.z.x*a.y.y,
             a.x.w*a.w.x*a.z.y + a.z.w*a.x.x*a.w.y + a.w.w*a.z.x*a.x.y - a.x.w*a.z.x*a.w.y - a.w.w*a.x.x*a.z.y - a.z.w*a.w.x*a.x.y,
             a.x.w*a.y.x*a.w.y + a.w.w*a.x.x*a.y.y + a.y.w*a.w.x*a.x.y - a.x.w*a.w.x*a.y.y - a.y.w*a.x.x*a.w.y - a.w.w*a.y.x*a.x.y,
             a.x.w*a.z.x*a.y.y + a.y.w*a.x.x*a.z.y + a.z.w*a.y.x*a.x.y - a.x.w*a.y.x*a.z.y - a.z.w*a.x.x*a.y.y - a.y.w*a.z.x*a.x.y},
            {a.y.x*a.w.y*a.z.z + a.z.x*a.y.y*a.w.z + a.w.x*a.z.y*a.y.z - a.y.x*a.z.y*a.w.z - a.w.x*a.y.y*a.z.z - a.z.x*a.w.y*a.y.z,
             a.x.x*a.z.y*a.w.z + a.w.x*a.x.y*a.z.z + a.z.x*a.w.y*a.x.z - a.x.x*a.w.y*a.z.z - a.z.x*a.x.y*a.w.z - a.w.x*a.z.y*a.x.z,
             a.x.x*a.w.y*a.y.z + a.y.x*a.x.y*a.w.z + a.w.x*a.y.y*a.x.z - a.x.x*a.y.y*a.w.z - a.w.x*a.x.y*a.y.z - a.y.x*a.w.y*a.x.z,
             a.x.x*a.y.y*a.z.z + a.z.x*a.x.y*a.y.z + a.y.x*a.z.y*a.x.z - a.x.x*a.z.y*a.y.z - a.y.x*a.x.y*a.z.z - a.z.x*a.y.y*a.x.z}}; 
}

template<class T> constexpr T linalg::determinant(const mat<T,4,4> & a) 
{ 
    return a.x.x*(a.y.y*a.z.z*a.w.w + a.w.y*a.y.z*a.z.w + a.z.y*a.w.z*a.y.w - a.y.y*a.w.z*a.z.w - a.z.y*a.y.z*a.w.w - a.w.y*a.z.z*a.y.w)
         + a.x.y*(a.y.z*a.w.w*a.z.x + a.z.z*a.y.w*a.w.x + a.w.z*a.z.w*a.y.x - a.y.z*a.z.w*a.w.x - a.w.z*a.y.w*a.z.x - a.z.z*a.w.w*a.y.x)
         + a.x.z*(a.y.w*a.z.x*a.w.y + a.w.w*a.y.x*a.z.y + a.z.w*a.w.x*a.y.y - a.y.w*a.w.x*a.z.y - a.z.w*a.y.x*a.w.y - a.w.w*a.z.x*a.y.y)
         + a.x.w*(a.y.x*a.w.y*a.z.z + a.z.x*a.y.y*a.w.z + a.w.x*a.z.y*a.y.z - a.y.x*a.z.y*a.w.z - a.w.x*a.y.y*a.z.z - a.z.x*a.w.y*a.y.z); 
}

template<class T> linalg::vec<T,4> linalg::rotation_quat(const mat<T,3,3> & m)
{
    const vec<T,4> q {m.x.x-m.y.y-m.z.z, m.y.y-m.x.x-m.z.z, m.z.z-m.x.x-m.y.y, m.x.x+m.y.y+m.z.z}, s[] {
        {1, m.x.y + m.y.x, m.z.x + m.x.z, m.y.z - m.z.y}, 
        {m.x.y + m.y.x, 1, m.y.z + m.z.y, m.z.x - m.x.z},
        {m.x.z + m.z.x, m.y.z + m.z.y, 1, m.x.y - m.y.x},
        {m.y.z - m.z.y, m.z.x - m.x.z, m.x.y - m.y.x, 1}};
    return copysign(normalize(sqrt(max(T(0), T(1)+q))), s[argmax(q)]);
}

template<class T> linalg::mat<T,4,4> linalg::frustum_matrix(T x0, T x1, T y0, T y1, T n, T f, fwd_axis a, z_range z) 
{ 
    const T s = a == pos_z ? T(1) : T(-1), o = z == neg_one_to_one ? n : 0;
    return {{2*n/(x1-x0),0,0,0}, {0,2*n/(y1-y0),0,0}, {-s*(x0+x1)/(x1-x0),-s*(y0+y1)/(y1-y0),s*(f+o)/(f-n),s}, {0,0,-(n+o)*f/(f-n),0}};
}

#endif
