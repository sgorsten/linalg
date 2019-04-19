// linalg-std.h - Addendum to linalg.h providing namespace std interop
//
// This file provides additional functionality which allows types defined
// in namespace linalg to interoperate with functionality from namespace std.
//
// The original author of this software is Sterling Orsten, and its permanent
// home is <http://github.com/sgorsten/linalg/>. If you find this software
// useful, an acknowledgement in your source text and/or product documentation
// is appreciated, but not required.



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

#include "linalg.h"
#include <array>
#include <iosfwd>

namespace linalg
{
    // Provide implicit conversion between linalg::vec<T,M> and std::array<T,M>
    template<class T> struct converter<vec<T,2>, std::array<T,2>> { vec<T,2> operator() (const std::array<T,2> & a) const { return {a[0], a[1]}; } };
    template<class T> struct converter<vec<T,3>, std::array<T,3>> { vec<T,3> operator() (const std::array<T,3> & a) const { return {a[0], a[1], a[2]}; } };
    template<class T> struct converter<vec<T,4>, std::array<T,4>> { vec<T,4> operator() (const std::array<T,4> & a) const { return {a[0], a[1], a[2], a[3]}; } };

    template<class T> struct converter<std::array<T,2>, vec<T,2>> { std::array<T,2> operator() (const vec<T,2> & a) const { return {a[0], a[1]}; } };
    template<class T> struct converter<std::array<T,3>, vec<T,3>> { std::array<T,3> operator() (const vec<T,3> & a) const { return {a[0], a[1], a[2]}; } };
    template<class T> struct converter<std::array<T,4>, vec<T,4>> { std::array<T,4> operator() (const vec<T,4> & a) const { return {a[0], a[1], a[2], a[3]}; } };

    // Provide output streaming operators, writing something that resembles an aggregate literal that could be used to construct the specified value
    namespace ostream_overloads
    {
        template<class C, class T> std::basic_ostream<C> & operator << (std::basic_ostream<C> & out, const vec<T,2> & v) { return out << '{' << v[0] << ',' << v[1] << '}'; }
        template<class C, class T> std::basic_ostream<C> & operator << (std::basic_ostream<C> & out, const vec<T,3> & v) { return out << '{' << v[0] << ',' << v[1] << ',' << v[2] << '}'; }
        template<class C, class T> std::basic_ostream<C> & operator << (std::basic_ostream<C> & out, const vec<T,4> & v) { return out << '{' << v[0] << ',' << v[1] << ',' << v[2] << ',' << v[3] << '}'; }

        template<class C, class T, int M> std::basic_ostream<C> & operator << (std::basic_ostream<C> & out, const mat<T,M,2> & m) { return out << '{' << m[0] << ',' << m[1] << '}'; }
        template<class C, class T, int M> std::basic_ostream<C> & operator << (std::basic_ostream<C> & out, const mat<T,M,3> & m) { return out << '{' << m[0] << ',' << m[1] << ',' << m[2] << '}'; }
        template<class C, class T, int M> std::basic_ostream<C> & operator << (std::basic_ostream<C> & out, const mat<T,M,4> & m) { return out << '{' << m[0] << ',' << m[1] << ',' << m[2] << ',' << m[3] << '}'; }
    }
}

namespace std 
{ 
    // Provide specializations for std::hash<...> with linalg types
    template<class T> struct hash<linalg::vec<T,2>> { std::size_t operator()(const linalg::vec<T,2> & v) const { std::hash<T> h; return h(v.x) ^ (h(v.y) << 1); } };
    template<class T> struct hash<linalg::vec<T,3>> { std::size_t operator()(const linalg::vec<T,3> & v) const { std::hash<T> h; return h(v.x) ^ (h(v.y) << 1) ^ (h(v.z) << 2); } };
    template<class T> struct hash<linalg::vec<T,4>> { std::size_t operator()(const linalg::vec<T,4> & v) const { std::hash<T> h; return h(v.x) ^ (h(v.y) << 1) ^ (h(v.z) << 2) ^ (h(v.w) << 3); } };

    template<class T, int M> struct hash<linalg::mat<T,M,2>> { std::size_t operator()(const linalg::mat<T,M,2> & v) const { std::hash<linalg::vec<T,M>> h; return h(v.x) ^ (h(v.y) << M); } };
    template<class T, int M> struct hash<linalg::mat<T,M,3>> { std::size_t operator()(const linalg::mat<T,M,3> & v) const { std::hash<linalg::vec<T,M>> h; return h(v.x) ^ (h(v.y) << M) ^ (h(v.z) << (M*2)); } };
    template<class T, int M> struct hash<linalg::mat<T,M,4>> { std::size_t operator()(const linalg::mat<T,M,4> & v) const { std::hash<linalg::vec<T,M>> h; return h(v.x) ^ (h(v.y) << M) ^ (h(v.z) << (M*2)) ^ (h(v.w) << (M*3)); } };
}