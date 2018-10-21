#pragma once
#ifndef LINALGX_H
#define LINALGX_H

// linalgx.h - Unstable extensions to linalg.h
//
// The intent of this file is to publish frequently used extensions to
// linalg.h, without offering a commitment to backwards compatibility.
//
// The original author of this software is Sterling Orsten, and its permanent
// home is <http://github.com/sgorsten/linalg/>. If you find this software
// useful, an acknowledgement in your source text and/or product documentation
// is appreciated, but not required.
//
// The author acknowledges significant insights and contributions by:
//     Dimitri Diakopoulos <http://github.com/ddiakopoulos/>
//     Stan Melax <http://github.com/melax/>



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
#include <cstdint>    // For std::uint8_t, std::uint16_t, std::int16_t, etc.
#include <functional> // For std::hash

namespace linalg
{
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

    namespace aliases
    {
         using byte1=vec<uint8_t,1>; using short1=vec<int16_t,1>; using ushort1=vec<uint16_t,1>; using uint1=vec<unsigned,1>;  
         using byte2=vec<uint8_t,2>; using short2=vec<int16_t,2>; using ushort2=vec<uint16_t,2>; using uint2=vec<unsigned,2>;  
         using byte3=vec<uint8_t,3>; using short3=vec<int16_t,3>; using ushort3=vec<uint16_t,3>; using uint3=vec<unsigned,3>;  
         using byte4=vec<uint8_t,4>; using short4=vec<int16_t,4>; using ushort4=vec<uint16_t,4>; using uint4=vec<unsigned,4>;      
    }
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

#endif