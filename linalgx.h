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
    // A value type representing an abstract direction vector in 3D space, independent of any coordinate system
    enum class coord_axis { forward, back, left, right, up, down };
    template<class T> constexpr T dot(coord_axis a, coord_axis b) { return a == b ? T(1) : (static_cast<int>(a) ^ static_cast<int>(b)) == 1 ? T(-1) : T(0); }

    // A concrete 3D coordinate system with defined x, y, and z axes
    struct coord_system
    {
        coord_axis x_axis, y_axis, z_axis;
        constexpr coord_system(coord_axis x, coord_axis y, coord_axis z) : x_axis{x}, y_axis{y}, z_axis{z} {}
        template<class T> constexpr vec<T,3> get(coord_axis axis) const { return {dot<T>(x_axis, axis), dot<T>(y_axis, axis), dot<T>(z_axis, axis)}; }
        template<class T> constexpr vec<T,3> cross(coord_axis a, coord_axis b) const { return linalg::cross(get<T>(a), get<T>(b)); }
        constexpr bool is_orthogonal() const { return dot<int>(x_axis, y_axis) == 0 && dot<int>(y_axis, z_axis) == 0 && dot<int>(z_axis, x_axis) == 0; }
        constexpr bool is_left_handed() const { return dot(cross<int>(coord_axis::forward, coord_axis::up), get<int>(coord_axis::left)) == 1; }
        constexpr bool is_right_handed() const { return dot(cross<int>(coord_axis::forward, coord_axis::up), get<int>(coord_axis::right)) == 1; }
    };

    template<class T> vec<T,2> rot(T a, const vec<T,2> & v) { const T s = std::sin(a), c = std::cos(a); return {v[0]*c - v[1]*s, v[0]*s + v[1]*c}; }

    // Produce a rotation quaternion representing a rotation of `angle` radians around the specified unit-length `axis`
    template<class T> quat<T> rotation_quat(const vec<T,3> & axis, T angle) { return {axis*std::sin(angle/2), std::cos(angle/2)}; }

    // Produce the rotation quaternion of shortest arc which rotates vector `from` to be parallel with vector `to`
    template<class T> quat<T> rotation_quat(const vec<T,3> & from, const vec<T,3> & to) { return rotation_quat(normalize(cross(from,to)), angle(from,to)); }

    // Produce the rotation quaternion which stores the equivalent rotation to the rotation matrix `rotation_matrix`
    template<class T> quat<T> rotation_quat(const mat<T,3,3> & rotation_matrix)
    {
        const auto & m = rotation_matrix;
        const vec<T,4> q {m[0][0]-m[1][1]-m[2][2], m[1][1]-m[0][0]-m[2][2], m[2][2]-m[0][0]-m[1][1], m[0][0]+m[1][1]+m[2][2]}, s[] {
            {1, m[0][1] + m[1][0], m[2][0] + m[0][2], m[1][2] - m[2][1]}, 
            {m[0][1] + m[1][0], 1, m[1][2] + m[2][1], m[2][0] - m[0][2]},
            {m[0][2] + m[2][0], m[1][2] + m[2][1], 1, m[0][1] - m[1][0]},
            {m[1][2] - m[2][1], m[2][0] - m[0][2], m[0][1] - m[1][0], 1}};
        return quat<T>{copysign(normalize(sqrt(max(T(0), T(1)+q))), s[argmax(q)])};    
    }

    // Transformation types
    enum fwd_axis { neg_z, pos_z }; // Should projection matrices be generated assuming forward is {0,0,-1} or {0,0,1}
    enum z_range { neg_one_to_one, zero_to_one }; // Should projection matrices map z into the range of [-1,1] or [0,1]?
    template<class T, int M> struct linear_transformation;
    template<class T, int M> struct homogeneous_transformation;

    // Helper to construct linear transformations in 3D space
    template<class T> struct linear_transformation<T,3>
    {
        using type = mat<T,3,3>;

        constexpr static type identity() { return linalg::identity; }
        constexpr static type coord_change(const coord_system & from, const coord_system & to)  { return {to.get<T>(from.x_axis), to.get<T>(from.y_axis), to.get<T>(from.z_axis)}; }    
        constexpr static type scaling(T uniform_scale)                                          { return {{uniform_scale,0,0},{0,uniform_scale,0},{0,0,uniform_scale}}; }
        constexpr static type scaling(const vec<T,3> & scale)                                   { return {{scale[0],0,0},{0,scale[1],0},{0,0,scale[2]}}; }
        constexpr static type rotation(const quat<T> & rotation)                                { return {qxdir(rotation), qydir(rotation), qzdir(rotation)}; }
    };

    // Helper to construct homogeneous transformations in 3D space
    template<class T> struct homogeneous_transformation<T,3>
    {
        using type = mat<T,4,4>;

        constexpr static type identity() { return linalg::identity; }
        constexpr static type linear(const mat<T,3,3> & transform)                              { return {{transform[0],0}, {transform[1],0}, {transform[2],0}, {0,0,0,1}}; }
        constexpr static type coord_change(const coord_system & from, const coord_system & to)  { return {{to.get<T>(from.x_axis),0}, {to.get<T>(from.y_axis),0}, {to.get<T>(from.z_axis),0}, {0,0,0,1}}; }    
        constexpr static type scaling(T uniform_scale)                                          { return {{uniform_scale,0,0,0},{0,uniform_scale,0,0},{0,0,uniform_scale,0},{0,0,0,1}}; }
        constexpr static type scaling(const vec<T,3> & scale)                                   { return {{scale[0],0,0,0},{0,scale[1],0,0},{0,0,scale[2],0},{0,0,0,1}}; }
        constexpr static type rotation(const quat<T> & rotation)                                { return {{qxdir(rotation),0}, {qydir(rotation),0}, {qzdir(rotation),0}, {0,0,0,1}}; }
        constexpr static type translation(const vec<T,3> & translation)                         { return {{1,0,0,0}, {0,1,0,0}, {0,0,1,0}, {translation,1}}; }
        constexpr static type pose(const quat<T> & rotation, const vec<T,3> & translation)      { return {{qxdir(rotation),0}, {qydir(rotation),0}, {qzdir(rotation),0}, {translation,1}}; }

        constexpr static type frustum(T x0, T x1, T y0, T y1, T n, T f, fwd_axis a = neg_z, z_range z = neg_one_to_one) 
        { 
            return (z == neg_one_to_one ? translation({0,0,-1}) * scaling({1,1,2}) : identity()) 
                * type{{2*n/(x1-x0),0,0,0}, {0,2*n/(y1-y0),0,0}, {-(x0+x1)/(x1-x0), -(y0+y1)/(y1-y0), f/(f-n), 1}, {0,0,-n*f/(f-n),0}}
                * (a == neg_z ? scaling({1,1,-1}) : identity()); 
        }
        static type perspective(T fovy, T aspect, T n, T f, fwd_axis a = neg_z, z_range z = neg_one_to_one)
        { 
            T y = n*std::tan(fovy / 2), x = y*aspect;
            return frustum(-x, x, -y, y, n, f, a, z); 
        }
    };

    template<class T> constexpr vec<T,3> project(const vec<T,4> & p) { return swizzle<0,1,2>(p) / p[3]; }

    // Transformation functions
    template<class T> constexpr vec<T,3> transform_vector   (const mat<T,3,3> & linear_transform, const vec<T,3> & vector   ) { return linear_transform * vector; }
    template<class T> constexpr vec<T,3> transform_direction(const mat<T,3,3> & linear_transform, const vec<T,3> & direction) { return normalize(transform_vector(linear_transform, direction)); }
    template<class T> constexpr vec<T,3> transform_bivector (const mat<T,3,3> & linear_transform, const vec<T,3> & bivector ) { return transform_vector(comatrix(linear_transform), bivector); }
    template<class T> constexpr vec<T,3> transform_normal   (const mat<T,3,3> & linear_transform, const vec<T,3> & normal   ) { return normalize(transform_bivector(linear_transform, normal)); }
    template<class T> constexpr vec<T,3> transform_point    (const mat<T,3,3> & linear_transform, const vec<T,3> & point    ) { return transform_vector(linear_transform, point); }

    template<class T> constexpr vec<T,3> transform_vector   (const mat<T,4,4> & homogeneous_transform, const vec<T,3> & vector   ) { return swizzle<0,1,2>(homogeneous_transform * vec<T,4>{vector,0}); }
    template<class T> constexpr vec<T,3> transform_direction(const mat<T,4,4> & homogeneous_transform, const vec<T,3> & direction) { return normalize(transform_vector(homogeneous_transform, direction)); }
    template<class T> constexpr vec<T,3> transform_bivector (const mat<T,4,4> & homogeneous_transform, const vec<T,3> & bivector ) { return transform_vector(comatrix(homogeneous_transform), bivector); }
    template<class T> constexpr vec<T,3> transform_normal   (const mat<T,4,4> & homogeneous_transform, const vec<T,3> & normal   ) { return normalize(transform_bivector(homogeneous_transform, normal)); }
    template<class T> constexpr vec<T,3> transform_point    (const mat<T,4,4> & homogeneous_transform, const vec<T,3> & point    ) { return project(homogeneous_transform * vec<T,4>{point,1}); }

    namespace aliases
    {
         using byte1=vec<uint8_t,1>; using short1=vec<int16_t,1>; using ushort1=vec<uint16_t,1>; using uint1=vec<unsigned,1>;  
         using byte2=vec<uint8_t,2>; using short2=vec<int16_t,2>; using ushort2=vec<uint16_t,2>; using uint2=vec<unsigned,2>;  
         using byte3=vec<uint8_t,3>; using short3=vec<int16_t,3>; using ushort3=vec<uint16_t,3>; using uint3=vec<unsigned,3>;  
         using byte4=vec<uint8_t,4>; using short4=vec<int16_t,4>; using ushort4=vec<uint16_t,4>; using uint4=vec<unsigned,4>;      
    }
}

////////////////////////////////////////////////////////
// Specializations of std::hash<...> for linalg types //
////////////////////////////////////////////////////////

namespace std 
{ 
    template<class T> struct hash<linalg::vec<T,1>> { std::size_t operator()(const linalg::vec<T,1> & v) const { std::hash<T> h; return h(v[0]); } };
    template<class T> struct hash<linalg::vec<T,2>> { std::size_t operator()(const linalg::vec<T,2> & v) const { std::hash<T> h; return h(v[0]) ^ (h(v[1]) << 1); } };
    template<class T> struct hash<linalg::vec<T,3>> { std::size_t operator()(const linalg::vec<T,3> & v) const { std::hash<T> h; return h(v[0]) ^ (h(v[1]) << 1) ^ (h(v[2]) << 2); } };
    template<class T> struct hash<linalg::vec<T,4>> { std::size_t operator()(const linalg::vec<T,4> & v) const { std::hash<T> h; return h(v[0]) ^ (h(v[1]) << 1) ^ (h(v[2]) << 2) ^ (h(v[3]) << 3); } };
    template<class T, int M> struct hash<linalg::mat<T,M,1>> { std::size_t operator()(const linalg::mat<T,M,1> & m) const { std::hash<linalg::vec<T,M>> h; return h(m[0]); } };
    template<class T, int M> struct hash<linalg::mat<T,M,2>> { std::size_t operator()(const linalg::mat<T,M,2> & m) const { std::hash<linalg::vec<T,M>> h; return h(m[0]) ^ (h(m[1]) << M); } };
    template<class T, int M> struct hash<linalg::mat<T,M,3>> { std::size_t operator()(const linalg::mat<T,M,3> & m) const { std::hash<linalg::vec<T,M>> h; return h(m[0]) ^ (h(m[1]) << M) ^ (h(m[2]) << (M*2)); } };
    template<class T, int M> struct hash<linalg::mat<T,M,4>> { std::size_t operator()(const linalg::mat<T,M,4> & m) const { std::hash<linalg::vec<T,M>> h; return h(m[0]) ^ (h(m[1]) << M) ^ (h(m[2]) << (M*2)) ^ (h(m[3]) << (M*3)); } };
    template<class T> struct hash<linalg::quat<T>> { std::size_t operator()(const linalg::quat<T> & q) const { std::hash<T> h; return h(q.x) ^ (h(q.y) << 1) ^ (h(q.z) << 2) ^ (h(q.w) << 3); } };
}

#endif