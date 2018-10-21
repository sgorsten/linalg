**NOTE: This is the development branch for `linalg.h v3.0`. Breaking changes may occur until the `v3.0` tag is published.**

[![Release is 3.0-beta](https://img.shields.io/badge/version-3.0--beta-blue.svg)](http://raw.githubusercontent.com/sgorsten/linalg/v3/linalg.h)
[![License is Unlicense](http://img.shields.io/badge/license-Unlicense-blue.svg?style=flat)](http://unlicense.org/)
[![Travis CI build status](http://travis-ci.org/sgorsten/linalg.svg?branch=v3)](https://travis-ci.org/sgorsten/linalg)
[![Appveyor build status](http://ci.appveyor.com/api/projects/status/l4bfv5omodkajuc9?svg=true)](https://ci.appveyor.com/project/sgorsten/linalg)

[`linalg.h`](/linalg.h) is a [single header](http://github.com/nothings/stb/blob/master/docs/other_libs.md), [public domain](http://unlicense.org/) linear algebra library for [C++11](http://en.cppreference.com/w/). It is inspired by the syntax of popular shading and compute languages and is intended to serve as a lightweight alternative to projects such as [GLM](http://glm.g-truc.net/0.9.7/), [Boost.QVM](https://www.boost.org/doc/libs/1_66_0/libs/qvm/doc/index.html) or [Eigen](http://eigen.tuxfamily.org/) in the domains such as computer graphics, computational geometry, and physical simulation. It allows you to easily write programs like the following:

```cpp
#include <linalg.h>
using namespace linalg::aliases;

// Compute the equation of a plane containing points a, b, and c
float4 compute_plane(float3 a, float3 b, float3 c)
{
    float3 n = cross(b-a, c-a);
    return {n, -dot(n,a)};
}
```

* [Data structures](#data-structures)
  * [Vectors](#vectors)
  * [Matrices](#matrices)
  * [Quaternions](#quaternions)
* [Function listing](#function-listing)
  * [Vector algebra](#vector-algebra)
  * [Matrix algebra](#matrix-algebra)
  * [Quaternion algebra](#quaternion-algebra)
  * [Component-wise operations](#component-wise-operations)
  * [Reductions](#reductions)
* [Optional features](#optional-features)
  * [Type aliases](#type-aliases)
  * [`ostream` overloads](#ostream-overloads)
  * [User-defined conversions](#user-defined-conversions)
  * [Extensions in `linalgx.h`](#extensions-in-linalgxh)
* [Higher order functions](#higher-order-functions)
* [Design rationale](#design-rationale)
* [Changes from v2.1](#changes-from-v21)

## Data structures

#### Vectors

`linalg::vec<T,M>` defines a fixed-length vector containing exactly `M` elements of type `T`. Convenience aliases such as `float3`, `float4`, or `int2` are provided in the [`linalg::aliases` namespace](#type-aliases). This data structure can be used to store a wide variety of types of data, including geometric vectors, points, homogeneous coordinates, plane equations, colors, texture coordinates, or any other situation where you need to manipulate a small sequence of numbers. As such, `vec<T,M>` is supported by a set of [algebraic](#vector-algebra) and [component-wise](#component-wise-operations) functions, as well as a set of standard [reductions](#reductions).

`vec<T,M>`:
* is default-constructible: `float3 v; // v contains 0,0,0`
* is implicitly constructible from `M` elements of type `T`: `float3 v {1,2,3}; // v contains 1,2,3`
* is explicitly constructible from a single element of type `T`: `float3 v {4}; // v contains 4,4,4`
* is explicitly constructible from a `vec<T,U>` of some other type `U`: `float3 v {int3{5,6,7}}; // v contains 5,6,7`
* is copy-constructible: `float3 v {1,2,3}, u {v}; // u and v contain 1,2,3`
* is copy-assignable: `float3 v {4,5,6}, u; u = v; // u and v contain 4,5,6`
* supports indexing: `float x = v[0]; // x contains first element of v`
* supports named accessors `x`,`y`,`z`,`w`: `float y = point.y; // y contains second element of point`
* supports named accessors `r`,`g`,`b`,`a`: `pixel.a = 0.5; // fourth element of pixel set to 0.5` 
* supports named accessors `s`,`t`,`p`,`q`: `float s = tc.s; // s contains first element of tc`
* supports swizzles: `float3 c = pixel.bgr; // c contains pixel[2],pixel[1],pixel[0]
* is [`EqualityComparable`](http://en.cppreference.com/w/cpp/concept/EqualityComparable): `bool b = (v == u); // b is true if v and u contain equal elements`
* is [`LessThanComparable`](http://en.cppreference.com/w/cpp/concept/LessThanComparable): `bool b = (v < u); // b is true if v precedes u lexicographically`
* supports unary operators `+`, `-`, `!` and `~` in component-wise fashion: `auto v = -float{2,3}; // v is float2{-2,-3}`
* supports binary operators `+`, `-`, `*`, `/`, `%`, `|`, `&`, `^`, `<<` and `>>` in component-wise fashion: `auto v = float2{1,1} + float2{2,3}; // v is float2{3,4}`
* supports mixed element types: `auto v = float3{1,2,3} + int3{4,5,6}; // v is float3{5,7,9}`
* supports binary operators with a scalar on the left or the right: `auto v = float3{1,2,3} * 2; // v is float3{2,4,6}`
* supports operators `+=`, `-=`, `*=`, `/=`, `%=`, `|=`, `&=`, `^=`, `<<=` and `>>=` with vectors or scalars on the right
* is iterable: `for(auto elem : float3{1,2,3}) cout << elem; // prints "123"`
* has a flat memory layout: `float3 v {1,2,3}; float * p = v.data(); // &v[i] == p+i`

#### Matrices

TODO: Explain `linalg::mat<T,M,N>`

#### Quaternions

TODO: Explain `linalg::quat<T>`

## Function listing

#### Vector algebra

#### Matrix algebra

#### Quaternion algebra

#### Component-wise operations

#### Reductions

## Optional features

#### Type aliases

TODO: Explain `namespace linalg::aliases`

#### `ostream` overloads

TODO: Explain `namespace linalg::ostream_overloads`

#### User-defined conversions

TODO: Explain `converter<T,U>`

#### Extensions in `linalgx.h`

A second header `linalgx.h` is under development, providing functionality that is not essential to the use of `linalg.h`. This header will be considered unstable even after the publication of `v3.0`, but should stabilize shortly thereafter.

## Higher order functions

TODO: Explain `apply` and `fold`

## Design rationale

This space exists for me to explain the design decisions that went into `linalg.h`.

## Changes from `v2.1`

TODO: Rework this to include a migration guide for users coming from previous versions of `linalg.h`.

#### Improvements in `v3.0`

* `mat<T,M,N>` now defines `operator *` as the matrix product
* New type `quat<T>` models quaternions, and defines `operator *` as the quaternion product
* `map(a,f)` and `zip(a,b,f)` subsumed by new `apply(f,a...)`
  * `apply(...)` supports unary, binary, and ternary operations for `vec`
  * `apply(...)` supports unary and binary operations for `mat` and `quat`
  * `apply(...)` can also be invoked exclusively with scalars, and supports arbitrary numbers of arguments
  * `apply(...)` supports mixed element types
  * Template type alias `apply_t<F,A...>` provides the return type of `apply(f,a...)`
* `vec<T,1>` and `mat<T,M,1>` specializations are now provided
* `compare(a,b)` provide three-way comparison between compatible types
* `clamp(a,b,c)` can be invoked with three distinct (but compatible) types
* `select(a,b,c)` provides the a component-wise equivalent to `a ? b : c`
* `lerp(a,b,t)` has been generalized to a component-wise operation where any of `a`, `b`, and `t` can be vectors or scalars
* `vec<T,M>` elements can be referred to via `x`,`y`,`z`,`w` or `r`,`g`,`b`,`a` or `s`,`t`,`p`,`q`
* `vec<T,M>` and `mat<T,M,N>` contain a member function `.data()` which returns a pointer to the elements
* Groups of `vec<T,M>` elements can be accessed via named swizzles, such as `xyz`, `bgr`, or `pq`
  * All swizzles can be used as rvalues
  * Only permutation swizzles (where no accessor is repeated) can be used as lvalues
* User can specialize `converter<T,U>` to enable implicit conversions from `U` to `T`, if either type is a `vec`, `mat`, or `quat`
  * `identity` is implemented using this facility to serve as an in-library example
* No undefined behavior according to the C++11 standard
* Almost all operations which do not internally call `<cmath>` functions are `constexpr`, except for `argmin` and `argmax`
* No lambdas are used in `linalg.h`, avoidng potential ODR violations

#### Breaking changes in `v3.0`

* `linalg.h` no longer supports Visual Studio 2013. However, it is known to work on GCC 4.9+, Clang 3.5+ in C++11 mode and Visual Studio 2015+.
* `vec<T,M>` and `mat<T,M,N>` may only be used with a `T` which is an [arithmetic type](https://en.cppreference.com/w/c/language/arithmetic_types)
  * This requirement will likely be relaxed by the official release of `v3.0` but will require specializing some trait type to indicate additional scalar types
* As stated above, the semantics of `operator *` for `mat<T,M,N>` has changed
  * Additional, the set of operator overloads on matrices has been reduced to those with an obvious arithmetic meaning.
  * `mat+mat`, `mat-mat`, `mat*scalar`, `scalar*mat`, and `mat/scalar` are implemented as elementwise operations
  * `mat*mat` and `mat*vec` computes the matrix product
  * All other binary matrix overloads have been deleted
  * Most componentwise and reduction functions no longer apply to matrices
  * However, componentwise operations can still be invoked manually with the `apply(...)` function
* `vec<T,M>` and `mat<T,M,N>` are now defined in terms of arrays, instead of explicit fields `x`, `y`, `z`, `w`
  * However, `vec<T,M>` can still be accessed with `x`,`y`,`z`,`w` due to special scalar accessor types which should silently convert to `T` in most scenarios
* Some constructors have been removed from `vec` and `mat`
  * `vec<T,M>` no longer has an implicit constructor from `std::array<T,M>`
  * `vec<T,M>` and `mat<T,M,N>` no longer have an explicit constructor from `const T *`
  * These capabilities can be added by specializing `converter<T,U>`, as shown in [test-user-defined-conversions.cpp](tests/test-user-defined-conversions.cpp) 
* The functions `vec::xy()` and `vec::xyz()` have been replaced by the swizzles `vec::xy` and `vec::xyz`
* Some functionality has been moved from `linalg.h` to optional `linalgx.h` header
  * quat/matrix factory functions for 3D transformations
  * std::hash<...> specializations

## Documentation checklist

Documentation needs to be provided for the following symbols.

- [x] `struct vec<T,M>`: `elems`, accessors, swizzles, constructors, `operator[]`
- [ ] `struct mat<T,M,N>`: `cols`, constructors, `operator[]`, `row`
- [ ] `struct quat<T>`: `x`, `y`, `z`, `w`, constructors, `xyz`
- [ ] user-defined conversions: `converter<T,U>`
- [ ] `identity`
- [ ] higher-order functions: `fold`, `apply`, `map`, `zip`, `apply_t`
- [ ] three-way comparison: `compare`
- [ ] [`EqualityComparable`](http://en.cppreference.com/w/cpp/concept/EqualityComparable): `operator ==, !=`
- [ ] [`LessThanComparable`](http://en.cppreference.com/w/cpp/concept/LessThanComparable): `operator <, >, <=, >=`
- [ ] component-wise operator overloads for `vec<T,M>`
- [ ] algebraic operator overloads for `mat<T,M,N>` and `quat<T>`
- [ ] reduction functions: `any`, `all`, `sum`, `product`, `minelem`, `maxelem`
- [ ] search functions: `argmin`, `argmax`
- [ ] `<cmath>` projections: `abs`, `floor`, `ceil`, `exp`, `log`, `log10`, `sqrt`, `sin`, `cos`, `tan`, `asin`, `acos`, `atan`, `sinh`, `cosh`, `tanh`, `round`, `fmod`, `pow`, `atan2`, `copysign`
- [ ] component-wise comparisons: `equal`, `nequal`, `less`, `greater`, `lequal`, `gequal`
- [ ] component-wise selection: `min`, `max`, `clamp`, `select`
- [ ] vector algebra: `cross`, `dot`, `length`, `length2`, `normalize`, `distance`, `distance2`, `angle`, `uangle`, `rot`, `lerp`, `nlerp`, `slerp`
- [ ] matrix algebra: `diagonal`, `outerprod`, `transpose`, `adjugate`, `determinant`, `trace`, `inverse`
- [ ] quaternion algebra: `conjugate`, `dot`, `length`, `length2`, `inverse`, `normalize`, `uangle`, `lerp`, `nlerp`, `slerp`, `qexp`, `qlog`, `qpow`
- [ ] rotation quaternion support: `qxdir`, `qydir`, `qzdir`, `qmat`, `qrot`, `qangle`, `qaxis`, `qnlerp`, `qslerp`
- [ ] iterators and ranges: `begin`, `end`
- [ ] `namespace linalg::aliases`
