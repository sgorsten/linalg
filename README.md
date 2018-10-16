# linalg.h

NOTE: This is the development branch for `linalg.h v3.0`. It should be considered **unstable** and breaking changes **will occur** until the `v3.0` tag is published. Furthermore, the documentation should be assumed to be out of date until this notice is removed. 

[![Release is 3.0-alpha](http://img.shields.io/badge/version-3.0--alpha-blue.svg)](http://raw.githubusercontent.com/sgorsten/linalg/v3/linalg.h)
[![License is Unlicense](http://img.shields.io/badge/license-Unlicense-blue.svg?style=flat)](http://unlicense.org/)

Platform | Build Status |
-------- | ------------ |
GCC 4.9-8.0 Clang 3.5-6.0 | [Travis CI](http://travis-ci.org): [![Build status](http://travis-ci.org/sgorsten/linalg.svg?branch=v3)](https://travis-ci.org/sgorsten/linalg) |
Visual Studio 2015, 2017 | [AppVeyor](http://ci.appveyor.com/): [![Build status](http://ci.appveyor.com/api/projects/status/l4bfv5omodkajuc9?svg=true)](https://ci.appveyor.com/project/sgorsten/linalg) |

[linalg.h](/linalg.h) is a [single header](http://github.com/nothings/stb/blob/master/docs/other_libs.md) [public domain](http://unlicense.org/) [linear algebra](http://en.wikipedia.org/wiki/Linear_algebra) library for [C++11](http://en.cppreference.com/w/). 

It is inspired by the syntax of popular shader languages and intended to serve as a lightweight (around 750 total lines of code) alternative to projects such as [GLM](http://glm.g-truc.net/0.9.7/), [Boost.QVM](https://www.boost.org/doc/libs/1_66_0/libs/qvm/doc/index.html) or [Eigen](http://eigen.tuxfamily.org/) in domains such as computer graphics, computational geometry, and physical simulation. It aims to be correct, complete, easy to use, readable, and quick to compile.

# Improvements in `v3.0`

* `mat<T,M,N>` now defines `operator *` as the matrix product
* New type `quat<T>` models quaternions, and defines `operator *` as the quaternion product
* `map(a,f)` and `zip(a,b,f)` subsumed by new `apply(f,a...)`
  * `apply(...)` supports unary, binary, and ternary operations for `vec`
  * `apply(...)` supports unary and binary operations for `mat` and `quat`
  * `apply(...)` can also be invoked exclusively with scalars, and supports arbitrary numbers of arguments
  * `apply(...)` supports mixed element types
  * Template type alias `apply_t<F,A...>` provides the return type of `apply(f,a...)`
* `compare(a,b)` provide three-way comparison between compatible types
* `clamp(a,b,c)` can be invoked with three distinct (but compatible) types
* `select(a,b,c)` provides the a component-wise equivalent to `a ? b : c`
* `lerp(a,b,t)` has been generalized to a component-wise operation where any of `a`, `b`, and `t` can be vectors or scalars
* `vec<T,M>` elements can be referred to via `x`,`y`,`z`,`w` or `r`,`g`,`b`,`a` or `s`,`t`,`p`,`q`
* Groups of `vec<T,M>` elements can be accessed via named swizzles, such as `xyz`, `bgr`, or `pq`
  * For now, all permutations of elements are supported, but no duplicates are allowed (such as `xyzz`)
  * This is to allow swizzles to be used as lvalues without introducing ambiguous meanings
* User can specialize `converter<T,U>` to enable implicit conversions from `U` to `T`, if either type is a `vec`, `mat`, or `quat`
  * `identity` is implemented using this facility to serve as an in-library example
* No undefined behavior according to the C++11 standard
* Almost all operations which do not internally call `<cmath>` functions are `constexpr`, except for `argmin` and `argmax`
* No lambdas are used in `linalg.h`, avoidng potential ODR violations


# Breaking changes in `v3.0`

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

# Documentation

This section has been removed pending a complete rewrite. Documentation needs to be provided for the following symbols.

- [ ] `struct vec<T,M>`: `elems`, accessors, swizzles, constructors, `operator[]`
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
- [ ] matrix algebra: `diagonal`, `outerprod`, `transpose`, `adjugate`, `determinant`, `inverse`
- [ ] quaternion algebra: `conjugate`, `dot`, `length`, `length2`, `inverse`, `normalize`, `uangle`, `lerp`, `nlerp`, `slerp`, `qexp`, `qlog`, `qpow`
- [ ] rotation quaternion support: `qxdir`, `qydir`, `qzdir`, `qmat`, `qrot`, `qangle`, `qaxis`, `qnlerp`, `qslerp`
- [ ] transformations support: `rotation_quat`, `translation_matrix`, `rotation_matrix`, `scaling_matrix`, `pose_matrix`, `frustum_matrix`, `perspective_matrix`, `fwd_axis`, `z_range`
- [ ] iterators and ranges: `begin`, `end`
- [ ] `std::hash<...>` specializations
- [ ] `namespace linalg::aliases`
