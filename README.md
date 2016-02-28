# linalg.h

![Release is 2.0](https://img.shields.io/badge/release-2.0-blue.svg?style=flat)
[![License is Unlicense](https://img.shields.io/badge/license-Unlicense-blue.svg?style=flat)](http://unlicense.org/)

Platform | Build Status |
-------- | ------------ |
Visual Studio 2013 | [![Build status](https://ci.appveyor.com/api/projects/status/l4bfv5omodkajuc9?svg=true)](https://ci.appveyor.com/project/sgorsten/linalg) |
GCC 4.8 | [![Build status](https://travis-ci.org/sgorsten/linalg.svg?branch=master)](https://travis-ci.org/sgorsten/linalg) |

linalg.h is a single header public domain linear algebra library for C++11. 

It is inspired by the syntax of popular shader languages and intended to serve as a lightweight (less than 500 total lines of code) alternative to projects such as [GLM](http://glm.g-truc.net/0.9.7/) or [Eigen](http://eigen.tuxfamily.org/) in domains such as computer graphics, computational geometry, and physical simulation. It aims to be correct, complete, easy to use, readable, and quick to compile.

# Documentation

* [Data Structures](#data-structures)
* [Operator Overloads](#operator-overloads)
* [Elementwise Functions](#elementwise-functions)
* [Reduction Functions](#reduction-functions)
* [Vector Algebra](#vector-algebra)
* [Quaternion Algebra](#quaternion-algebra)
* [Matrix Algebra](#matrix-algebra)
* [Factory Functions](#factory-functions)
* [Higher Order Functions](#higher-order-functions)

## Data Structures

The library is built on two fundamental template types, `linalg::vec<T,M>` and `linalg::mat<T,M,N>`, and provides a large set of `typedef`s of commonly used types in the `linalg::aliases` namespace, including the familiar `float3`, `float4x4`, `int2`, `bool4` etc. Library support, and the convenience aliases, are currently provided for vectors of length `2` to `4`, and matrices of between `2` to `4` columns and `2` to `4` rows.

### `vec<T,M>`

`vec<T,M>` represents a fixed-length vector containing exactly `M` elements of type `T`. By convention, it is assumed to have column semantics. The following operations are available:

* `vec<T,M>()` default constructs all elements of the vector
* `vec<T,M>(T, ...)` constructs vector from exactly `M` instances of type `T`
* `explicit vec<T,M>(T s)` constructs all elements of the vector to the scalar `s`
* `explicit vec<T,M>(const T * p)` constructs a vector by copying elements from an array which begins at address `p`
* `explicit vec<T,M>(vec<U,M> v)` constructs a vector by casting all elements of `v` from `U` to `T`
* `T & operator[] (int i)` retrieves the element from the `i`th row of the vector

### `mat<T,M,N>`

`mat<T,M,N>` represents a fixed-sized `M`x`N` matrix, consisting of exactly `N` columns, each represented as an `M` length vector. The following operations are available:

* `mat<T,M,N>()` default constructs all elements of the matrix
* `mat<T,M,N>(vec<T,M>, ...)` constructs matrix from exactly `N` column vectors
* `explicit mat<T,M,N>(T s)` constructs all elements of the matrix to the scalar `s`
* `explicit mat<T,M,N>(const T * p)` constructs a matrix by copying elements in column major order from an array which begins at address `p`
* `explicit mat<T,M,N>(mat<U,M,N> m)` constructs a matrix by casting all elements of `m` from `U` to `T`
* `vec<T,M> & operator[] (int j)` retrieves the `j`th column vector of the matrix
* `vec<T,N> row (int i)` retrieves the `i`th row of the matrix, as a vector

### Convenience Aliases

A variety of useful `typedef`s are provided in `namespace linalg::aliases`, which can be brought into scope with a `using` declaration.

* `floatM` for `M` in `2`,`3`,`4` is an alias for `linalg::vec<float,M>`
* `doubleM` for `M` in `2`,`3`,`4` is an alias for `linalg::vec<double,M>`
* `intM` for `M` in `2`,`3`,`4` is an alias for `linalg::vec<int,M>`
* `uintM` for `M` in `2`,`3`,`4` is an alias for `linalg::vec<unsigned,M>` 
* `shortM` for `M` in `2`,`3`,`4` is an alias for `linalg::vec<int16_t,M>`
* `ushortM` for `M` in `2`,`3`,`4` is an alias for `linalg::vec<uint16_t,M>` 
* `byteM` for `M` in `2`,`3`,`4` is an alias for `linalg::vec<uint8_t,M>`
* `boolM` for `M` in `2`,`3`,`4` is an alias for `linalg::vec<bool,M>`
* `floatMxN` for `M` and `N` in `2`,`3`,`4` is an alias for `linalg::mat<float,M,N>`
* `doubleMxN` for `M` and `N` in `2`,`3`,`4` is an alias for `linalg::mat<double,M,N>`
* `intMxN` for `M` and `N` in `2`,`3`,`4` is an alias for `linalg::mat<int,M,N>`
* `boolMxN` for `M` and `N` in `2`,`3`,`4` is an alias for `linalg::mat<bool,M,N>`

## Operator Overloads

For any operator `$` in the set `+`, `-`, `~`, and `!`, the following operation is supported:

* `vec<T,M> operator $ (const vec<T,M> & a)`

For any operator `$` in the set `+`, `-`, `*`, `/`, `%`, `|`, `^`, `&`, `<<`, and `>>`, the following operations are supported:

* `vec<T,M> operator $ (const vec<T,M> & a, const vec<T,M> & b)`
* `vec<T,M> operator $ (const vec<T,M> & a, T b)`
* `vec<T,M> operator $ (T a, const vec<T,M> & b)`
* `vec<T,M> & operator $= (vec<T,M> & a, const vec<T,M> & b)`
* `vec<T,M> & operator $= (vec<T,M> & a, T b)`

In all cases, the behavior is equivalent to applying the operator `$` to componentwise pairs of elements from `a` and `b`, and producing a new vector or matrix from the results. For the overloads which accept a scalar on the left or right hand sides, the behavior is as though the scalar were replaced by a vector or matrix of compatible size whose elements are all equal to the provided scalar.

The `==` and `!=` operators are defined in terms of exact equivalence, two vectors or matrices compare equal if and only if all componentwise pairs of elements compare equal.

The `<`, `>`, `<=`, and `>=` operators compare two vectors or matrices by lexicographically comparing their elements. This allows for vectors or matrices to be passed to `std::sort` or used as the key type in `std::set` and `std::map`.

## Elementwise Functions

The unary functions `abs`, `floor`, `ceil`, `exp`, `log`, `log10`, `sqrt`, `sin`, `cos`, `tan`, `asin`, `acos`, `atan`, `sinh`, `cosh`, `tanh`, and `round` can be called with any vector type. The result is equivalent to producing a vector by calling the equivalent function from the `std::` namespace on every element of a vector.

The binary functions `min`, `max`, `fmod`, `pow`, and `atan2` can be called with any two matching vectors. The result is equivalent to producing a vector or matrix by calling the equivalent function from the `std::` namespace on componentwise pairs of elements from the vectors.

* `clamp(a,b,c)` clamps the elements of `a` to the lower bound `b` and the upper bound `c`

The binary functions `equal`, `nequal`, `less`, `greater`, `lequal`, and `gequal` can be called with any two matching vectors. The result is a vector or matrix of `bool` type formed by comparing the elements of the vectors with operator `==`, `!=`, `<`, `>`, `<=`, and `>=` respectively.

## Reduction Functions

These functions take a vector type and return a scalar value.

* `any(a)` returns true if any element of `a` is true
* `all(a)` returns true if all elements of `a` are true
* `sum(a)` returns the scalar sum of all elements in `a`, as if written `a[0] + a[1] + ... a[M-1]`
* `product(a)` returns the scalar product of all elements in `a`, as if written `a[0] * a[1] * ... a[M-1]`
* `minelem(a)` returns the value of the smallest element in `a`
* `maxelem(a)` returns the value of the largest element in `a`
* `argmin(a)` returns the zero-based index of the smallest element in `a`
* `argmax(a)` returns the zero-based index of the largest element in `a`

## Vector Algebra

These functions assume that a `vec<T,M>` represents a mathematical vector within an `M`-dimensional vector space.

* `cross(a,b)` computes the cross product of vectors `a` and `b`
* `dot(a,b)` computes the dot product (also known as the inner or scalar product) of vectors `a` and `b`
* `length(a)` computes the length (magnitude) of vector `a`
* `length2(a)` computes the square of the length of vector `a`
* `normalize(a)` computes a vector of unit length with the same direction as `a`
* `distance(a,b)` computes the Euclidean distance between two points `a` and `b`
* `distance2(a,b)` computes the square of the Euclidean distance between two points `a` and `b`
* `lerp(a,b,t)` linearly interpolates between `a` and `b` using parameter `t`
* `nlerp(a,b,t)` is equivalent to `normalize(lerp(a,b,t))`

## Quaternion Algebra

These functions assume that a `vec<T,4>` represents a quaternion, expressed as `xi + yj + zk + w`. Note that quaternion multiplication is explicitly denoted via the function `qmul`, as `operator *` already refers to elementwise multiplication of two vectors.

* `qmul(a,b)` computes the product `ab` of quaternions `a` and `b`
* `qinv(q)` computes the multiplicative inverse of quaternion `q`
* `qconj(q)` computes `q*`, the conjugate of quaternion `q`

Additionally, there are several functions which assume that a quaternion `q` represents a spatial rotation in 3D space, which transforms a vector `v` via the formula `qvq*`.

* `qangle(q)` computes the angle of rotation for quaternion `q`, in radians
* `qaxis(q)` computes the axis of rotation for quaternion `q`
* `qlerp(a,b,t)` interpolates between the spatial rotations represented by `a` and `b` using parameter `t`
* `qrot(q,v)` computes the result of rotating the vector `v` by quaternion `q`
* `qxdir(q)` computes the result of rotating the vector `{1,0,0}` by quaternion `q`
* `qydir(q)` computes the result of rotating the vector `{0,1,0}` by quaternion `q`
* `qzdir(q)` computes the result of rotating the vector `{0,0,1}` by quaternion `q`
* `qmat(q)` computes a `3`x`3` rotation matrix with the same effect as rotating by quaternion `q`

## Matrix Algebra

These functions assume that a `mat<T,M,N>` represents an `M`x`N` matrix, and a `vec<T,M>` represents an `M`x`1` matrix.

* `mul(a,b)` computes the product `ab` of matrices `a` and `b`
* `transpose(a)` computes the transpose of matrix `a`
* `inverse(a)` computes the inverse of matrix `a`
* `determinant(a)` computes the determinant of matrix `a`
* `adjugate(a)` computes the adjugate of matrix `a`, which is the transpose of the cofactor matrix

## Factory Functions

These functions exist for easy interoperability with 3D APIs, which frequently use `4`x`4` homogeneous matrices to represent general 3D transformations, and quaternions to represent 3D rotations.

* `rotation_quat(axis,angle)` constructs a quaternion of `angle` radians about the `axis` vector
* `translation_matrix(translation)` constructs a transformation matrix which translates by vector `translation`
* `rotation_matrix(rotation)` constructs a transformation matrix which rotates by quaternion `rotation`
* `pose_matrix(q,p)` constructs a transformation matrix which rotates by quaternion `q` and translates by vector `p`
* `frustum_matrix(l,r,b,t,n,f)` constructs a transformation matrix which projects by a specified frustum
* `perspective_matrix(fovy,aspect,n,f)` constructs a transformation matrix for a right handed perspective projection

## Higher Order Functions

The following higher order functions are provided by the library:

* `fold(a, f)` combines the elements of `a` using the binary function `f` in left-to-right order
* `zip(a, b, f)` produces a vector or matrix by passing componentwise pairs of elements from `a` and `b` to binary function `f`
