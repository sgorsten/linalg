# linalg.h

[![Release is 2.0](http://img.shields.io/badge/release-2.0-blue.svg?style=flat)](https://raw.githubusercontent.com/sgorsten/linalg/master/linalg.h)
[![License is Unlicense](http://img.shields.io/badge/license-Unlicense-blue.svg?style=flat)](http://unlicense.org/)

Platform | Build Status |
-------- | ------------ |
Visual Studio 2013, 2015, and 2017 | [AppVeyor](http://ci.appveyor.com/): [![Build status](http://ci.appveyor.com/api/projects/status/l4bfv5omodkajuc9?svg=true)](https://ci.appveyor.com/project/sgorsten/linalg) |
GCC 4.9 and Clang 3.7 | [Travis CI](http://travis-ci.org): [![Build status](http://travis-ci.org/sgorsten/linalg.svg?branch=master)](https://travis-ci.org/sgorsten/linalg) |

[linalg.h](/linalg.h) is a [single header](http://github.com/nothings/stb/blob/master/docs/other_libs.md) [public domain](http://unlicense.org/) [linear algebra](http://en.wikipedia.org/wiki/Linear_algebra) library for [C++11](http://en.cppreference.com/w/). 

It is inspired by the syntax of popular shader languages and intended to serve as a lightweight (less than 400 total lines of code) alternative to projects such as [GLM](http://glm.g-truc.net/0.9.7/) or [Eigen](http://eigen.tuxfamily.org/) in domains such as computer graphics, computational geometry, and physical simulation. It aims to be correct, complete, easy to use, readable, and quick to compile.

# FAQ

###### Why another linear algebra library?

Existing linear algebra libraries are good but most are rather large, slowing down compile times and complicating inclusion into projects. `linalg.h` is a single file designed to be dropped directly into your source tree, and imposes no restrictions on your software from a technical, architectural, or legal standpoint.

###### Why C++11?

Mostly due to broad availability of mostly compliant C++11 compilers. Earlier versions of C++ lack the features (lambdas, decltype, braced initializer lists, etc.) needed to implement `linalg.h` as generically and tersely as it has been. Later versions of C++ do provide useful features which could be used to make `linalg.h` even smaller and cleaner (generic lambdas and auto return types in particular), but do not appreciably improve the functionality of the library in its current form.

###### Why doesn't `operator *` perform matrix multiplication?

Most operator overloads and many function definitions in `linalg.h` use only a single line of code to define vector/vector, vector/scalar, scalar/vector, matrix/matrix, matrix/scalar, and scalar/matrix variations, for all possible element types and all dimensions of vector and matrix, and provide the behavior of applying the given operation to corresponding pairs of elements to produce a vector or matrix valued result. I chose to implement `operator *` in terms of elementwise multiplication for consistency with the rest of the library, and defined the `mul` function to provide matrix multiplication, alongside `dot`, `cross`, and `qmul`.

###### In that case, why do `operator ==`, `operator <` return a `bool` instead of a vector or matrix?

I wanted all instances of `linalg.h` types to model value semantics, and satisfy the [`EqualityComparable`](http://en.cppreference.com/w/cpp/concept/EqualityComparable) and [`LessThanComparable`](http://en.cppreference.com/w/cpp/concept/LessThanComparable) concepts from the C++ standard library. This means that code like `if(a == b) ...` behaves as it typically would, and data structures like `std::set` and `std::map` and functions like `std::find` and `std::sort` can be used with `linalg.h` types without modification.

# Documentation

* [Data Structures](#data-structures)
* [Relational Operators](#relational-operators)
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

A variety of useful `typedef`s are provided in `namespace linalg::aliases`, which can be brought into scope with a `using` declaration. The typedefs for `float` based vectors and matrices are shown below:

|          | vector   | 2 columns  | 3 columns  | 4 columns  |
|----------|----------|------------|------------|------------|
| *2 rows* | `float2` | `float2x2` | `float2x3` | `float2x4` |
| *3 rows* | `float3` | `float3x2` | `float3x3` | `float3x4` |
| *4 rows* | `float4` | `float4x2` | `float4x3` | `float4x4` |

The general pattern for vectors and matrices of other types are shown in the following table:

| underlying type | `vec<T,M>` typedef | `mat<T,M,N>` typedef |
|-----------------|--------------------|----------------------|
| `float`         | `floatM`           | `floatMxN`           |
| `double`        | `doubleM`          | `doubleMxN`          |
| `int`           | `intM`             | `intMxN`             |
| `bool`          | `boolM `           | `boolMxN`            |
| `unsigned`      | `uintM`            |                      |
| `int16_t`       | `shortM`           |                      |
| `uint16_t`      | `ushortM`          |                      |
| `uint8_t`       | `byteM`            |                      |

## Relational Operators

The equivalence and relational operators on `vec<T,M>` are defined as though it were a `std::array<T,M>`. The equivalence and relational operators on `mat<T,M,N>` are defined as though it were a `std::array<T,M*N>`, with the elements laid out in column-major order. Therefore, both types satisfy the [`EqualityComparable`](http://en.cppreference.com/w/cpp/concept/EqualityComparable) and [`LessThanComparable`](http://en.cppreference.com/w/cpp/concept/LessThanComparable) concepts from the C++ standard library, and are suitable for use as the key type in `std::set`, `std::map`, etc.

Additionally, specializations are provided for `std::hash<T>` for both `vec<T,M>` and `mat<T,M,N>`, making them suitable for use as the key type in `std::unordered_set` and `std::unordered_map`.

## Elementwise Functions

A large number of functions and operator overloads exist which interpret vectors and matrices simply as fixed-size blocks of numbers. These operations will simply apply operations to each component of a vector or matrix, or to componentwise pairs of vectors and matrices of compatible size, and produce a new vector or matrix of the same size containing the results.

### Unary Operations

For any vector or matrix `a`, the following unary operations will result in a vector or matrix of the same size:

* `+a` applies the unary `operator +` to each element from `a`
* `-a` applies the unary `operator -` to each element from `a`
* `~a` applies the `operator ~` to each element from `a`
* `!a` applies the `operator !` to each element from `a`, producing a vector or matrix of bools
* `abs(a)` applies `std::abs(...)` to each element from `a`
* `floor(a)` applies `std::floor(...)` to each element from `a`
* `ceil(a)` applies `std::ceil(...)` to each element from `a`
* `exp(a)` applies `std::exp(...)` to each element from `a`
* `log(a)` applies `std::log(...)` to each element from `a`
* `log10(a)` applies `std::log10(...)` to each element from `a`
* `sqrt(a)` applies `std::sqrt(...)` to each element from `a`
* `sin(a)` applies `std::sin(...)` to each element from `a`
* `cos(a)` applies `std::cos(...)` to each element from `a`
* `tan(a)` applies `std::tan(...)` to each element from `a`
* `asin(a)` applies `std::asin(...)` to each element from `a`
* `acos(a)` applies `std::acos(...)` to each element from `a`
* `atan(a)` applies `std::atan(...)` to each element from `a`
* `sinh(a)` applies `std::sinh(...)` to each element from `a`
* `cosh(a)` applies `std::cosh(...)` to each element from `a`
* `tanh(a)` applies `std::tanh(...)` to each element from `a`
* `round(a)` applies `std::round(...)` to each element from `a`

### Binary Operations

For values `a` and `b`, there are a number of binary operations available, which produce vectors or matrices by performing the operation on elementwise pairs from `a` and `b`. If either `a` or `b` (but not both) are a scalar, the scalar is paired with each element from the other value, as described in the following table:

type of `a`  | type of `b`  | `f(a,b)` yields |
------------ | ------------ | ----------------|
`vec<T,M>`   | `vec<T,M>`   | `vec<T,M> { f(a[0], b[0]), f(a[1], b[1]), ... }` |
`vec<T,M>`   | `T`          | `vec<T,M> { f(a[0], b), f(a[1], b), ... }` |
`T`          | `vec<T,M>`   | `vec<T,M> { f(a, b[0]), f(a, b[1]), ... }` |
`mat<T,M,N>` | `mat<T,M,N>` | `mat<T,M,N> { {f(a[0][0], b[0][0]), f(a[0][1], b[0][1]), ...}, ... }` |
`mat<T,M,N>` | `T`          | `mat<T,M,N> { {f(a[0][0], b), f(a[0][1], b), ...}, ... }` |
`T`          | `mat<T,M,N>` | `mat<T,M,N> { {f(a, b[0][0]), f(a, b[0][1]), ...}, ... }` |

The following operations are available:

* `a+b` applies the binary `operator +` to componentwise pairs of elements from `a` and `b`
* `a-b` applies the binary `operator -` to componentwise pairs of elements from `a` and `b`
* `a*b` applies the `operator *` to componentwise pairs of elements from `a` and `b`
* `a/b` applies the `operator /` to componentwise pairs of elements from `a` and `b`
* `a%b` applies the `operator %` to componentwise pairs of elements from `a` and `b`
* `a|b` applies the `operator |` to componentwise pairs of elements from `a` and `b`
* `a^b` applies the `operator ^` to componentwise pairs of elements from `a` and `b`
* `a&b` applies the `operator &` to componentwise pairs of elements from `a` and `b`
* `a<<b` applies the `operator <<` to componentwise pairs of elements from `a` and `b`
* `a>>b` applies the `operator >>` to componentwise pairs of elements from `a` and `b`
* `min(a,b)` is equivalent to applying `std::min(...)` to componentwise pairs of elements from `a` and `b`
* `max(a,b)` is equivalent to applying `std::max(...)` to componentwise pairs of elements from `a` and `b`
* `fmod(a,b)` applies `std::fmod(...)` to componentwise pairs of elements from `a` and `b`
* `pow(a,b)` applies `std::pow(...)` to componentwise pairs of elements from `a` and `b`
* `atan2(a,b)` applies `std::atan2(...)` to componentwise pairs of elements from `a` and `b`
* `copysign(a,b)` applies `std::copysign(...)` to componentwise pairs of elements from `a` and `b`

* `equal(a,b)` applies the `operator ==` to componentwise pairs of elements from `a` and `b`, producing a vector or matrix of bools
* `nequal(a,b)` applies the `operator !=` to componentwise pairs of elements from `a` and `b`, producing a vector or matrix of bools
* `less(a,b)` applies the `operator <` to componentwise pairs of elements from `a` and `b`, producing a vector or matrix of bools
* `greater(a,b)` applies the `operator >` to componentwise pairs of elements from `a` and `b`, producing a vector or matrix of bools
* `lequal(a,b)` applies the `operator <=` to componentwise pairs of elements from `a` and `b`, producing a vector or matrix of bools
* `gequal(a,b)` applies the `operator >=` to componentwise pairs of elements from `a` and `b`, producing a vector or matrix of bools

### Ternary operations

* `clamp(a,b,c)` clamps the elements of `a` to the lower bound `b` and the upper bound `c`

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
* `uangle(a,b)` computes the angle, in radians, between unit length vectors `a` and `b`
* `angle(a,b)` computes the angle, in radians, between nonzero length vectors `a` and `b`
* `lerp(a,b,t)` linearly interpolates between `a` and `b` using parameter `t`
* `nlerp(a,b,t)` is equivalent to `normalize(lerp(a,b,t))`
* `slerp(a,b,t)` performs spherical linear interpolation between unit length vectors `a` and `b` using parameter `t`
* `outerprod(a,b)` compute the outer product of `a` and `b`

## Quaternion Algebra

These functions assume that a `vec<T,4>` represents a quaternion, expressed as `xi + yj + zk + w`. Note that quaternion multiplication is explicitly denoted via the function `qmul`, as `operator *` already refers to elementwise multiplication of two vectors.

* `qmul(a,b)` computes the product `ab` of quaternions `a` and `b`
* `qinv(q)` computes the multiplicative inverse of quaternion `q`
* `qconj(q)` computes `q*`, the conjugate of quaternion `q`

Additionally, there are several functions which assume that a quaternion `q` represents a spatial rotation in 3D space, which transforms a vector `v` via the formula `qvq*`. The unit length quaternions form a double cover over spatial rotations. Therefore, the following functions assume quaternion parameters are of unit length and treat `q` as logically identical to `-q`.

* `qangle(q)` computes the angle of rotation for quaternion `q`, in radians
* `qaxis(q)` computes the axis of rotation for quaternion `q`
* `qnlerp(a,b,t)` performs normalized linear interpolation between the spatial rotations represented by `a` and `b` using parameter `t`
* `qslerp(a,b,t)` performs spherical linear interpolation between the spatial rotations represented by `a` and `b` using parameter `t`
* `qrot(q,v)` computes the result of rotating the vector `v` by quaternion `q`
* `qxdir(q)` computes the result of rotating the vector `{1,0,0}` by quaternion `q`
* `qydir(q)` computes the result of rotating the vector `{0,1,0}` by quaternion `q`
* `qzdir(q)` computes the result of rotating the vector `{0,0,1}` by quaternion `q`
* `qmat(q)` computes a `3`x`3` rotation matrix with the same effect as rotating by quaternion `q`

## Matrix Algebra

These functions assume that a `mat<T,M,N>` represents an `M`x`N` matrix, and a `vec<T,M>` represents an `M`x`1` matrix. Note that matrix multiplication is explicitly denoted via the function `mul`, as `operator *` already refers to elementwise multiplication of two matrices.

* `mul(a,b)` computes the product `ab` of matrices `a` and `b`
* `diagonal(a)` returns the diagonal of square matrix `a` as a vector
* `transpose(a)` computes the transpose of matrix `a`
* `inverse(a)` computes the inverse of matrix `a`
* `determinant(a)` computes the determinant of matrix `a`
* `adjugate(a)` computes the adjugate of matrix `a`, which is the transpose of the cofactor matrix

## Factory Functions

These functions exist for easy interoperability with 3D APIs, which frequently use `4`x`4` homogeneous matrices to represent general 3D transformations, and quaternions to represent 3D rotations.

* `rotation_quat(axis,angle)` constructs a rotation quaternion of `angle` radians about the `axis` vector
* `rotation_quat(matrix)` constructs a rotation quaternion from a `3`x`3` rotation matrix
* `translation_matrix(translation)` constructs a transformation matrix which translates by vector `translation`
* `rotation_matrix(rotation)` constructs a transformation matrix which rotates by quaternion `rotation`
* `scaling_matrix(scaling)` constructs a transformation matrix which scales on the x, y, and z axes by the components of vector `scaling`
* `pose_matrix(q,p)` constructs a transformation matrix which rotates by quaternion `q` and translates by vector `p`
* `frustum_matrix(l,r,b,t,n,f)` constructs a transformation matrix which projects by a specified frustum
* `perspective_matrix(fovy,aspect,n,f)` constructs a transformation matrix for a right handed perspective projection

## Higher Order Functions

The following higher order functions are provided by the library:

* `fold(a, f)` combines the elements of `a` using the binary function `f` in left-to-right order
* `map(a, f)` produces a vector or matrix by passing elements from `a` to unary function `f`
* `zip(a, b, f)` produces a vector or matrix by passing componentwise pairs of elements from `a` and `b` to binary function `f`
