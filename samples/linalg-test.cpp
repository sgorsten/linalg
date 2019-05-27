#include "../linalg.h"
using namespace linalg::aliases;

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "thirdparty/doctest.h"

#include <random>
#include <algorithm>
#include <typeinfo>

#define FLOATING_POINT_TYPES double, float
#define INTEGRAL_TYPES int, short, unsigned int, unsigned short
#define ARITHMETIC_TYPES double, float, int, short, unsigned int, unsigned short

// Facility for retrieving random numbers
class random_number_generator
{
    std::mt19937 rng;
    std::normal_distribution<double> dist_double;
    std::normal_distribution<float> dist_float;
    std::uniform_int_distribution<int> dist_int;
    std::uniform_int_distribution<short> dist_short;
    std::uniform_int_distribution<unsigned> dist_uint;
    std::uniform_int_distribution<unsigned short> dist_ushort;
public:
    random_number_generator() : dist_int(-1000, 1000), dist_short(-100, 100), dist_uint(0, 1000), dist_ushort(0, 100) {}

    operator double () { return dist_double(rng); }
    operator float () { return dist_float(rng); }
    operator int () { return dist_int(rng); }
    operator short () { return dist_short(rng); }
    operator unsigned int () { return dist_uint(rng); }
    operator unsigned short () { return dist_ushort(rng); }
    template<class T> operator linalg::vec<T,1> () { return linalg::vec<T,1>((T)*this); }
    template<class T> operator linalg::vec<T,2> () { return linalg::vec<T,2>((T)*this, (T)*this); }
    template<class T> operator linalg::vec<T,3> () { return linalg::vec<T,3>((T)*this, (T)*this, (T)*this); }
    template<class T> operator linalg::vec<T,4> () { return linalg::vec<T,4>((T)*this, (T)*this, (T)*this, *this); }
    template<class T, int M> operator linalg::mat<T,M,1> () { return linalg::mat<T,M,1>((linalg::vec<T,M>)*this); }
    template<class T, int M> operator linalg::mat<T,M,2> () { return linalg::mat<T,M,2>((linalg::vec<T,M>)*this, (linalg::vec<T,M>)*this); }
    template<class T, int M> operator linalg::mat<T,M,3> () { return linalg::mat<T,M,3>((linalg::vec<T,M>)*this, (linalg::vec<T,M>)*this, (linalg::vec<T,M>)*this); }
    template<class T, int M> operator linalg::mat<T,M,4> () { return linalg::mat<T,M,4>((linalg::vec<T,M>)*this, (linalg::vec<T,M>)*this, (linalg::vec<T,M>)*this, (linalg::vec<T,M>)*this); }
};
static const int reps = 3; // Tests which use random data will be repeated this many times

//////////////////////////////////////////////////////////
// Test semantics of vec<T,M> element-wise constructors //
//////////////////////////////////////////////////////////

TEST_CASE_TEMPLATE("vec<T,1> can be constructed from 1 element of type T", T, ARITHMETIC_TYPES) 
{
    random_number_generator rng;
    for(int i=0; i<reps; ++i)
    {
        const T e0 = rng;
        const linalg::vec<T,1> v(e0);
        CHECK(v.x == e0); CHECK(v[0] == e0);
    }
}

TEST_CASE_TEMPLATE("vec<T,2> can be constructed from 2 elements of type T", T, ARITHMETIC_TYPES) 
{
    random_number_generator rng;
    for(int i=0; i<reps; ++i)
    {
        const T e0 = rng, e1 = rng;
        const linalg::vec<T,2> v(e0,e1);
        CHECK(v.x == e0); CHECK(v[0] == e0);
        CHECK(v.y == e1); CHECK(v[1] == e1);
    }
}

TEST_CASE_TEMPLATE("vec<T,3> can be constructed from 3 elements of type T", T, ARITHMETIC_TYPES) 
{
    random_number_generator rng;
    for(int i=0; i<reps; ++i)
    {
        const T e0 = rng, e1 = rng, e2 = rng;
        const linalg::vec<T,3> v(e0,e1,e2);
        CHECK(v.x == e0); CHECK(v[0] == e0);
        CHECK(v.y == e1); CHECK(v[1] == e1);
        CHECK(v.z == e2); CHECK(v[2] == e2);
    }
}

TEST_CASE_TEMPLATE("vec<T,4> can be constructed from 4 elements of type T", T, ARITHMETIC_TYPES)
{
    random_number_generator rng;
    for(int i=0; i<reps; ++i)
    {
        const T e0 = rng, e1 = rng, e2 = rng, e3 = rng;
        const linalg::vec<T,4> v(e0,e1,e2,e3);
        CHECK(v.x == e0); CHECK(v[0] == e0);
        CHECK(v.y == e1); CHECK(v[1] == e1);
        CHECK(v.z == e2); CHECK(v[2] == e2);
        CHECK(v.w == e3); CHECK(v[3] == e3);
    }
}

///////////////////////////////////////////////////////////
// Test semantics of mat<T,M,N> column-wise constructors //
///////////////////////////////////////////////////////////

TEST_CASE_TEMPLATE("mat<T,1,1> can be constructed from 1 column of type vec<T,1>", T, ARITHMETIC_TYPES)
{
    random_number_generator rng;
    for(int i=0; i<reps; ++i)
    {
        const T m00=rng;
        const auto m = linalg::mat<T,1,1>(
            linalg::vec<T,1>(m00)
        );
        CHECK(m.x.x == m00); CHECK(m[0][0] == m00); CHECK(m.row(0)[0] == m00);
    }
}

TEST_CASE_TEMPLATE("mat<T,2,2> can be constructed from 2 columns of type vec<T,2>", T, ARITHMETIC_TYPES)
{
    random_number_generator rng;
    for(int i=0; i<reps; ++i)
    {
        const T m00=rng, m01=rng;
        const T m10=rng, m11=rng;
        const linalg::mat<T,2,2> m(
            linalg::vec<T,2>(m00,m10),    
            linalg::vec<T,2>(m01,m11)
        );
        CHECK(m.x.x == m00); CHECK(m[0][0] == m00); CHECK(m.row(0)[0] == m00);
        CHECK(m.x.y == m10); CHECK(m[0][1] == m10); CHECK(m.row(1)[0] == m10);
        CHECK(m.y.x == m01); CHECK(m[1][0] == m01); CHECK(m.row(0)[1] == m01);
        CHECK(m.y.y == m11); CHECK(m[1][1] == m11); CHECK(m.row(1)[1] == m11);
    }
}

TEST_CASE_TEMPLATE("mat<T,2,3> can be constructed from 3 columns of type vec<T,2>", T, ARITHMETIC_TYPES)
{
    random_number_generator rng;
    for(int i=0; i<reps; ++i)
    {
        const T m00=rng, m01=rng, m02=rng;
        const T m10=rng, m11=rng, m12=rng;
        const linalg::mat<T,2,3> m(
            linalg::vec<T,2>(m00,m10),    
            linalg::vec<T,2>(m01,m11),
            linalg::vec<T,2>(m02,m12)
        );
        CHECK(m.x.x == m00); CHECK(m[0][0] == m00); CHECK(m.row(0)[0] == m00);
        CHECK(m.x.y == m10); CHECK(m[0][1] == m10); CHECK(m.row(1)[0] == m10);
        CHECK(m.y.x == m01); CHECK(m[1][0] == m01); CHECK(m.row(0)[1] == m01);
        CHECK(m.y.y == m11); CHECK(m[1][1] == m11); CHECK(m.row(1)[1] == m11);
        CHECK(m.z.x == m02); CHECK(m[2][0] == m02); CHECK(m.row(0)[2] == m02);
        CHECK(m.z.y == m12); CHECK(m[2][1] == m12); CHECK(m.row(1)[2] == m12);
    }
}

TEST_CASE_TEMPLATE("mat<T,2,4> can be constructed from 4 columns of type vec<T,2>", T, ARITHMETIC_TYPES)
{
    random_number_generator rng;
    for(int i=0; i<reps; ++i)
    {
        const T m00=rng, m01=rng, m02=rng, m03=rng;
        const T m10=rng, m11=rng, m12=rng, m13=rng;
        const linalg::mat<T,2,4> m(
            linalg::vec<T,2>(m00,m10),    
            linalg::vec<T,2>(m01,m11),
            linalg::vec<T,2>(m02,m12),
            linalg::vec<T,2>(m03,m13)
        );
        CHECK(m.x.x == m00); CHECK(m[0][0] == m00); CHECK(m.row(0)[0] == m00);
        CHECK(m.x.y == m10); CHECK(m[0][1] == m10); CHECK(m.row(1)[0] == m10);
        CHECK(m.y.x == m01); CHECK(m[1][0] == m01); CHECK(m.row(0)[1] == m01);
        CHECK(m.y.y == m11); CHECK(m[1][1] == m11); CHECK(m.row(1)[1] == m11);
        CHECK(m.z.x == m02); CHECK(m[2][0] == m02); CHECK(m.row(0)[2] == m02);
        CHECK(m.z.y == m12); CHECK(m[2][1] == m12); CHECK(m.row(1)[2] == m12);
        CHECK(m.w.x == m03); CHECK(m[3][0] == m03); CHECK(m.row(0)[3] == m03);
        CHECK(m.w.y == m13); CHECK(m[3][1] == m13); CHECK(m.row(1)[3] == m13);
    }
}

TEST_CASE_TEMPLATE("mat<T,3,2> can be constructed from 2 columns of type vec<T,3>", T, ARITHMETIC_TYPES)
{
    random_number_generator rng;
    for(int i=0; i<reps; ++i)
    {
        const T m00=rng, m01=rng;
        const T m10=rng, m11=rng;
        const T m20=rng, m21=rng;
        const linalg::mat<T,3,2> m(
            linalg::vec<T,3>(m00,m10,m20),    
            linalg::vec<T,3>(m01,m11,m21)
        );
        CHECK(m.x.x == m00); CHECK(m[0][0] == m00); CHECK(m.row(0)[0] == m00);
        CHECK(m.x.y == m10); CHECK(m[0][1] == m10); CHECK(m.row(1)[0] == m10);
        CHECK(m.x.z == m20); CHECK(m[0][2] == m20); CHECK(m.row(2)[0] == m20);
        CHECK(m.y.x == m01); CHECK(m[1][0] == m01); CHECK(m.row(0)[1] == m01);
        CHECK(m.y.y == m11); CHECK(m[1][1] == m11); CHECK(m.row(1)[1] == m11);
        CHECK(m.y.z == m21); CHECK(m[1][2] == m21); CHECK(m.row(2)[1] == m21);
    }
}

TEST_CASE_TEMPLATE("mat<T,3,3> can be constructed from 3 columns of type vec<T,3>", T, ARITHMETIC_TYPES)
{
    random_number_generator rng;
    for(int i=0; i<reps; ++i)
    {
        const T m00=rng, m01=rng, m02=rng;
        const T m10=rng, m11=rng, m12=rng;
        const T m20=rng, m21=rng, m22=rng;
        const linalg::mat<T,3,3> m(
            linalg::vec<T,3>(m00,m10,m20),    
            linalg::vec<T,3>(m01,m11,m21),
            linalg::vec<T,3>(m02,m12,m22)
        );
        CHECK(m.x.x == m00); CHECK(m[0][0] == m00); CHECK(m.row(0)[0] == m00);
        CHECK(m.x.y == m10); CHECK(m[0][1] == m10); CHECK(m.row(1)[0] == m10);
        CHECK(m.x.z == m20); CHECK(m[0][2] == m20); CHECK(m.row(2)[0] == m20);
        CHECK(m.y.x == m01); CHECK(m[1][0] == m01); CHECK(m.row(0)[1] == m01);
        CHECK(m.y.y == m11); CHECK(m[1][1] == m11); CHECK(m.row(1)[1] == m11);
        CHECK(m.y.z == m21); CHECK(m[1][2] == m21); CHECK(m.row(2)[1] == m21);
        CHECK(m.z.x == m02); CHECK(m[2][0] == m02); CHECK(m.row(0)[2] == m02);
        CHECK(m.z.y == m12); CHECK(m[2][1] == m12); CHECK(m.row(1)[2] == m12);
        CHECK(m.z.z == m22); CHECK(m[2][2] == m22); CHECK(m.row(2)[2] == m22);
    }
}

TEST_CASE_TEMPLATE("mat<T,3,4> can be constructed from 4 columns of type vec<T,3>", T, ARITHMETIC_TYPES)
{
    random_number_generator rng;
    for(int i=0; i<reps; ++i)
    {
        const T m00=rng, m01=rng, m02=rng, m03=rng;
        const T m10=rng, m11=rng, m12=rng, m13=rng;
        const T m20=rng, m21=rng, m22=rng, m23=rng;
        const linalg::mat<T,3,4> m(
            linalg::vec<T,3>(m00,m10,m20),    
            linalg::vec<T,3>(m01,m11,m21),
            linalg::vec<T,3>(m02,m12,m22),
            linalg::vec<T,3>(m03,m13,m23)
        );
        CHECK(m.x.x == m00); CHECK(m[0][0] == m00); CHECK(m.row(0)[0] == m00);
        CHECK(m.x.y == m10); CHECK(m[0][1] == m10); CHECK(m.row(1)[0] == m10);
        CHECK(m.x.z == m20); CHECK(m[0][2] == m20); CHECK(m.row(2)[0] == m20);
        CHECK(m.y.x == m01); CHECK(m[1][0] == m01); CHECK(m.row(0)[1] == m01);
        CHECK(m.y.y == m11); CHECK(m[1][1] == m11); CHECK(m.row(1)[1] == m11);
        CHECK(m.y.z == m21); CHECK(m[1][2] == m21); CHECK(m.row(2)[1] == m21);
        CHECK(m.z.x == m02); CHECK(m[2][0] == m02); CHECK(m.row(0)[2] == m02);
        CHECK(m.z.y == m12); CHECK(m[2][1] == m12); CHECK(m.row(1)[2] == m12);
        CHECK(m.z.z == m22); CHECK(m[2][2] == m22); CHECK(m.row(2)[2] == m22);
        CHECK(m.w.x == m03); CHECK(m[3][0] == m03); CHECK(m.row(0)[3] == m03);
        CHECK(m.w.y == m13); CHECK(m[3][1] == m13); CHECK(m.row(1)[3] == m13);
        CHECK(m.w.z == m23); CHECK(m[3][2] == m23); CHECK(m.row(2)[3] == m23);
    }
}

TEST_CASE_TEMPLATE("mat<T,4,2> can be constructed from 2 columns of type vec<T,4>", T, ARITHMETIC_TYPES)
{
    random_number_generator rng;
    for(int i=0; i<reps; ++i)
    {
        const T m00=rng, m01=rng;
        const T m10=rng, m11=rng;
        const T m20=rng, m21=rng;
        const T m30=rng, m31=rng;
        const linalg::mat<T,4,2> m(
            linalg::vec<T,4>(m00,m10,m20,m30),    
            linalg::vec<T,4>(m01,m11,m21,m31)
        );
        CHECK(m.x.x == m00); CHECK(m[0][0] == m00); CHECK(m.row(0)[0] == m00);
        CHECK(m.x.y == m10); CHECK(m[0][1] == m10); CHECK(m.row(1)[0] == m10);
        CHECK(m.x.z == m20); CHECK(m[0][2] == m20); CHECK(m.row(2)[0] == m20);
        CHECK(m.x.w == m30); CHECK(m[0][3] == m30); CHECK(m.row(3)[0] == m30);
        CHECK(m.y.x == m01); CHECK(m[1][0] == m01); CHECK(m.row(0)[1] == m01);
        CHECK(m.y.y == m11); CHECK(m[1][1] == m11); CHECK(m.row(1)[1] == m11);
        CHECK(m.y.z == m21); CHECK(m[1][2] == m21); CHECK(m.row(2)[1] == m21);
        CHECK(m.y.w == m31); CHECK(m[1][3] == m31); CHECK(m.row(3)[1] == m31);
    }
}

TEST_CASE_TEMPLATE("mat<T,4,3> can be constructed from 3 columns of type vec<T,4>", T, ARITHMETIC_TYPES)
{
    random_number_generator rng;
    for(int i=0; i<reps; ++i)
    {
        const T m00=rng, m01=rng, m02=rng;
        const T m10=rng, m11=rng, m12=rng;
        const T m20=rng, m21=rng, m22=rng;
        const T m30=rng, m31=rng, m32=rng;
        const linalg::mat<T,4,3> m(
            linalg::vec<T,4>(m00,m10,m20,m30),    
            linalg::vec<T,4>(m01,m11,m21,m31),
            linalg::vec<T,4>(m02,m12,m22,m32)
        );
        CHECK(m.x.x == m00); CHECK(m[0][0] == m00); CHECK(m.row(0)[0] == m00);
        CHECK(m.x.y == m10); CHECK(m[0][1] == m10); CHECK(m.row(1)[0] == m10);
        CHECK(m.x.z == m20); CHECK(m[0][2] == m20); CHECK(m.row(2)[0] == m20);
        CHECK(m.x.w == m30); CHECK(m[0][3] == m30); CHECK(m.row(3)[0] == m30);
        CHECK(m.y.x == m01); CHECK(m[1][0] == m01); CHECK(m.row(0)[1] == m01);
        CHECK(m.y.y == m11); CHECK(m[1][1] == m11); CHECK(m.row(1)[1] == m11);
        CHECK(m.y.z == m21); CHECK(m[1][2] == m21); CHECK(m.row(2)[1] == m21);
        CHECK(m.y.w == m31); CHECK(m[1][3] == m31); CHECK(m.row(3)[1] == m31);
        CHECK(m.z.x == m02); CHECK(m[2][0] == m02); CHECK(m.row(0)[2] == m02);
        CHECK(m.z.y == m12); CHECK(m[2][1] == m12); CHECK(m.row(1)[2] == m12);
        CHECK(m.z.z == m22); CHECK(m[2][2] == m22); CHECK(m.row(2)[2] == m22);
        CHECK(m.z.w == m32); CHECK(m[2][3] == m32); CHECK(m.row(3)[2] == m32);
    }
}

TEST_CASE_TEMPLATE("mat<T,4,4> can be constructed from 4 columns of type vec<T,4>", T, ARITHMETIC_TYPES)
{
    random_number_generator rng;
    for(int i=0; i<reps; ++i)
    {
        const T m00=rng, m01=rng, m02=rng, m03=rng;
        const T m10=rng, m11=rng, m12=rng, m13=rng;
        const T m20=rng, m21=rng, m22=rng, m23=rng;
        const T m30=rng, m31=rng, m32=rng, m33=rng;
        const linalg::mat<T,4,4> m(
            linalg::vec<T,4>(m00,m10,m20,m30),    
            linalg::vec<T,4>(m01,m11,m21,m31),
            linalg::vec<T,4>(m02,m12,m22,m32),
            linalg::vec<T,4>(m03,m13,m23,m33)
        );
        CHECK(m.x.x == m00); CHECK(m[0][0] == m00); CHECK(m.row(0)[0] == m00);
        CHECK(m.x.y == m10); CHECK(m[0][1] == m10); CHECK(m.row(1)[0] == m10);
        CHECK(m.x.z == m20); CHECK(m[0][2] == m20); CHECK(m.row(2)[0] == m20);
        CHECK(m.x.w == m30); CHECK(m[0][3] == m30); CHECK(m.row(3)[0] == m30);
        CHECK(m.y.x == m01); CHECK(m[1][0] == m01); CHECK(m.row(0)[1] == m01);
        CHECK(m.y.y == m11); CHECK(m[1][1] == m11); CHECK(m.row(1)[1] == m11);
        CHECK(m.y.z == m21); CHECK(m[1][2] == m21); CHECK(m.row(2)[1] == m21);
        CHECK(m.y.w == m31); CHECK(m[1][3] == m31); CHECK(m.row(3)[1] == m31);
        CHECK(m.z.x == m02); CHECK(m[2][0] == m02); CHECK(m.row(0)[2] == m02);
        CHECK(m.z.y == m12); CHECK(m[2][1] == m12); CHECK(m.row(1)[2] == m12);
        CHECK(m.z.z == m22); CHECK(m[2][2] == m22); CHECK(m.row(2)[2] == m22);
        CHECK(m.z.w == m32); CHECK(m[2][3] == m32); CHECK(m.row(3)[2] == m32);
        CHECK(m.w.x == m03); CHECK(m[3][0] == m03); CHECK(m.row(0)[3] == m03);
        CHECK(m.w.y == m13); CHECK(m[3][1] == m13); CHECK(m.row(1)[3] == m13);
        CHECK(m.w.z == m23); CHECK(m[3][2] == m23); CHECK(m.row(2)[3] == m23);
        CHECK(m.w.w == m33); CHECK(m[3][3] == m33); CHECK(m.row(3)[3] == m33);
    }
}

///////////////////////////////////////////////////
// Test semantics of operator == and operator != //
///////////////////////////////////////////////////

TEST_CASE_TEMPLATE("vec<T,M> a and b compare equal if they contain exactly the same elements", T, ARITHMETIC_TYPES)
{
    random_number_generator rng;
    for(int i=0; i<reps; ++i)
    {
        const T a=rng, b=rng, c=rng, d=rng;

        CHECK(linalg::vec<T,1>(a      ) == linalg::vec<T,1>(a      ));
        CHECK(linalg::vec<T,2>(a,b    ) == linalg::vec<T,2>(a,b    ));
        CHECK(linalg::vec<T,3>(a,b,c  ) == linalg::vec<T,3>(a,b,c  ));
        CHECK(linalg::vec<T,4>(a,b,c,d) == linalg::vec<T,4>(a,b,c,d));

        CHECK_FALSE(linalg::vec<T,1>(a      ) == linalg::vec<T,1>(a+1      ));
        CHECK_FALSE(linalg::vec<T,2>(a,b    ) == linalg::vec<T,2>(a+1,b    ));
        CHECK_FALSE(linalg::vec<T,2>(a,b    ) == linalg::vec<T,2>(a,b+1    ));
        CHECK_FALSE(linalg::vec<T,3>(a,b,c  ) == linalg::vec<T,3>(a+1,b,c  ));
        CHECK_FALSE(linalg::vec<T,3>(a,b,c  ) == linalg::vec<T,3>(a,b+1,c  ));
        CHECK_FALSE(linalg::vec<T,3>(a,b,c  ) == linalg::vec<T,3>(a,b,c+1  ));
        CHECK_FALSE(linalg::vec<T,4>(a,b,c,d) == linalg::vec<T,4>(a+1,b,c,d));
        CHECK_FALSE(linalg::vec<T,4>(a,b,c,d) == linalg::vec<T,4>(a,b+1,c,d));
        CHECK_FALSE(linalg::vec<T,4>(a,b,c,d) == linalg::vec<T,4>(a,b,c+1,d));
        CHECK_FALSE(linalg::vec<T,4>(a,b,c,d) == linalg::vec<T,4>(a,b,c,d+1));
    }
}

TEST_CASE_TEMPLATE("vec<T,M> a and b compare unequal if at least one element differs between them", T, ARITHMETIC_TYPES)
{
    random_number_generator rng;
    for(int i=0; i<reps; ++i)
    {
        const T a=rng, b=rng, c=rng, d=rng;

        CHECK_FALSE(linalg::vec<T,1>(a      ) != linalg::vec<T,1>(a      ));
        CHECK_FALSE(linalg::vec<T,2>(a,b    ) != linalg::vec<T,2>(a,b    ));
        CHECK_FALSE(linalg::vec<T,3>(a,b,c  ) != linalg::vec<T,3>(a,b,c  ));
        CHECK_FALSE(linalg::vec<T,4>(a,b,c,d) != linalg::vec<T,4>(a,b,c,d));

        CHECK(linalg::vec<T,1>(a      ) != linalg::vec<T,1>(a+1      ));
        CHECK(linalg::vec<T,2>(a,b    ) != linalg::vec<T,2>(a+1,b    ));
        CHECK(linalg::vec<T,2>(a,b    ) != linalg::vec<T,2>(a,b+1    ));
        CHECK(linalg::vec<T,3>(a,b,c  ) != linalg::vec<T,3>(a+1,b,c  ));
        CHECK(linalg::vec<T,3>(a,b,c  ) != linalg::vec<T,3>(a,b+1,c  ));
        CHECK(linalg::vec<T,3>(a,b,c  ) != linalg::vec<T,3>(a,b,c+1  ));
        CHECK(linalg::vec<T,4>(a,b,c,d) != linalg::vec<T,4>(a+1,b,c,d));
        CHECK(linalg::vec<T,4>(a,b,c,d) != linalg::vec<T,4>(a,b+1,c,d));
        CHECK(linalg::vec<T,4>(a,b,c,d) != linalg::vec<T,4>(a,b,c+1,d));
        CHECK(linalg::vec<T,4>(a,b,c,d) != linalg::vec<T,4>(a,b,c,d+1));
    }
}

////////////////////////////////////////////
// Test semantics of default constructors //
////////////////////////////////////////////

TEST_CASE_TEMPLATE("vec<T,M>'s default constructor zero-initializes its elements", T, ARITHMETIC_TYPES) 
{
    const linalg::vec<T,1> v1; CHECK(v1 == linalg::vec<T,1>(0));
    const linalg::vec<T,2> v2; CHECK(v2 == linalg::vec<T,2>(0,0));
    const linalg::vec<T,3> v3; CHECK(v3 == linalg::vec<T,3>(0,0,0));
    const linalg::vec<T,4> v4; CHECK(v4 == linalg::vec<T,4>(0,0,0,0));
}


TEST_CASE_TEMPLATE("mat<T,M,N>'s default constructor zero-initializes its columns", T, ARITHMETIC_TYPES) 
{
    const linalg::vec<T,2> z2(0,0); 
    const linalg::vec<T,3> z3(0,0,0); 
    const linalg::vec<T,4> z4(0,0,0,0); 

    const linalg::mat<T,2,2> m2x2; CHECK(m2x2 == linalg::mat<T,2,2>(z2,z2));
    const linalg::mat<T,2,3> m2x3; CHECK(m2x3 == linalg::mat<T,2,3>(z2,z2,z2));
    const linalg::mat<T,2,4> m2x4; CHECK(m2x4 == linalg::mat<T,2,4>(z2,z2,z2,z2));

    const linalg::mat<T,3,2> m3x2; CHECK(m3x2 == linalg::mat<T,3,2>(z3,z3));
    const linalg::mat<T,3,3> m3x3; CHECK(m3x3 == linalg::mat<T,3,3>(z3,z3,z3));
    const linalg::mat<T,3,4> m3x4; CHECK(m3x4 == linalg::mat<T,3,4>(z3,z3,z3,z3));

    const linalg::mat<T,4,2> m4x2; CHECK(m4x2 == linalg::mat<T,4,2>(z4,z4));
    const linalg::mat<T,4,3> m4x3; CHECK(m4x3 == linalg::mat<T,4,3>(z4,z4,z4));
    const linalg::mat<T,4,4> m4x4; CHECK(m4x4 == linalg::mat<T,4,4>(z4,z4,z4,z4)); 
}

///////////////////////////////////////////
// Test semantics of scalar constructors //
///////////////////////////////////////////

TEST_CASE_TEMPLATE("vec<T,M>'s scalar constructor initializes its elements to the specified scalar", T, ARITHMETIC_TYPES) 
{
    random_number_generator rng;
    for(int i=0; i<reps; ++i)
    {
        const T s = rng;
        const linalg::vec<T,1> v1(s); CHECK(v1 == linalg::vec<T,1>(s));
        const linalg::vec<T,2> v2(s); CHECK(v2 == linalg::vec<T,2>(s,s));
        const linalg::vec<T,3> v3(s); CHECK(v3 == linalg::vec<T,3>(s,s,s));
        const linalg::vec<T,4> v4(s); CHECK(v4 == linalg::vec<T,4>(s,s,s,s));
    }
}

TEST_CASE_TEMPLATE("mat<T,M,N>'s scalar constructor initializes its columns to the specified scalar", T, ARITHMETIC_TYPES) 
{
    random_number_generator rng;
    for(int i=0; i<reps; ++i)
    {
        const T s = rng;
        const linalg::vec<T,2> s2(s,s); 
        const linalg::vec<T,3> s3(s,s,s); 
        const linalg::vec<T,4> s4(s,s,s,s); 

        const linalg::mat<T,2,2> m2x2(s); CHECK(m2x2 == linalg::mat<T,2,2>(s2,s2));
        const linalg::mat<T,2,3> m2x3(s); CHECK(m2x3 == linalg::mat<T,2,3>(s2,s2,s2));
        const linalg::mat<T,2,4> m2x4(s); CHECK(m2x4 == linalg::mat<T,2,4>(s2,s2,s2,s2));

        const linalg::mat<T,3,2> m3x2(s); CHECK(m3x2 == linalg::mat<T,3,2>(s3,s3));
        const linalg::mat<T,3,3> m3x3(s); CHECK(m3x3 == linalg::mat<T,3,3>(s3,s3,s3));
        const linalg::mat<T,3,4> m3x4(s); CHECK(m3x4 == linalg::mat<T,3,4>(s3,s3,s3,s3));

        const linalg::mat<T,4,2> m4x2(s); CHECK(m4x2 == linalg::mat<T,4,2>(s4,s4));
        const linalg::mat<T,4,3> m4x3(s); CHECK(m4x3 == linalg::mat<T,4,3>(s4,s4,s4));
        const linalg::mat<T,4,4> m4x4(s); CHECK(m4x4 == linalg::mat<T,4,4>(s4,s4,s4,s4));
    }
}

////////////////////////////////////////////
// Test semantics of pointer constructors //
////////////////////////////////////////////

TEST_CASE_TEMPLATE("vec<T,M>'s pointer constructor initializes its elements in order from a pointer to contiguous elements in memory", T, ARITHMETIC_TYPES) 
{
    random_number_generator rng;
    for(int i=0; i<reps; ++i)
    {
        const T a[4] {rng, rng, rng, rng};
        const T * const p = a;
        CHECK(linalg::vec<T,2>(p) == linalg::vec<T,2>(a[0], a[1]));
        CHECK(linalg::vec<T,3>(p) == linalg::vec<T,3>(a[0], a[1], a[2]));
        CHECK(linalg::vec<T,4>(p) == linalg::vec<T,4>(a[0], a[1], a[2], a[3]));
    }
}

TEST_CASE_TEMPLATE("mat<T,M,N>'s pointer constructor initializes its elements in column-major order from a pointer to contiguous elements in memory", T, ARITHMETIC_TYPES) 
{
    random_number_generator rng;
    for(int i=0; i<reps; ++i)
    {
        const T a[16] {rng, rng, rng, rng, rng, rng, rng, rng, rng, rng, rng, rng, rng, rng, rng, rng};
        const T * const p = a;

        CHECK(linalg::mat<T,2,2>(p) == linalg::mat<T,2,2>({a[0],a[1]}, {a[2],a[3]}));
        CHECK(linalg::mat<T,2,3>(p) == linalg::mat<T,2,3>({a[0],a[1]}, {a[2],a[3]}, {a[4],a[5]}));
        CHECK(linalg::mat<T,2,4>(p) == linalg::mat<T,2,4>({a[0],a[1]}, {a[2],a[3]}, {a[4],a[5]}, {a[6],a[7]}));

        CHECK(linalg::mat<T,3,2>(p) == linalg::mat<T,3,2>({a[0],a[1],a[2]}, {a[3],a[4],a[5]}));
        CHECK(linalg::mat<T,3,3>(p) == linalg::mat<T,3,3>({a[0],a[1],a[2]}, {a[3],a[4],a[5]}, {a[6],a[7],a[8]}));
        CHECK(linalg::mat<T,3,4>(p) == linalg::mat<T,3,4>({a[0],a[1],a[2]}, {a[3],a[4],a[5]}, {a[6],a[7],a[8]}, {a[9],a[10],a[11]}));

        CHECK(linalg::mat<T,4,2>(p) == linalg::mat<T,4,2>({a[0],a[1],a[2],a[3]}, {a[4],a[5],a[6],a[7]}));
        CHECK(linalg::mat<T,4,3>(p) == linalg::mat<T,4,3>({a[0],a[1],a[2],a[3]}, {a[4],a[5],a[6],a[7]}, {a[8],a[9],a[10],a[11]}));
        CHECK(linalg::mat<T,4,4>(p) == linalg::mat<T,4,4>({a[0],a[1],a[2],a[3]}, {a[4],a[5],a[6],a[7]}, {a[8],a[9],a[10],a[11]}, {a[12],a[13],a[14],a[15]}));
    }
}

// TODO: Appending constructors

//////////////////////////////////////////
// Test semantics of operator overloads //
//////////////////////////////////////////

TEST_CASE_TEMPLATE("arithmetic unary operator overloads on vec<T,M> are defined elementwise", T, ARITHMETIC_TYPES) 
{
    using U = decltype(+T());
    random_number_generator rng;
    for(int i=0; i<reps; ++i)
    {
        const T a=rng, b=rng, c=rng, d=rng;
        
        CHECK(+linalg::vec<T,2>(a,b    ) == linalg::vec<U,2>(+a, +b        ));
        CHECK(+linalg::vec<T,3>(a,b,c  ) == linalg::vec<U,3>(+a, +b, +c    ));
        CHECK(+linalg::vec<T,4>(a,b,c,d) == linalg::vec<U,4>(+a, +b, +c, +d));

        // NOTE: Will likely generate a warning about operator- applied to unsigned type. This is probably desirable.
        CHECK(-linalg::vec<T,2>(a,b    ) == linalg::vec<U,2>(-a, -b        ));
        CHECK(-linalg::vec<T,3>(a,b,c  ) == linalg::vec<U,3>(-a, -b, -c    ));
        CHECK(-linalg::vec<T,4>(a,b,c,d) == linalg::vec<U,4>(-a, -b, -c, -d));

        CHECK(!linalg::vec<T,2>(a,b    ) == linalg::vec<bool,2>(!a, !b        ));
        CHECK(!linalg::vec<T,3>(a,b,c  ) == linalg::vec<bool,3>(!a, !b, !c    ));
        CHECK(!linalg::vec<T,4>(a,b,c,d) == linalg::vec<bool,4>(!a, !b, !c, !d));
    }
}

TEST_CASE_TEMPLATE("arithmetic binary operator overloads on vec<T,M> are defined elementwise", T, ARITHMETIC_TYPES) 
{
    using U = decltype(T()+T());
    random_number_generator rng;
    for(int i=0; i<reps; ++i)
    {
        const T a=rng, b=rng, c=rng, d=rng, e=rng, f=rng, g=rng, h=rng;

        CHECK(linalg::vec<T,2>(a,b    ) + linalg::vec<T,2>(e,f    ) == linalg::vec<U,2>(a+e, b+f          ));
        CHECK(linalg::vec<T,3>(a,b,c  ) + linalg::vec<T,3>(e,f,g  ) == linalg::vec<U,3>(a+e, b+f, c+g     ));
        CHECK(linalg::vec<T,4>(a,b,c,d) + linalg::vec<T,4>(e,f,g,h) == linalg::vec<U,4>(a+e, b+f, c+g, d+h));

        CHECK(linalg::vec<T,2>(a,b    ) - linalg::vec<T,2>(e,f    ) == linalg::vec<U,2>(a-e, b-f          ));
        CHECK(linalg::vec<T,3>(a,b,c  ) - linalg::vec<T,3>(e,f,g  ) == linalg::vec<U,3>(a-e, b-f, c-g     ));
        CHECK(linalg::vec<T,4>(a,b,c,d) - linalg::vec<T,4>(e,f,g,h) == linalg::vec<U,4>(a-e, b-f, c-g, d-h));

        CHECK(linalg::vec<T,2>(a,b    ) * linalg::vec<T,2>(e,f    ) == linalg::vec<U,2>(a*e, b*f          ));
        CHECK(linalg::vec<T,3>(a,b,c  ) * linalg::vec<T,3>(e,f,g  ) == linalg::vec<U,3>(a*e, b*f, c*g     ));
        CHECK(linalg::vec<T,4>(a,b,c,d) * linalg::vec<T,4>(e,f,g,h) == linalg::vec<U,4>(a*e, b*f, c*g, d*h));

        CHECK(linalg::vec<T,2>(a,b    ) / linalg::vec<T,2>(e,f    ) == linalg::vec<U,2>(a/e, b/f          ));
        CHECK(linalg::vec<T,3>(a,b,c  ) / linalg::vec<T,3>(e,f,g  ) == linalg::vec<U,3>(a/e, b/f, c/g     ));
        CHECK(linalg::vec<T,4>(a,b,c,d) / linalg::vec<T,4>(e,f,g,h) == linalg::vec<U,4>(a/e, b/f, c/g, d/h));
    }
}

TEST_CASE_TEMPLATE("integral unary operator overloads on vec<T,M> are defined elementwise", T, INTEGRAL_TYPES) 
{
    using U = decltype(+T());
    random_number_generator rng;
    for(int i=0; i<reps; ++i)
    {
        const T a=rng, b=rng, c=rng, d=rng;

        CHECK(~linalg::vec<T,2>(a,b    ) == linalg::vec<U,2>(~a, ~b        ));
        CHECK(~linalg::vec<T,3>(a,b,c  ) == linalg::vec<U,3>(~a, ~b, ~c    ));
        CHECK(~linalg::vec<T,4>(a,b,c,d) == linalg::vec<U,4>(~a, ~b, ~c, ~d));
    }
}

TEST_CASE_TEMPLATE("integral binary operator overloads on vec<T,M> are defined elementwise", T, INTEGRAL_TYPES) 
{
    using U = decltype(T()+T());
    random_number_generator rng;
    for(int i=0; i<reps; ++i)
    {
        const T a=rng, b=rng, c=rng, d=rng, e=rng, f=rng, g=rng, h=rng;       

        CHECK( linalg::vec<T,2>(a,b    ) % linalg::vec<T,2>(e,f    )  == linalg::vec<U,2>(a%e, b%f          ));
        CHECK( linalg::vec<T,3>(a,b,c  ) % linalg::vec<T,3>(e,f,g  )  == linalg::vec<U,3>(a%e, b%f, c%g     ));
        CHECK( linalg::vec<T,4>(a,b,c,d) % linalg::vec<T,4>(e,f,g,h)  == linalg::vec<U,4>(a%e, b%f, c%g, d%h));

        CHECK((linalg::vec<T,2>(a,b    ) | linalg::vec<T,2>(e,f    )) == linalg::vec<U,2>(a|e, b|f          ));
        CHECK((linalg::vec<T,3>(a,b,c  ) | linalg::vec<T,3>(e,f,g  )) == linalg::vec<U,3>(a|e, b|f, c|g     ));
        CHECK((linalg::vec<T,4>(a,b,c,d) | linalg::vec<T,4>(e,f,g,h)) == linalg::vec<U,4>(a|e, b|f, c|g, d|h));

        CHECK((linalg::vec<T,2>(a,b    ) & linalg::vec<T,2>(e,f    )) == linalg::vec<U,2>(a&e, b&f          ));
        CHECK((linalg::vec<T,3>(a,b,c  ) & linalg::vec<T,3>(e,f,g  )) == linalg::vec<U,3>(a&e, b&f, c&g     ));
        CHECK((linalg::vec<T,4>(a,b,c,d) & linalg::vec<T,4>(e,f,g,h)) == linalg::vec<U,4>(a&e, b&f, c&g, d&h));

        CHECK((linalg::vec<T,2>(a,b    ) ^ linalg::vec<T,2>(e,f    )) == linalg::vec<U,2>(a^e, b^f          ));
        CHECK((linalg::vec<T,3>(a,b,c  ) ^ linalg::vec<T,3>(e,f,g  )) == linalg::vec<U,3>(a^e, b^f, c^g     ));
        CHECK((linalg::vec<T,4>(a,b,c,d) ^ linalg::vec<T,4>(e,f,g,h)) == linalg::vec<U,4>(a^e, b^f, c^g, d^h));

        CHECK((linalg::vec<T,2>(a,b    ) << linalg::vec<T,2>(e,f    )) == linalg::vec<U,2>(a<<e, b<<f            ));
        CHECK((linalg::vec<T,3>(a,b,c  ) << linalg::vec<T,3>(e,f,g  )) == linalg::vec<U,3>(a<<e, b<<f, c<<g      ));
        CHECK((linalg::vec<T,4>(a,b,c,d) << linalg::vec<T,4>(e,f,g,h)) == linalg::vec<U,4>(a<<e, b<<f, c<<g, d<<h));

        CHECK((linalg::vec<T,2>(a,b    ) >> linalg::vec<T,2>(e,f    )) == linalg::vec<U,2>(a>>e, b>>f            ));
        CHECK((linalg::vec<T,3>(a,b,c  ) >> linalg::vec<T,3>(e,f,g  )) == linalg::vec<U,3>(a>>e, b>>f, c>>g      ));
        CHECK((linalg::vec<T,4>(a,b,c,d) >> linalg::vec<T,4>(e,f,g,h)) == linalg::vec<U,4>(a>>e, b>>f, c>>g, d>>h));
    }
}

TEST_CASE_TEMPLATE("vec<T,M> does not have unintended argument dependent lookup on operator +=", T, ARITHMETIC_TYPES) 
{
    std::vector<linalg::vec<T,3>> a, b = {{0,1,2}, {0,2,3}, {0,3,4}};
    CHECK(a.size() == 0);
    CHECK(b.size() == 3);
    a = std::move(b); // This line is known to cause problems if linalg::operator+= is allowed to match too broadly.
    CHECK(a.size() == 3);
    CHECK(b.size() == 0);
}

/////////////////////

template<class T, int M> void require_approx_equal(const linalg::vec<T,M> & a, const linalg::vec<T,M> & b) { for(int j=0; j<M; ++j) REQUIRE( a[j] == doctest::Approx(b[j]) ); }
template<class T, int M, int N> void require_approx_equal(const linalg::mat<T,M,N> & a, const linalg::mat<T,M,N> & b) { for(int i=0; i<N; ++i) require_approx_equal(a[i], b[i]); }

TEST_CASE( "elementwise comparison functions produce correct results" )
{
    REQUIRE( equal  (float3(1,2,3), float3(4,-2,3)) == bool3(false, false, true ) );
    REQUIRE( nequal (float3(1,2,3), float3(4,-2,3)) == bool3(true,  true,  false) );
    REQUIRE( less   (float3(1,2,3), float3(4,-2,3)) == bool3(true,  false, false) );
    REQUIRE( greater(float3(1,2,3), float3(4,-2,3)) == bool3(false, true,  false) );
    REQUIRE( lequal (float3(1,2,3), float3(4,-2,3)) == bool3(true,  false, true ) );
    REQUIRE( gequal (float3(1,2,3), float3(4,-2,3)) == bool3(false, true,  true ) );
}

TEST_CASE( "no unintended ADL on operator +=" )
{
    std::vector<int3> tris_a = {{0,1,2}, {0,2,3}, {0,3,4}}, tris_b;
    tris_b = std::move(tris_a); // This line is known to cause problems if linalg::operator+= is allowed to match too broadly.
    REQUIRE( tris_b.size() == 3 );
    REQUIRE( tris_a.size() == 0 );
}

TEST_CASE( "fold functions behave as intended" )
{
    REQUIRE( any(bool1(false)) == false );
    REQUIRE( any(bool1(true)) == true );
    REQUIRE( any(bool3(false,false,false)) == false );
    REQUIRE( any(bool3(true,false,false)) == true );
    REQUIRE( any(bool3(false,true,false)) == true );
    REQUIRE( any(bool3(false,false,true)) == true );
    REQUIRE( all(bool2x2(bool2(true,true),bool2(true,true))) == true );
    REQUIRE( all(bool2x2(bool2(false,true),bool2(true,true))) == false );
    REQUIRE( all(bool2x2(bool2(true,false),bool2(true,true))) == false );
    REQUIRE( all(bool2x2(bool2(true,true),bool2(false,true))) == false );
    REQUIRE( all(bool2x2(bool2(true,true),bool2(true,false))) == false );
    REQUIRE( sum(int2(2,3)) == 5 );
    REQUIRE( sum(float3(2,3,4.1f)) == 9.1f );
    REQUIRE( sum(double4(2,3,4.1,5.2)) == 14.3 );
    REQUIRE( sum(double4x4(double4(1,2,3,4),double4(5,6,7,8),double4(9,10,11,12),double4(13,14,15,16))) == 136 );
}

TEST_CASE( "unary functions behave as intended" )
{
    // Unary functions should apply elementwise to their arguments
    REQUIRE( abs  (float3(1.1f, -2.3f, 3.5f)) == float3(std::abs  (1.1f), std::abs  (-2.3f), std::abs  (3.5f)) );
    REQUIRE( floor(float3(1.1f, -2.3f, 3.5f)) == float3(std::floor(1.1f), std::floor(-2.3f), std::floor(3.5f)) );
    REQUIRE( ceil (float3(1.1f, -2.3f, 3.5f)) == float3(std::ceil (1.1f), std::ceil (-2.3f), std::ceil (3.5f)) );
    REQUIRE( exp  (float3(1.1f, -2.3f, 3.5f)) == float3(std::exp  (1.1f), std::exp  (-2.3f), std::exp  (3.5f)) );
    REQUIRE( log  (float3(1.1f, +2.3f, 3.5f)) == float3(std::log  (1.1f), std::log  (+2.3f), std::log  (3.5f)) );
    REQUIRE( log10(float3(1.1f, +2.3f, 3.5f)) == float3(std::log10(1.1f), std::log10(+2.3f), std::log10(3.5f)) );
    REQUIRE( sqrt (float3(1.1f, +2.3f, 3.5f)) == float3(std::sqrt (1.1f), std::sqrt (+2.3f), std::sqrt (3.5f)) );
    REQUIRE( sin  (float3(1.1f, -2.3f, 3.5f)) == float3(std::sin  (1.1f), std::sin  (-2.3f), std::sin  (3.5f)) );
    REQUIRE( cos  (float3(1.1f, -2.3f, 3.5f)) == float3(std::cos  (1.1f), std::cos  (-2.3f), std::cos  (3.5f)) );
    REQUIRE( tan  (float3(1.1f, -2.3f, 3.5f)) == float3(std::tan  (1.1f), std::tan  (-2.3f), std::tan  (3.5f)) );
    REQUIRE( asin (float3(0.3f, -0.6f, 0.9f)) == float3(std::asin (0.3f), std::asin (-0.6f), std::asin (0.9f)) );
    REQUIRE( acos (float3(0.3f, -0.6f, 0.9f)) == float3(std::acos (0.3f), std::acos (-0.6f), std::acos (0.9f)) );
    REQUIRE( atan (float3(1.1f, -2.3f, 3.5f)) == float3(std::atan (1.1f), std::atan (-2.3f), std::atan (3.5f)) );
    REQUIRE( sinh (float3(1.1f, -2.3f, 3.5f)) == float3(std::sinh (1.1f), std::sinh (-2.3f), std::sinh (3.5f)) );
    REQUIRE( cosh (float3(1.1f, -2.3f, 3.5f)) == float3(std::cosh (1.1f), std::cosh (-2.3f), std::cosh (3.5f)) );
    REQUIRE( tanh (float3(1.1f, -2.3f, 3.5f)) == float3(std::tanh (1.1f), std::tanh (-2.3f), std::tanh (3.5f)) );
    REQUIRE( round(float3(1.1f, -2.3f, 3.5f)) == float3(std::round(1.1f), std::round(-2.3f), std::round(3.5f)) );

    // Unary functions should retain element type and vector/matrix size
    REQUIRE( abs(int4(-5)) == int4(5) );
    REQUIRE( floor(float2(7.7f)) == float2(std::floor(7.7f)) );
    REQUIRE( ceil (float3(7.7f)) == float3(std::ceil (7.7f)) );
    REQUIRE( exp  (float4(7.7f)) == float4(std::exp  (7.7f)) );
    REQUIRE( log  (double2(7.7)) == double2(std::log  (7.7)) );
    REQUIRE( log10(double3(7.7)) == double3(std::log10(7.7)) );
    REQUIRE( sqrt (double4(7.7)) == double4(std::sqrt (7.7)) );
    REQUIRE( sin  (float2x2(7.7f)) == float2x2(std::sin (7.7f)) );
    REQUIRE( cos  (float2x3(7.7f)) == float2x3(std::cos (7.7f)) );
    REQUIRE( tan  (float2x4(7.7f)) == float2x4(std::tan (7.7f)) );
    REQUIRE( asin (float3x2(0.5f)) == float3x2(std::asin(0.5f)) );
    REQUIRE( acos (float3x3(0.5f)) == float3x3(std::acos(0.5f)) );
    REQUIRE( atan (float3x4(7.7f)) == float3x4(std::atan(7.7f)) );
    REQUIRE( sinh (float4x2(7.7f)) == float4x2(std::sinh(7.7f)) );
    REQUIRE( cosh (float4x3(7.7f)) == float4x3(std::cosh(7.7f)) );
    REQUIRE( tanh (float4x4(7.7f)) == float4x4(std::tanh(7.7f)) );
}

TEST_CASE( "relational operators model LessThanComparable" )
{
    // Should not compare equal to begin with
    std::vector<int3> points = {{3,5,2}, {1,2,6}, {3,2,2}, {8,2,5}, {4,5,8}, {3,5,8}, {1,2,2}, {4,2,9}};
    std::vector<int3> ordered_points = {{1,2,2}, {1,2,6}, {3,2,2}, {3,5,2}, {3,5,8}, {4,2,9}, {4,5,8}, {8,2,5}};
    REQUIRE( points != ordered_points );

    // After sorting, points should be in ascending order, and match ordered_points
    std::sort(begin(points), end(points));
    REQUIRE( points == ordered_points );

    // Sorting in descending order should produce the reverse ordering
    std::sort(begin(points), end(points), [](const int3 & a, const int3 & b) { return a > b; });
    std::reverse(begin(ordered_points), end(ordered_points));
    REQUIRE( points == ordered_points );
}

TEST_CASE( "matrix multiplication produces correct result dimensions" )
{
    REQUIRE( mul(float2x2(), float2()) == float2() );
    REQUIRE( mul(float2x3(), float3()) == float2() );
    REQUIRE( mul(float2x4(), float4()) == float2() );
    REQUIRE( mul(float3x2(), float2()) == float3() );
    REQUIRE( mul(float3x3(), float3()) == float3() );
    REQUIRE( mul(float3x4(), float4()) == float3() );
    REQUIRE( mul(float4x2(), float2()) == float4() );
    REQUIRE( mul(float4x3(), float3()) == float4() );
    REQUIRE( mul(float4x4(), float4()) == float4() );

    REQUIRE( mul(float2x2(), float2x2()) == float2x2() );
    REQUIRE( mul(float2x3(), float3x2()) == float2x2() );
    REQUIRE( mul(float2x4(), float4x2()) == float2x2() );
    REQUIRE( mul(float2x2(), float2x3()) == float2x3() );
    REQUIRE( mul(float2x3(), float3x3()) == float2x3() );
    REQUIRE( mul(float2x4(), float4x3()) == float2x3() );
    REQUIRE( mul(float2x2(), float2x4()) == float2x4() );
    REQUIRE( mul(float2x3(), float3x4()) == float2x4() );
    REQUIRE( mul(float2x4(), float4x4()) == float2x4() );
    REQUIRE( mul(float3x2(), float2x2()) == float3x2() );
    REQUIRE( mul(float3x3(), float3x2()) == float3x2() );
    REQUIRE( mul(float3x4(), float4x2()) == float3x2() );
    REQUIRE( mul(float3x2(), float2x3()) == float3x3() );
    REQUIRE( mul(float3x3(), float3x3()) == float3x3() );
    REQUIRE( mul(float3x4(), float4x3()) == float3x3() );
    REQUIRE( mul(float3x2(), float2x4()) == float3x4() );
    REQUIRE( mul(float3x3(), float3x4()) == float3x4() );
    REQUIRE( mul(float3x4(), float4x4()) == float3x4() );
    REQUIRE( mul(float4x2(), float2x2()) == float4x2() );
    REQUIRE( mul(float4x3(), float3x2()) == float4x2() );
    REQUIRE( mul(float4x4(), float4x2()) == float4x2() );
    REQUIRE( mul(float4x2(), float2x3()) == float4x3() );
    REQUIRE( mul(float4x3(), float3x3()) == float4x3() );
    REQUIRE( mul(float4x4(), float4x3()) == float4x3() );
    REQUIRE( mul(float4x2(), float2x4()) == float4x4() );
    REQUIRE( mul(float4x3(), float3x4()) == float4x4() );
    REQUIRE( mul(float4x4(), float4x4()) == float4x4() );

    // Outer product of vec<T,M> and vec<T,N> is equivalent to product of Mx1 and 1xN matrices
    REQUIRE( outerprod(float2(), float2()) == float2x2() );
    REQUIRE( outerprod(float2(), float3()) == float2x3() );
    REQUIRE( outerprod(float2(), float4()) == float2x4() );
    REQUIRE( outerprod(float3(), float2()) == float3x2() );
    REQUIRE( outerprod(float3(), float3()) == float3x3() );
    REQUIRE( outerprod(float3(), float4()) == float3x4() );
    REQUIRE( outerprod(float4(), float2()) == float4x2() );
    REQUIRE( outerprod(float4(), float3()) == float4x3() );
    REQUIRE( outerprod(float4(), float4()) == float4x4() );
}

TEST_CASE( "matrix inverse is correct for trivial cases" )
{
    const float2x2 id2 {{1,0},{0,1}};
    const float3x3 id3 {{1,0,0},{0,1,0},{0,0,1}};
    const float4x4 id4 {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
    REQUIRE( diagonal(id2) == float2(1,1) );
    REQUIRE( diagonal(id3) == float3(1,1,1) );
    REQUIRE( diagonal(id4) == float4(1,1,1,1) );
    REQUIRE( transpose(id2) == id2 );
    REQUIRE( transpose(id3) == id3 );
    REQUIRE( transpose(id4) == id4 );
    REQUIRE( inverse(id2) == id2 );
    REQUIRE( inverse(id3) == id3 );
    REQUIRE( inverse(id4) == id4 );
    REQUIRE( adjugate(id2) == id2 );
    REQUIRE( adjugate(id3) == id3 );
    REQUIRE( adjugate(id4) == id4 );
    REQUIRE( determinant(id2) == 1.0f );
    REQUIRE( determinant(id3) == 1.0f );
    REQUIRE( determinant(id4) == 1.0f );
}

TEST_CASE_TEMPLATE( "matrix inverse is correct for general case", T, FLOATING_POINT_TYPES )
{
    const linalg::mat<T,4,4> mat {{1,2,3,4}, {5,-6,7,8}, {9,10,-11,12}, {13,14,15,-16}};
    const linalg::mat<T,4,4> inv = inverse(mat);
    const linalg::mat<T,4,4> id = mul(mat, inv);
    for(int j=0; j<4; ++j)
    {
        for(int i=0; i<4; ++i)
        {
            if(i == j) REQUIRE( id[j][i] == doctest::Approx(1.0f) );
            else REQUIRE( id[j][i] == doctest::Approx(0.0f) );
        }
    }
}

TEST_CASE_TEMPLATE( "linalg::identity functions correctly", T, ARITHMETIC_TYPES )
{
    const linalg::mat<T,2,2> a2 {linalg::identity}, b2 {{1,0},{0,1}}, c2 {};
    const linalg::mat<T,3,3> a3 {linalg::identity}, b3 {{1,0,0},{0,1,0},{0,0,1}}, c3 {};
    const linalg::mat<T,4,4> a4 {linalg::identity}, b4 {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}}, c4 {};

    REQUIRE(a2 == b2);
    REQUIRE(a3 == b3);
    REQUIRE(a4 == b4);
    REQUIRE(a2 != c2);
    REQUIRE(a3 != c3);
    REQUIRE(a4 != c4);
}

TEST_CASE_TEMPLATE( "rotation quaternions roundtrip with rotation matrices", T, FLOATING_POINT_TYPES )
{
    std::mt19937 engine;
    std::normal_distribution<T> dist;

    for(int i=0; i<1000; ++i)
    {
        linalg::vec<T,4> q = normalize(linalg::vec<T,4>(dist(engine), dist(engine), dist(engine), dist(engine)));
        linalg::vec<T,4> q2 = rotation_quat(qmat(q));
        if(dot(q, q2) > 0) // q2 should either equal q or -q
        {
            REQUIRE( std::abs(q.x - q2.x) < 0.0001 );
            REQUIRE( std::abs(q.y - q2.y) < 0.0001 );
            REQUIRE( std::abs(q.z - q2.z) < 0.0001 );
            REQUIRE( std::abs(q.w - q2.w) < 0.0001 );
        }
        else
        {
            REQUIRE( std::abs(q.x + q2.x) < 0.0001 );
            REQUIRE( std::abs(q.y + q2.y) < 0.0001 );
            REQUIRE( std::abs(q.z + q2.z) < 0.0001 );
            REQUIRE( std::abs(q.w + q2.w) < 0.0001 );        
        }
    }
}

TEST_CASE( "90 degree rotation matrices round trip with rotation quaternions" )
{
    const float3x3 matrices[] {
        {{+1,0,0},{0,+1,0},{0,0,+1}},
        {{+1,0,0},{0,-1,0},{0,0,-1}},
        {{+1,0,0},{0,0,+1},{0,-1,0}},
        {{+1,0,0},{0,0,-1},{0,+1,0}},
        {{-1,0,0},{0,+1,0},{0,0,-1}},
        {{-1,0,0},{0,-1,0},{0,0,+1}},
        {{-1,0,0},{0,0,+1},{0,+1,0}},
        {{-1,0,0},{0,0,-1},{0,-1,0}},
        {{0,+1,0},{+1,0,0},{0,0,-1}},
        {{0,+1,0},{-1,0,0},{0,0,+1}},
        {{0,+1,0},{0,0,+1},{+1,0,0}},
        {{0,+1,0},{0,0,-1},{-1,0,0}},
        {{0,-1,0},{+1,0,0},{0,0,+1}},
        {{0,-1,0},{-1,0,0},{0,0,-1}},
        {{0,-1,0},{0,0,+1},{-1,0,0}},
        {{0,-1,0},{0,0,-1},{+1,0,0}},
        {{0,0,+1},{+1,0,0},{0,+1,0}},
        {{0,0,+1},{-1,0,0},{0,-1,0}},
        {{0,0,+1},{0,+1,0},{-1,0,0}},
        {{0,0,+1},{0,-1,0},{+1,0,0}},
        {{0,0,-1},{+1,0,0},{0,-1,0}},
        {{0,0,-1},{-1,0,0},{0,+1,0}},
        {{0,0,-1},{0,+1,0},{+1,0,0}},
        {{0,0,-1},{0,-1,0},{-1,0,0}},
    };
    for(auto & m : matrices)
    {
        require_approx_equal(m, qmat(rotation_quat(m)));
    }
}

TEST_CASE( "hashing works as expected" )
{
    // std::hash specializations should take their specified type and return size_t
    REQUIRE( typeid(std::hash<int2     >()(int2     {})) == typeid(size_t) );
    REQUIRE( typeid(std::hash<float3   >()(float3   {})) == typeid(size_t) );
    REQUIRE( typeid(std::hash<double4  >()(double4  {})) == typeid(size_t) );
    REQUIRE( typeid(std::hash<int2x4   >()(int2x4   {})) == typeid(size_t) );
    REQUIRE( typeid(std::hash<float3x2 >()(float3x2 {})) == typeid(size_t) );
    REQUIRE( typeid(std::hash<double4x3>()(double4x3{})) == typeid(size_t) );

    // Small list of items which are known to have no duplicate hashes
    const float3 points[] = {
        {1,2,3}, {1,3,2}, {2,1,3}, {2,3,1}, {3,1,2}, {3,2,1},
        {1,1,1}, {2,2,2}, {3,3,3}, {4,4,4}, {5,5,5}, {6,6,6},
        {0,0,0}, {0,0,1}, {0,1,0}, {1,0,0}, {0,1,1}, {1,1,0}
    };
    const std::hash<float3> h = {};
    for(auto & a : points)
    {
        for(auto & b : points)
        {
            if(a == b) REQUIRE( h(a) == h(b) );
            else REQUIRE( h(a) != h(b) ); // Note: Not required to be different in the general case
        }
    }
}

TEST_CASE( "special quaternion functions behave as expected" )
{
    // qexp of a scalar should simply be the exp of that scalar
    require_approx_equal( qexp(float4(0,0,0,5)), float4(0,0,0,std::exp(5.0f)) );

    // e^(tau*i) == 1
    require_approx_equal( qexp(float4(6.28318531f,0,0,0)), float4(0,0,0,1) );

    // e^(tau*j) == 1
    require_approx_equal( qexp(float4(0,6.28318531f,0,0)), float4(0,0,0,1) );

    // qlog of a scalar should simply be the log of that scalar
    require_approx_equal( qlog(float4(0,0,0,5)), float4(0,0,0,std::log(5.0f)) );

    // qexp(qlog(q)) == q
    require_approx_equal( qexp(qlog(float4(1,2,3,4))), float4(1,2,3,4) );

    // qpow of a scalar should simply be the pow of that scalar
    require_approx_equal( qpow(float4(0,0,0,5), 3.14f), float4(0,0,0,std::pow(5.0f, 3.14f)) );

    // qpow(q,2) == qmul(q,q)
    require_approx_equal( qpow(float4(1,2,3,4), 2.0f), qmul(float4(1,2,3,4), float4(1,2,3,4)) );

    // qpow(q,3) == qmul(q,q,q)
    require_approx_equal( qpow(float4(1,2,3,4), 3.0f), qmul(float4(1,2,3,4), float4(1,2,3,4), float4(1,2,3,4)) );

    // qpow(qpow(q,2),3) == qpow(q,2*3)
    require_approx_equal( qpow(qpow(float4(1,2,3,4), 2.0f), 3.0f), qpow(float4(1,2,3,4), 2.0f*3.0f) );
}

float3 transform_point(const float4x4 & m, const float3 & p) { const auto r = mul(m,float4(p,1)); return r.xyz()/r.w; }

TEST_CASE( "Projection matrices behave as intended" )
{
    const float n = 0.1f, f = 10.0f;
    const float nx0 = -0.9f*n, ny0 = -0.6f*n, nx1 = 0.8f*n, ny1 = 0.7f*n, ncx = (nx0+nx1)/2, ncy = (ny0+ny1)/2;
    const float fx0 = -0.9f*f, fy0 = -0.6f*f, fx1 = 0.8f*f, fy1 = 0.7f*f, fcx = (fx0+fx1)/2, fcy = (fy0+fy1)/2;

    // Right handed OpenGL convention, x-right, y-up, z-back
    const float4x4 gl_rh = frustum_matrix(nx0, nx1, ny0, ny1, n, f, linalg::neg_z, linalg::neg_one_to_one); 
    require_approx_equal( transform_point(gl_rh, float3(ncx, ncy, -n)), float3( 0,  0, -1) );
    require_approx_equal( transform_point(gl_rh, float3(ncx, ny0, -n)), float3( 0, -1, -1) );
    require_approx_equal( transform_point(gl_rh, float3(ncx, ny1, -n)), float3( 0, +1, -1) );
    require_approx_equal( transform_point(gl_rh, float3(nx0, ncy, -n)), float3(-1,  0, -1) );
    require_approx_equal( transform_point(gl_rh, float3(nx1, ncy, -n)), float3(+1,  0, -1) );
    require_approx_equal( transform_point(gl_rh, float3(fcx, fcy, -f)), float3( 0,  0, +1) );
    require_approx_equal( transform_point(gl_rh, float3(fcx, fy0, -f)), float3( 0, -1, +1) );
    require_approx_equal( transform_point(gl_rh, float3(fcx, fy1, -f)), float3( 0, +1, +1) );
    require_approx_equal( transform_point(gl_rh, float3(fx0, fcy, -f)), float3(-1,  0, +1) );
    require_approx_equal( transform_point(gl_rh, float3(fx1, fcy, -f)), float3(+1,  0, +1) );

    // Left handed OpenGL convention, x-right, y-up, z-forward
    const float4x4 gl_lh = frustum_matrix(nx0, nx1, ny0, ny1, n, f, linalg::pos_z, linalg::neg_one_to_one);
    require_approx_equal( transform_point(gl_lh, float3(ncx, ncy, +n)), float3( 0,  0, -1) );
    require_approx_equal( transform_point(gl_lh, float3(ncx, ny0, +n)), float3( 0, -1, -1) );
    require_approx_equal( transform_point(gl_lh, float3(ncx, ny1, +n)), float3( 0, +1, -1) );
    require_approx_equal( transform_point(gl_lh, float3(nx0, ncy, +n)), float3(-1,  0, -1) );
    require_approx_equal( transform_point(gl_lh, float3(nx1, ncy, +n)), float3(+1,  0, -1) );
    require_approx_equal( transform_point(gl_lh, float3(fcx, fcy, +f)), float3( 0,  0, +1) );
    require_approx_equal( transform_point(gl_lh, float3(fcx, fy0, +f)), float3( 0, -1, +1) );
    require_approx_equal( transform_point(gl_lh, float3(fcx, fy1, +f)), float3( 0, +1, +1) );
    require_approx_equal( transform_point(gl_lh, float3(fx0, fcy, +f)), float3(-1,  0, +1) );
    require_approx_equal( transform_point(gl_lh, float3(fx1, fcy, +f)), float3(+1,  0, +1) );

    // Right handed Vulkan convention, x-right, y-down, z-forward
    const float4x4 vk_rh = frustum_matrix(nx0, nx1, ny0, ny1, n, f, linalg::pos_z, linalg::zero_to_one);
    require_approx_equal( transform_point(vk_rh, float3(ncx, ncy, +n)), float3( 0,  0, 0) );
    require_approx_equal( transform_point(vk_rh, float3(ncx, ny0, +n)), float3( 0, -1, 0) );
    require_approx_equal( transform_point(vk_rh, float3(ncx, ny1, +n)), float3( 0, +1, 0) );
    require_approx_equal( transform_point(vk_rh, float3(nx0, ncy, +n)), float3(-1,  0, 0) );
    require_approx_equal( transform_point(vk_rh, float3(nx1, ncy, +n)), float3(+1,  0, 0) );
    require_approx_equal( transform_point(vk_rh, float3(fcx, fcy, +f)), float3( 0,  0, 1) );
    require_approx_equal( transform_point(vk_rh, float3(fcx, fy0, +f)), float3( 0, -1, 1) );
    require_approx_equal( transform_point(vk_rh, float3(fcx, fy1, +f)), float3( 0, +1, 1) );
    require_approx_equal( transform_point(vk_rh, float3(fx0, fcy, +f)), float3(-1,  0, 1) );
    require_approx_equal( transform_point(vk_rh, float3(fx1, fcy, +f)), float3(+1,  0, 1) );

    // Left handed Vulkan convention, x-right, y-down, z-back
    const float4x4 vk_lh = frustum_matrix(nx0, nx1, ny0, ny1, n, f, linalg::neg_z, linalg::zero_to_one); 
    require_approx_equal( transform_point(vk_lh, float3(ncx, ncy, -n)), float3( 0,  0, 0) );
    require_approx_equal( transform_point(vk_lh, float3(ncx, ny0, -n)), float3( 0, -1, 0) );
    require_approx_equal( transform_point(vk_lh, float3(ncx, ny1, -n)), float3( 0, +1, 0) );
    require_approx_equal( transform_point(vk_lh, float3(nx0, ncy, -n)), float3(-1,  0, 0) );
    require_approx_equal( transform_point(vk_lh, float3(nx1, ncy, -n)), float3(+1,  0, 0) );
    require_approx_equal( transform_point(vk_lh, float3(fcx, fcy, -f)), float3( 0,  0, 1) );
    require_approx_equal( transform_point(vk_lh, float3(fcx, fy0, -f)), float3( 0, -1, 1) );
    require_approx_equal( transform_point(vk_lh, float3(fcx, fy1, -f)), float3( 0, +1, 1) );
    require_approx_equal( transform_point(vk_lh, float3(fx0, fcy, -f)), float3(-1,  0, 1) );
    require_approx_equal( transform_point(vk_lh, float3(fx1, fcy, -f)), float3(+1,  0, 1) );
}

template<class T> void take(const T &) {}
#define MATCH(TYPE, ...) static_assert(std::is_same<TYPE, decltype(__VA_ARGS__)>::value, #TYPE " != " #__VA_ARGS__); take(__VA_ARGS__)

TEST_CASE( "templates instantiate correctly" ) 
{
    // Declare some variables to test functions requiring an lvalue
    const float2 cf2; const float3 cf3; const float4 cf4;
    float2 f2; float3 f3; float4 f4; int2 i2; int3 i3; int4 i4;
    float fs[] = {1,2,3,4};

    // Exercise vec<T,2>
    MATCH(float2, float2());
    MATCH(float2, float2(1,2) );
    MATCH(float2, float2(5) );
    MATCH(float2, float2(fs) );
    MATCH(float2, float2(int2(3,4)) );
    MATCH(const float&, cf2[1] );
    MATCH(float&, f2[1] );

    // Exercise vec<T,3>
    MATCH(float3, float3() );
    MATCH(float3, float3(1,2,3) );
    MATCH(float3, float3(float2(),4) );
    MATCH(float3, float3(5) );
    MATCH(float3, float3(fs) );
    MATCH(float3, float3(int3(3,4,5)) );
    MATCH(const float&, cf3[1] );
    MATCH(float&, f3[1] );
    MATCH(float2 &, f3.xy() );

    // Exercise vec<T,4>
    MATCH(float4, float4() );
    MATCH(float4, float4(1,2,3,4) );
    MATCH(float4, float4(float3(),4) );
    MATCH(float4, float4(5) );
    MATCH(float4, float4(fs) );
    MATCH(float4, float4(int4(3,4,5,6)) );
    MATCH(const float&, cf4[1] );
    MATCH(float&, f4[1] );
    MATCH(float3 &, f4.xyz() );

    // TODO: Exercise mat<T,M,N> for N=2,3,4

    // Exercise sequence functions
    for(float & f : f4) take(f);
    for(float f : float4()) take(f);
    for(float4 & f : float4x4()) take(f);

    // Exercise relational operators
    MATCH(bool, int2() == int2() );
    MATCH(bool, float3() == float3() );
    MATCH(bool, double4() == double4() );
    MATCH(bool, short2() != short2() );
    MATCH(bool, int2() < int2() );
    MATCH(bool, float3() < float3() );
    MATCH(bool, double4() < double4() );
    MATCH(bool, int2() > int2() );
    MATCH(bool, float3() <= float3() );
    MATCH(bool, double4() >= double4() );

    // Exercise binary operators
    MATCH(float2 , float2 () +  float2 () );
    MATCH(float3 , float3 () -  float3 () );
    MATCH(double4, double4() *  double4() );
    MATCH(double2, double2() /  double2() );
    MATCH(int3   , int3   () %  int3   (1) );
    MATCH(int4   , int4   () |  int4   () );
    MATCH(int2   , short2 () ^  short2 () );
    MATCH(int3   , short3 () &  short3 () );
    MATCH(int3   , int3   () << int3   () );
    MATCH(int4   , int4   () >> int4   () );

    MATCH(float2 , float2 () +  float   () );
    MATCH(float3 , float3 () -  float   () );
    MATCH(double4, double4() *  double  () );
    MATCH(double2, double2() /  double  () );
    MATCH(int3   , int3   () %  int     (1) );
    MATCH(int4   , int4   () |  int     () );
    MATCH(int2   , short2 () ^  short   () );
    MATCH(int3   , short3 () &  short   () );
    MATCH(int3   , int3   () << int     () );
    MATCH(uint4  , uint4  () >> unsigned() );

    MATCH(float2 , float   () +  float2 () );
    MATCH(float3 , float   () -  float3 () );
    MATCH(double4, double  () *  double4() );
    MATCH(double2, double  () /  double2() );
    MATCH(int3   , int     () %  int3   (1) );
    MATCH(int4   , int     () |  int4   () );
    MATCH(int2   , short   () ^  short2 () );
    MATCH(int3   , short   () &  short3 () );
    MATCH(int3   , int     () << int3   () );
    MATCH(uint4  , unsigned() >> uint4  () );

    MATCH(float2&, f2 += float2() );
    MATCH(float3&, f3 -= float3() );
    MATCH(float4&, f4 *= float4() );
    MATCH(float2&, f2 /= float2() );
    MATCH(int2&  , i2 %= int2(1) );
    MATCH(int3&  , i3 |= int3() );
    MATCH(int4&  , i4 ^= int4() );
    MATCH(int2&  , i2 &= int2() );
    MATCH(int3&  , i3 <<= int3() );
    MATCH(int4&  , i4 >>= int4() );

    MATCH(float2&, f2 += float() );
    MATCH(float3&, f3 -= float() );
    MATCH(float4&, f4 *= float() );
    MATCH(float2&, f2 /= float() );
    MATCH(int2&  , i2 %= int(1) );
    MATCH(int3&  , i3 |= int() );
    MATCH(int4&  , i4 ^= int() );
    MATCH(int2&  , i2 &= int() );
    MATCH(int3&  , i3 <<= int() );
    MATCH(int4&  , i4 >>= int() );

    MATCH(float2 , min  (float2(), float2()) );
    MATCH(float3 , max  (float3(), float3()) );
    MATCH(float2 , fmod (float2(), float2()) );
    MATCH(float3 , pow  (float3(), float3()) );
    MATCH(float4 , atan2(float4(), float4()) );
    MATCH(float3 , clamp(float3(), float3(), float3()) );

    // Exercise componentwise relational ops
    MATCH(bool2, equal  (float2(), float2()) );
    MATCH(bool3, nequal (float3(), float3()) );
    MATCH(bool4, less   (double4(), double4()) );
    MATCH(bool2, greater(double2(), double2()) );
    MATCH(bool3, lequal (int3(), int3()) );
    MATCH(bool4, gequal (int4(), int4()) );

    // Exercise reduction functions
    MATCH(bool, any(bool3()) );
    MATCH(bool, all(bool4()) );
    MATCH(int, sum(int2()) );
    MATCH(float, product(float4()) );

    // Exercise selection functions
    MATCH(int, argmin(float2()) );
    MATCH(int, argmax(float3()) );
    MATCH(float, minelem(float4()) );
    MATCH(float, maxelem(float2()) );

    // Exercise vector algebra functions
    MATCH(float, cross(float2(), float2()) );
    MATCH(float3, cross(float3(), float3()) );
    MATCH(float, dot(float4(), float4()) );
    MATCH(float, length2(float2()) );
    MATCH(float, length(float3()) );
    MATCH(float, distance2(float4(), float4()) );
    MATCH(float, distance(float2(), float2()) );
    MATCH(float3, normalize(float3()) );
    MATCH(float4, lerp(float4(), float4(), float()) );
    MATCH(float4, nlerp(float4(), float4(), float()) );
    MATCH(float4, slerp(float4(), float4(), float()) );

    // Exercise quaternion algebra functions
    MATCH(float4, qconj(float4()) );
    MATCH(float4, qinv(float4()) );
    MATCH(float4, qmul(float4(), float4()) );
    MATCH(float3, qxdir(float4()) );
    MATCH(float3, qydir(float4()) );
    MATCH(float3, qzdir(float4()) );
    MATCH(float3, qrot(float4(), float3()) );
    MATCH(float , qangle(float4()) );
    MATCH(float3, qaxis(float4()) );
    MATCH(float4, qnlerp(float4(), float4(), float()) );
    MATCH(float4, qslerp(float4(), float4(), float()) );

    // Exercise factory functions
    MATCH(float4, rotation_quat(float3(), float()) );
    MATCH(float4x4, translation_matrix(float3()) );
    MATCH(float4x4, rotation_matrix(float4()) );
    MATCH(float4x4, scaling_matrix(float3()) );
    MATCH(float4x4, pose_matrix(float4(), float3()) );
    MATCH(float4x4, linalg::frustum_matrix(float(), float(), float(), float(), float(), float()) );
    MATCH(float4x4, linalg::perspective_matrix(float(), float(), float(), float()) );
}