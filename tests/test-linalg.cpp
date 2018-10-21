#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "test-linalg.h"

#include <algorithm>
#include <typeinfo>

//////////////////////////////////////////////////////////
// Test semantics of vec<T,M> element-wise constructors //
//////////////////////////////////////////////////////////

TEST_CASE_TEMPLATE("vec<T,2> can be constructed from 2 elements of type T", T, double, float, int, short, unsigned int, unsigned short) 
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

TEST_CASE_TEMPLATE("vec<T,3> can be constructed from 3 elements of type T", T, double, float, int, short, unsigned int, unsigned short) 
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

TEST_CASE_TEMPLATE("vec<T,4> can be constructed from 4 elements of type T", T, double, float, int, short, unsigned int, unsigned short)
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

TEST_CASE_TEMPLATE("mat<T,2,2> can be constructed from 2 columns of type vec<T,2>", T, double, float, int, short, unsigned int, unsigned short)
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
        CHECK(m[0][0] == m00); CHECK(m[0][0] == m00); CHECK(m.row(0)[0] == m00);
        CHECK(m[0][1] == m10); CHECK(m[0][1] == m10); CHECK(m.row(1)[0] == m10);
        CHECK(m[1][0] == m01); CHECK(m[1][0] == m01); CHECK(m.row(0)[1] == m01);
        CHECK(m[1][1] == m11); CHECK(m[1][1] == m11); CHECK(m.row(1)[1] == m11);
    }
}

TEST_CASE_TEMPLATE("mat<T,2,3> can be constructed from 3 columns of type vec<T,2>", T, double, float, int, short, unsigned int, unsigned short)
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
        CHECK(m[0][0] == m00); CHECK(m[0][0] == m00); CHECK(m.row(0)[0] == m00);
        CHECK(m[0][1] == m10); CHECK(m[0][1] == m10); CHECK(m.row(1)[0] == m10);
        CHECK(m[1][0] == m01); CHECK(m[1][0] == m01); CHECK(m.row(0)[1] == m01);
        CHECK(m[1][1] == m11); CHECK(m[1][1] == m11); CHECK(m.row(1)[1] == m11);
        CHECK(m[2][0] == m02); CHECK(m[2][0] == m02); CHECK(m.row(0)[2] == m02);
        CHECK(m[2][1] == m12); CHECK(m[2][1] == m12); CHECK(m.row(1)[2] == m12);
    }
}

TEST_CASE_TEMPLATE("mat<T,2,4> can be constructed from 4 columns of type vec<T,2>", T, double, float, int, short, unsigned int, unsigned short)
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
        CHECK(m[0][0] == m00); CHECK(m[0][0] == m00); CHECK(m.row(0)[0] == m00);
        CHECK(m[0][1] == m10); CHECK(m[0][1] == m10); CHECK(m.row(1)[0] == m10);
        CHECK(m[1][0] == m01); CHECK(m[1][0] == m01); CHECK(m.row(0)[1] == m01);
        CHECK(m[1][1] == m11); CHECK(m[1][1] == m11); CHECK(m.row(1)[1] == m11);
        CHECK(m[2][0] == m02); CHECK(m[2][0] == m02); CHECK(m.row(0)[2] == m02);
        CHECK(m[2][1] == m12); CHECK(m[2][1] == m12); CHECK(m.row(1)[2] == m12);
        CHECK(m[3][0] == m03); CHECK(m[3][0] == m03); CHECK(m.row(0)[3] == m03);
        CHECK(m[3][1] == m13); CHECK(m[3][1] == m13); CHECK(m.row(1)[3] == m13);
    }
}

TEST_CASE_TEMPLATE("mat<T,3,2> can be constructed from 2 columns of type vec<T,3>", T, double, float, int, short, unsigned int, unsigned short)
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
        CHECK(m[0][0] == m00); CHECK(m[0][0] == m00); CHECK(m.row(0)[0] == m00);
        CHECK(m[0][1] == m10); CHECK(m[0][1] == m10); CHECK(m.row(1)[0] == m10);
        CHECK(m[0][2] == m20); CHECK(m[0][2] == m20); CHECK(m.row(2)[0] == m20);
        CHECK(m[1][0] == m01); CHECK(m[1][0] == m01); CHECK(m.row(0)[1] == m01);
        CHECK(m[1][1] == m11); CHECK(m[1][1] == m11); CHECK(m.row(1)[1] == m11);
        CHECK(m[1][2] == m21); CHECK(m[1][2] == m21); CHECK(m.row(2)[1] == m21);
    }
}

TEST_CASE_TEMPLATE("mat<T,3,3> can be constructed from 3 columns of type vec<T,3>", T, double, float, int, short, unsigned int, unsigned short)
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
        CHECK(m[0][0] == m00); CHECK(m[0][0] == m00); CHECK(m.row(0)[0] == m00);
        CHECK(m[0][1] == m10); CHECK(m[0][1] == m10); CHECK(m.row(1)[0] == m10);
        CHECK(m[0][2] == m20); CHECK(m[0][2] == m20); CHECK(m.row(2)[0] == m20);
        CHECK(m[1][0] == m01); CHECK(m[1][0] == m01); CHECK(m.row(0)[1] == m01);
        CHECK(m[1][1] == m11); CHECK(m[1][1] == m11); CHECK(m.row(1)[1] == m11);
        CHECK(m[1][2] == m21); CHECK(m[1][2] == m21); CHECK(m.row(2)[1] == m21);
        CHECK(m[2][0] == m02); CHECK(m[2][0] == m02); CHECK(m.row(0)[2] == m02);
        CHECK(m[2][1] == m12); CHECK(m[2][1] == m12); CHECK(m.row(1)[2] == m12);
        CHECK(m[2][2] == m22); CHECK(m[2][2] == m22); CHECK(m.row(2)[2] == m22);
    }
}

TEST_CASE_TEMPLATE("mat<T,3,4> can be constructed from 4 columns of type vec<T,3>", T, double, float, int, short, unsigned int, unsigned short)
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
        CHECK(m[0][0] == m00); CHECK(m[0][0] == m00); CHECK(m.row(0)[0] == m00);
        CHECK(m[0][1] == m10); CHECK(m[0][1] == m10); CHECK(m.row(1)[0] == m10);
        CHECK(m[0][2] == m20); CHECK(m[0][2] == m20); CHECK(m.row(2)[0] == m20);
        CHECK(m[1][0] == m01); CHECK(m[1][0] == m01); CHECK(m.row(0)[1] == m01);
        CHECK(m[1][1] == m11); CHECK(m[1][1] == m11); CHECK(m.row(1)[1] == m11);
        CHECK(m[1][2] == m21); CHECK(m[1][2] == m21); CHECK(m.row(2)[1] == m21);
        CHECK(m[2][0] == m02); CHECK(m[2][0] == m02); CHECK(m.row(0)[2] == m02);
        CHECK(m[2][1] == m12); CHECK(m[2][1] == m12); CHECK(m.row(1)[2] == m12);
        CHECK(m[2][2] == m22); CHECK(m[2][2] == m22); CHECK(m.row(2)[2] == m22);
        CHECK(m[3][0] == m03); CHECK(m[3][0] == m03); CHECK(m.row(0)[3] == m03);
        CHECK(m[3][1] == m13); CHECK(m[3][1] == m13); CHECK(m.row(1)[3] == m13);
        CHECK(m[3][2] == m23); CHECK(m[3][2] == m23); CHECK(m.row(2)[3] == m23);
    }
}

TEST_CASE_TEMPLATE("mat<T,4,2> can be constructed from 2 columns of type vec<T,4>", T, double, float, int, short, unsigned int, unsigned short)
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
        CHECK(m[0][0] == m00); CHECK(m[0][0] == m00); CHECK(m.row(0)[0] == m00);
        CHECK(m[0][1] == m10); CHECK(m[0][1] == m10); CHECK(m.row(1)[0] == m10);
        CHECK(m[0][2] == m20); CHECK(m[0][2] == m20); CHECK(m.row(2)[0] == m20);
        CHECK(m[0][3] == m30); CHECK(m[0][3] == m30); CHECK(m.row(3)[0] == m30);
        CHECK(m[1][0] == m01); CHECK(m[1][0] == m01); CHECK(m.row(0)[1] == m01);
        CHECK(m[1][1] == m11); CHECK(m[1][1] == m11); CHECK(m.row(1)[1] == m11);
        CHECK(m[1][2] == m21); CHECK(m[1][2] == m21); CHECK(m.row(2)[1] == m21);
        CHECK(m[1][3] == m31); CHECK(m[1][3] == m31); CHECK(m.row(3)[1] == m31);
    }
}

TEST_CASE_TEMPLATE("mat<T,4,3> can be constructed from 3 columns of type vec<T,4>", T, double, float, int, short, unsigned int, unsigned short)
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
        CHECK(m[0][0] == m00); CHECK(m[0][0] == m00); CHECK(m.row(0)[0] == m00);
        CHECK(m[0][1] == m10); CHECK(m[0][1] == m10); CHECK(m.row(1)[0] == m10);
        CHECK(m[0][2] == m20); CHECK(m[0][2] == m20); CHECK(m.row(2)[0] == m20);
        CHECK(m[0][3] == m30); CHECK(m[0][3] == m30); CHECK(m.row(3)[0] == m30);
        CHECK(m[1][0] == m01); CHECK(m[1][0] == m01); CHECK(m.row(0)[1] == m01);
        CHECK(m[1][1] == m11); CHECK(m[1][1] == m11); CHECK(m.row(1)[1] == m11);
        CHECK(m[1][2] == m21); CHECK(m[1][2] == m21); CHECK(m.row(2)[1] == m21);
        CHECK(m[1][3] == m31); CHECK(m[1][3] == m31); CHECK(m.row(3)[1] == m31);
        CHECK(m[2][0] == m02); CHECK(m[2][0] == m02); CHECK(m.row(0)[2] == m02);
        CHECK(m[2][1] == m12); CHECK(m[2][1] == m12); CHECK(m.row(1)[2] == m12);
        CHECK(m[2][2] == m22); CHECK(m[2][2] == m22); CHECK(m.row(2)[2] == m22);
        CHECK(m[2][3] == m32); CHECK(m[2][3] == m32); CHECK(m.row(3)[2] == m32);
    }
}

TEST_CASE_TEMPLATE("mat<T,4,4> can be constructed from 4 columns of type vec<T,4>", T, double, float, int, short, unsigned int, unsigned short)
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
        CHECK(m[0][0] == m00); CHECK(m[0][0] == m00); CHECK(m.row(0)[0] == m00);
        CHECK(m[0][1] == m10); CHECK(m[0][1] == m10); CHECK(m.row(1)[0] == m10);
        CHECK(m[0][2] == m20); CHECK(m[0][2] == m20); CHECK(m.row(2)[0] == m20);
        CHECK(m[0][3] == m30); CHECK(m[0][3] == m30); CHECK(m.row(3)[0] == m30);
        CHECK(m[1][0] == m01); CHECK(m[1][0] == m01); CHECK(m.row(0)[1] == m01);
        CHECK(m[1][1] == m11); CHECK(m[1][1] == m11); CHECK(m.row(1)[1] == m11);
        CHECK(m[1][2] == m21); CHECK(m[1][2] == m21); CHECK(m.row(2)[1] == m21);
        CHECK(m[1][3] == m31); CHECK(m[1][3] == m31); CHECK(m.row(3)[1] == m31);
        CHECK(m[2][0] == m02); CHECK(m[2][0] == m02); CHECK(m.row(0)[2] == m02);
        CHECK(m[2][1] == m12); CHECK(m[2][1] == m12); CHECK(m.row(1)[2] == m12);
        CHECK(m[2][2] == m22); CHECK(m[2][2] == m22); CHECK(m.row(2)[2] == m22);
        CHECK(m[2][3] == m32); CHECK(m[2][3] == m32); CHECK(m.row(3)[2] == m32);
        CHECK(m[3][0] == m03); CHECK(m[3][0] == m03); CHECK(m.row(0)[3] == m03);
        CHECK(m[3][1] == m13); CHECK(m[3][1] == m13); CHECK(m.row(1)[3] == m13);
        CHECK(m[3][2] == m23); CHECK(m[3][2] == m23); CHECK(m.row(2)[3] == m23);
        CHECK(m[3][3] == m33); CHECK(m[3][3] == m33); CHECK(m.row(3)[3] == m33);
    }
}

///////////////////////////////////////////////////
// Test semantics of operator == and operator != //
///////////////////////////////////////////////////

TEST_CASE_TEMPLATE("vec<T,M> a and b compare equal if they contain exactly the same elements", T, double, float, int, short, unsigned int, unsigned short)
{
    random_number_generator rng;
    for(int i=0; i<reps; ++i)
    {
        const T a=rng, b=rng, c=rng, d=rng;

        CHECK(linalg::vec<T,2>(a,b    ) == linalg::vec<T,2>(a,b    ));
        CHECK(linalg::vec<T,3>(a,b,c  ) == linalg::vec<T,3>(a,b,c  ));
        CHECK(linalg::vec<T,4>(a,b,c,d) == linalg::vec<T,4>(a,b,c,d));

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

TEST_CASE_TEMPLATE("vec<T,M> a and b compare unequal if at least one element differs between them", T, double, float, int, short, unsigned int, unsigned short)
{
    random_number_generator rng;
    for(int i=0; i<reps; ++i)
    {
        const T a=rng, b=rng, c=rng, d=rng;

        CHECK_FALSE(linalg::vec<T,2>(a,b    ) != linalg::vec<T,2>(a,b    ));
        CHECK_FALSE(linalg::vec<T,3>(a,b,c  ) != linalg::vec<T,3>(a,b,c  ));
        CHECK_FALSE(linalg::vec<T,4>(a,b,c,d) != linalg::vec<T,4>(a,b,c,d));

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

TEST_CASE_TEMPLATE("vec<T,M>'s default constructor zero-initializes its elements", T, double, float, int, short, unsigned int, unsigned short) 
{
    const linalg::vec<T,2> v2; CHECK(v2 == linalg::vec<T,2>(0,0));
    const linalg::vec<T,3> v3; CHECK(v3 == linalg::vec<T,3>(0,0,0));
    const linalg::vec<T,4> v4; CHECK(v4 == linalg::vec<T,4>(0,0,0,0));
}


TEST_CASE_TEMPLATE("mat<T,M,N>'s default constructor zero-initializes its columns", T, double, float, int, short, unsigned int, unsigned short) 
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

TEST_CASE_TEMPLATE("vec<T,M>'s scalar constructor initializes its elements to the specified scalar", T, double, float, int, short, unsigned int, unsigned short) 
{
    random_number_generator rng;
    for(int i=0; i<reps; ++i)
    {
        const T s = rng;
        const linalg::vec<T,2> v2(s); CHECK(v2 == linalg::vec<T,2>(s,s));
        const linalg::vec<T,3> v3(s); CHECK(v3 == linalg::vec<T,3>(s,s,s));
        const linalg::vec<T,4> v4(s); CHECK(v4 == linalg::vec<T,4>(s,s,s,s));
    }
}

TEST_CASE_TEMPLATE("mat<T,M,N>'s scalar constructor initializes its columns to the specified scalar", T, double, float, int, short, unsigned int, unsigned short) 
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

//////////////////////////////////////////
// Test semantics of operator overloads //
//////////////////////////////////////////

TEST_CASE_TEMPLATE("arithmetic unary operator overloads on vec<T,M> are defined elementwise", T, double, float, int, short, unsigned int, unsigned short) 
{
    random_number_generator rng;
    for(int i=0; i<reps; ++i)
    {
        const T a=rng, b=rng, c=rng, d=rng;
        
        CHECK(+linalg::vec<T,2>(a,b    ) == linalg::vec<decltype(+T()),2>(+a, +b        ));
        CHECK(+linalg::vec<T,3>(a,b,c  ) == linalg::vec<decltype(+T()),3>(+a, +b, +c    ));
        CHECK(+linalg::vec<T,4>(a,b,c,d) == linalg::vec<decltype(+T()),4>(+a, +b, +c, +d));

        CHECK(!linalg::vec<T,2>(a,b    ) == linalg::vec<bool,2>(!a, !b        ));
        CHECK(!linalg::vec<T,3>(a,b,c  ) == linalg::vec<bool,3>(!a, !b, !c    ));
        CHECK(!linalg::vec<T,4>(a,b,c,d) == linalg::vec<bool,4>(!a, !b, !c, !d));
    }
}

TEST_CASE_TEMPLATE("arithmetic unary operator - on vec<T,M> is defined elementwise", T, double, float, int, short) 
{
    random_number_generator rng;
    for(int i=0; i<reps; ++i)
    {
        const T a=rng, b=rng, c=rng, d=rng;

        CHECK(-linalg::vec<T,2>(a,b    ) == linalg::vec<decltype(-T()),2>(-a, -b        ));
        CHECK(-linalg::vec<T,3>(a,b,c  ) == linalg::vec<decltype(-T()),3>(-a, -b, -c    ));
        CHECK(-linalg::vec<T,4>(a,b,c,d) == linalg::vec<decltype(-T()),4>(-a, -b, -c, -d));
    }
}

TEST_CASE_TEMPLATE("arithmetic binary operator overloads on vec<T,M> are defined elementwise", T, double, float, int, short, unsigned int, unsigned short) 
{
    using U = decltype(T()+T());
    random_number_generator rng;
    for(int i=0; i<reps; ++i)
    {
        T a=rng, b=rng, c=rng, d=rng, e=rng, f=rng, g=rng, h=rng;

        CHECK(linalg::vec<T,2>(a,b    ) + linalg::vec<T,2>(e,f    ) == linalg::vec<U,2>(a+e, b+f          ));
        CHECK(linalg::vec<T,3>(a,b,c  ) + linalg::vec<T,3>(e,f,g  ) == linalg::vec<U,3>(a+e, b+f, c+g     ));
        CHECK(linalg::vec<T,4>(a,b,c,d) + linalg::vec<T,4>(e,f,g,h) == linalg::vec<U,4>(a+e, b+f, c+g, d+h));

        CHECK(linalg::vec<T,2>(a,b    ) - linalg::vec<T,2>(e,f    ) == linalg::vec<U,2>(a-e, b-f          ));
        CHECK(linalg::vec<T,3>(a,b,c  ) - linalg::vec<T,3>(e,f,g  ) == linalg::vec<U,3>(a-e, b-f, c-g     ));
        CHECK(linalg::vec<T,4>(a,b,c,d) - linalg::vec<T,4>(e,f,g,h) == linalg::vec<U,4>(a-e, b-f, c-g, d-h));

        CHECK(linalg::vec<T,2>(a,b    ) * linalg::vec<T,2>(e,f    ) == linalg::vec<U,2>(a*e, b*f          ));
        CHECK(linalg::vec<T,3>(a,b,c  ) * linalg::vec<T,3>(e,f,g  ) == linalg::vec<U,3>(a*e, b*f, c*g     ));
        CHECK(linalg::vec<T,4>(a,b,c,d) * linalg::vec<T,4>(e,f,g,h) == linalg::vec<U,4>(a*e, b*f, c*g, d*h));

        // Ensure nonzero denominator
        e = e ? e : 2;
        f = f ? f : 2;
        g = g ? g : 2;
        h = h ? h : 2;

        CHECK(linalg::vec<T,2>(a,b    ) / linalg::vec<T,2>(e,f    ) == linalg::vec<U,2>(a/e, b/f          ));
        CHECK(linalg::vec<T,3>(a,b,c  ) / linalg::vec<T,3>(e,f,g  ) == linalg::vec<U,3>(a/e, b/f, c/g     ));
        CHECK(linalg::vec<T,4>(a,b,c,d) / linalg::vec<T,4>(e,f,g,h) == linalg::vec<U,4>(a/e, b/f, c/g, d/h));
    }
}

TEST_CASE_TEMPLATE("integral unary operator overloads on vec<T,M> are defined elementwise", T, int, short, unsigned int, unsigned short) 
{
    using U = decltype(~T());
    random_number_generator rng;
    for(int i=0; i<reps; ++i)
    {
        const T a=rng, b=rng, c=rng, d=rng;

        CHECK(~linalg::vec<T,2>(a,b    ) == linalg::vec<U,2>(~a, ~b        ));
        CHECK(~linalg::vec<T,3>(a,b,c  ) == linalg::vec<U,3>(~a, ~b, ~c    ));
        CHECK(~linalg::vec<T,4>(a,b,c,d) == linalg::vec<U,4>(~a, ~b, ~c, ~d));
    }
}

TEST_CASE_TEMPLATE("integral binary operator overloads on vec<T,M> are defined elementwise", T, int, short, unsigned int, unsigned short) 
{
    using U = decltype(T()+T());
    random_number_generator rng;
    for(int i=0; i<reps; ++i)
    {
        T a=rng, b=rng, c=rng, d=rng, e=rng, f=rng, g=rng, h=rng;       

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

        // Ensure nonzero denominator
        e = e ? e : 2;
        f = f ? f : 2;
        g = g ? g : 2;
        h = h ? h : 2;

        CHECK( linalg::vec<T,2>(a,b    ) % linalg::vec<T,2>(e,f    )  == linalg::vec<U,2>(a%e, b%f          ));
        CHECK( linalg::vec<T,3>(a,b,c  ) % linalg::vec<T,3>(e,f,g  )  == linalg::vec<U,3>(a%e, b%f, c%g     ));
        CHECK( linalg::vec<T,4>(a,b,c,d) % linalg::vec<T,4>(e,f,g,h)  == linalg::vec<U,4>(a%e, b%f, c%g, d%h));
    }
}

TEST_CASE_TEMPLATE("elementwise comparison functions on vec<T,M> are defined elementwise", T, double, float, int, short, unsigned int, unsigned short) 
{
    random_number_generator rng;
    for(int i=0; i<reps; ++i)
    {
        const T a=rng, b=rng, c=rng, d=rng, e=rng, f=rng, g=rng, h=rng;       

        CHECK(equal  (linalg::vec<T,2>(a,b    ), linalg::vec<T,2>(e,f    )) == linalg::vec<bool,2>(a==e, b==f            ));
        CHECK(equal  (linalg::vec<T,3>(a,b,c  ), linalg::vec<T,3>(e,f,g  )) == linalg::vec<bool,3>(a==e, b==f, c==g      ));
        CHECK(equal  (linalg::vec<T,4>(a,b,c,d), linalg::vec<T,4>(e,f,g,h)) == linalg::vec<bool,4>(a==e, b==f, c==g, d==h));

        CHECK(nequal (linalg::vec<T,2>(a,b    ), linalg::vec<T,2>(e,f    )) == linalg::vec<bool,2>(a!=e, b!=f            ));
        CHECK(nequal (linalg::vec<T,3>(a,b,c  ), linalg::vec<T,3>(e,f,g  )) == linalg::vec<bool,3>(a!=e, b!=f, c!=g      ));
        CHECK(nequal (linalg::vec<T,4>(a,b,c,d), linalg::vec<T,4>(e,f,g,h)) == linalg::vec<bool,4>(a!=e, b!=f, c!=g, d!=h));

        CHECK(less   (linalg::vec<T,2>(a,b    ), linalg::vec<T,2>(e,f    )) == linalg::vec<bool,2>(a<e, b<f          ));
        CHECK(less   (linalg::vec<T,3>(a,b,c  ), linalg::vec<T,3>(e,f,g  )) == linalg::vec<bool,3>(a<e, b<f, c<g     ));
        CHECK(less   (linalg::vec<T,4>(a,b,c,d), linalg::vec<T,4>(e,f,g,h)) == linalg::vec<bool,4>(a<e, b<f, c<g, d<h));

        CHECK(greater(linalg::vec<T,2>(a,b    ), linalg::vec<T,2>(e,f    )) == linalg::vec<bool,2>(a>e, b>f          ));
        CHECK(greater(linalg::vec<T,3>(a,b,c  ), linalg::vec<T,3>(e,f,g  )) == linalg::vec<bool,3>(a>e, b>f, c>g     ));
        CHECK(greater(linalg::vec<T,4>(a,b,c,d), linalg::vec<T,4>(e,f,g,h)) == linalg::vec<bool,4>(a>e, b>f, c>g, d>h));

        CHECK(lequal (linalg::vec<T,2>(a,b    ), linalg::vec<T,2>(e,f    )) == linalg::vec<bool,2>(a<=e, b<=f            ));
        CHECK(lequal (linalg::vec<T,3>(a,b,c  ), linalg::vec<T,3>(e,f,g  )) <= linalg::vec<bool,3>(a<=e, b<=f, c<=g      ));
        CHECK(lequal (linalg::vec<T,4>(a,b,c,d), linalg::vec<T,4>(e,f,g,h)) <= linalg::vec<bool,4>(a<=e, b<=f, c<=g, d<=h));

        CHECK(gequal (linalg::vec<T,2>(a,b    ), linalg::vec<T,2>(e,f    )) == linalg::vec<bool,2>(a>=e, b>=f            ));
        CHECK(gequal (linalg::vec<T,3>(a,b,c  ), linalg::vec<T,3>(e,f,g  )) == linalg::vec<bool,3>(a>=e, b>=f, c>=g      ));
        CHECK(gequal (linalg::vec<T,4>(a,b,c,d), linalg::vec<T,4>(e,f,g,h)) == linalg::vec<bool,4>(a>=e, b>=f, c>=g, d>=h));
    }
}

TEST_CASE_TEMPLATE("vec<T,M> does not have unintended argument dependent lookup on operator +=", T, double, float, int, short, unsigned int, unsigned short) 
{
    std::vector<linalg::vec<T,3>> a, b = {{0,1,2}, {0,2,3}, {0,3,4}};
    CHECK(a.size() == 0);
    CHECK(b.size() == 3);
    a = std::move(b); // This line is known to cause problems if linalg::operator+= is allowed to match too broadly.
    CHECK(a.size() == 3);
    CHECK(b.size() == 0);
}

/////////////////////

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
    std::vector<int3> tris_a = {{0,1,2}, {0,2,3}, {0,3,4}}, tris_b; // This line is known to cause problems if linalg::operator!= is allowed to match too broadly
    tris_b = std::move(tris_a); // This line is known to cause problems if linalg::operator+= is allowed to match too broadly.
    REQUIRE( tris_b.size() == 3 );
    REQUIRE( tris_a.size() == 0 );
}

TEST_CASE( "fold functions behave as intended" )
{
    REQUIRE( any(bool3(false,false,false)) == false );
    REQUIRE( any(bool3(true,false,false)) == true );
    REQUIRE( any(bool3(false,true,false)) == true );
    REQUIRE( any(bool3(false,false,true)) == true );
    REQUIRE( sum(int2(2,3)) == 5 );
    REQUIRE( sum(float3(2,3,4.1f)) == 9.1f );
    REQUIRE( sum(double4(2,3,4.1,5.2)) == 14.3 );
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
    REQUIRE( floor(float2 (7.7f)) == float2 (std::floor(7.7f)) );
    REQUIRE( ceil (float3 (7.7f)) == float3 (std::ceil (7.7f)) );
    REQUIRE( exp  (float4 (7.7f)) == float4 (std::exp  (7.7f)) );
    REQUIRE( log  (double2(7.7 )) == double2(std::log  (7.7 )) );
    REQUIRE( log10(double3(7.7 )) == double3(std::log10(7.7 )) );
    REQUIRE( sqrt (double4(7.7 )) == double4(std::sqrt (7.7 )) );
    REQUIRE( sin  (float2 (7.7f)) == float2 (std::sin  (7.7f)) );
    REQUIRE( cos  (float3 (7.7f)) == float3 (std::cos  (7.7f)) );
    REQUIRE( tan  (float4 (7.7f)) == float4 (std::tan  (7.7f)) );
    REQUIRE( asin (double2(0.5 )) == double2(std::asin (0.5 )) );
    REQUIRE( acos (double3(0.5 )) == double3(std::acos (0.5 )) );
    REQUIRE( atan (double4(7.7 )) == double4(std::atan (7.7 )) );
    REQUIRE( sinh (float2 (7.7f)) == float2 (std::sinh (7.7f)) );
    REQUIRE( cosh (float3 (7.7f)) == float3 (std::cosh (7.7f)) );
    REQUIRE( tanh (float4 (7.7f)) == float4 (std::tanh (7.7f)) );
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

TEST_CASE_TEMPLATE( "3D cross products behave as intended", T, float, double )
{
    random_number_generator rng;
    for(int i=0; i<reps; ++i)
    {
        const vec3<T> a = rng, b = rng;
        check_approx_equal( cross(a,b), -cross(b,a) );
        CHECK( dot(cross(a,b),a) == doctest::Approx(0) );
        CHECK( dot(cross(a,b),b) == doctest::Approx(0) );
        if(angle(a,b) > 0) CHECK( length2(cross(a,b)) > 0 );
    }
}

TEST_CASE_TEMPLATE( "2D cross products behave as intended", T, int, float, double )
{
    random_number_generator rng;
    for(int i=0; i<reps; ++i)
    {
        vec2<T> a = rng, b = rng;
        T c = rng;
        CHECK( vec3<T>{0,0,cross(a,b)} == cross(vec3<T>{a,0}, vec3<T>{b,0}) );
        CHECK( vec3<T>{0,0,cross(b,a)} == cross(vec3<T>{b,0}, vec3<T>{a,0}) );
        CHECK( vec3<T>{cross(a,c),0} == cross(vec3<T>{a,0}, vec3<T>{0,0,c}) );
        CHECK( vec3<T>{cross(b,c),0} == cross(vec3<T>{b,0}, vec3<T>{0,0,c}) );
        CHECK( vec3<T>{cross(c,a),0} == cross(vec3<T>{0,0,c}, vec3<T>{a,0}) );
        CHECK( vec3<T>{cross(c,b),0} == cross(vec3<T>{0,0,c}, vec3<T>{b,0}) );
    }
}

template<class T> void take(const T &) {}
#define MATCH(TYPE, ...) static_assert(std::is_same<TYPE, decltype(__VA_ARGS__)>::value, #TYPE " != " #__VA_ARGS__); take(__VA_ARGS__)

TEST_CASE( "templates instantiate correctly" ) 
{
    // Declare some variables to test functions requiring an lvalue
    const float2 cf2; const float3 cf3; const float4 cf4;
    float2 f2; float3 f3; float4 f4; int2 i2; int3 i3; int4 i4;

    // Exercise vec<T,2>
    MATCH(float2, float2());
    MATCH(float2, float2(1,2) );
    MATCH(float2, float2(5) );
    MATCH(float2, float2(int2(3,4)) );
    MATCH(const float&, cf2[1] );
    MATCH(float&, f2[1] );

    // Exercise vec<T,3>
    MATCH(float3, float3() );
    MATCH(float3, float3(1,2,3) );
    MATCH(float3, float3(float2(),4) );
    MATCH(float3, float3(5) );
    MATCH(float3, float3(int3(3,4,5)) );
    MATCH(const float&, cf3[1] );
    MATCH(float&, f3[1] );

    // Exercise vec<T,4>
    MATCH(float4, float4() );
    MATCH(float4, float4(1,2,3,4) );
    MATCH(float4, float4(float3(),4) );
    MATCH(float4, float4(5) );
    MATCH(float4, float4(int4(3,4,5,6)) );
    MATCH(const float&, cf4[1] );
    MATCH(float&, f4[1] );

    // TODO: Exercise mat<T,M,N> for N=2,3,4

    // Exercise sequence functions
    for(float & f : f4) take(f);
    for(float f : float4()) take(f);
    for(float4 & f : float4x4()) take(f);

    // Exercise relational operators
    MATCH(bool, int2() == int2() );
    MATCH(bool, float3() == float3() );
    MATCH(bool, double4() == double4() );
    MATCH(bool, int2() != int2() );
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
    MATCH(int2   , int2   () ^  int2   () );
    MATCH(int3   , int3   () &  int3   () );
    MATCH(int3   , int3   () << int3   () );
    MATCH(int4   , int4   () >> int4   () );

    MATCH(float2 , float2 () +  float   () );
    MATCH(float3 , float3 () -  float   () );
    MATCH(double4, double4() *  double  () );
    MATCH(double2, double2() /  double  () );
    MATCH(int3   , int3   () %  int     (1) );
    MATCH(int4   , int4   () |  int     () );
    MATCH(int2   , int2   () ^  short   () );
    MATCH(int3   , int3   () &  short   () );
    MATCH(int3   , int3   () << int     () );
    MATCH(int4   , int4   () >> unsigned() );

    MATCH(float2 , float   () +  float2 () );
    MATCH(float3 , float   () -  float3 () );
    MATCH(double4, double  () *  double4() );
    MATCH(double2, double  () /  double2() );
    MATCH(int3   , int     () %  int3   (1) );
    MATCH(int4   , int     () |  int4   () );
    MATCH(int2   , short   () ^  int2   () );
    MATCH(int3   , short   () &  int3   () );
    MATCH(int3   , int     () << int3   () );
    MATCH(int4   , int     () >> int4   () );

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
    MATCH(quatf, conjugate(quatf()) );
    MATCH(quatf, inverse(quatf()) );
    MATCH(quatf, quatf() * quatf() );
    MATCH(float3, qxdir(quatf()) );
    MATCH(float3, qydir(quatf()) );
    MATCH(float3, qzdir(quatf()) );
    MATCH(float3, qrot(quatf(), float3()) );
    MATCH(float , qangle(quatf()) );
    MATCH(float3, qaxis(quatf()) );
}