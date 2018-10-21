#include "test-linalg.h"

        struct op_pos { template<class A> constexpr auto operator() (A a) const -> decltype(+a) { return +a; } };
        struct op_neg { template<class A> constexpr auto operator() (A a) const -> decltype(-a) { return -a; } };
        struct op_not { template<class A> constexpr auto operator() (A a) const -> decltype(!a) { return !a; } };
        struct op_cmp { template<class A> constexpr auto operator() (A a) const -> decltype(~(a)) { return ~a; } };
        struct op_mul { template<class A, class B> constexpr auto operator() (A a, B b) const -> decltype(a * b)  { return a * b; } };
        struct op_div { template<class A, class B> constexpr auto operator() (A a, B b) const -> decltype(a / b)  { return a / b; } };
        struct op_mod { template<class A, class B> constexpr auto operator() (A a, B b) const -> decltype(a % b)  { return a % b; } };
        struct op_add { template<class A, class B> constexpr auto operator() (A a, B b) const -> decltype(a + b)  { return a + b; } };
        struct op_sub { template<class A, class B> constexpr auto operator() (A a, B b) const -> decltype(a - b)  { return a - b; } };
        struct op_lsh { template<class A, class B> constexpr auto operator() (A a, B b) const -> decltype(a << b) { return a << b; } };
        struct op_rsh { template<class A, class B> constexpr auto operator() (A a, B b) const -> decltype(a >> b) { return a >> b; } };
        struct op_lt  { template<class A, class B> constexpr auto operator() (A a, B b) const -> decltype(a < b)  { return a < b; } };
        struct op_gt  { template<class A, class B> constexpr auto operator() (A a, B b) const -> decltype(a > b)  { return a > b; } };
        struct op_le  { template<class A, class B> constexpr auto operator() (A a, B b) const -> decltype(a <= b) { return a <= b; } };
        struct op_ge  { template<class A, class B> constexpr auto operator() (A a, B b) const -> decltype(a >= b) { return a >= b; } };
        struct op_eq  { template<class A, class B> constexpr auto operator() (A a, B b) const -> decltype(a == b) { return a == b; } };
        struct op_ne  { template<class A, class B> constexpr auto operator() (A a, B b) const -> decltype(a != b) { return a != b; } };
        struct op_int { template<class A, class B> constexpr auto operator() (A a, B b) const -> decltype(a & b)  { return a & b; } };        
        struct op_xor { template<class A, class B> constexpr auto operator() (A a, B b) const -> decltype(a ^ b)  { return a ^ b; } };
        struct op_un  { template<class A, class B> constexpr auto operator() (A a, B b) const -> decltype(a | b)  { return a | b; } };
        struct op_and { template<class A, class B> constexpr auto operator() (A a, B b) const -> decltype(a && b) { return a && b; } };
        struct op_or  { template<class A, class B> constexpr auto operator() (A a, B b) const -> decltype(a || b) { return a || b; } };


// SFINAE based expression validity helpers
namespace detail
{
    // These function overloads exist if specific the expression T{} op U{} is valid, and return true_type
    template<class T> inline decltype(+std::declval<T>(), std::true_type{}) has_op_pos(int) { return {}; }
    template<class T> inline decltype(-std::declval<T>(), std::true_type{}) has_op_neg(int) { return {}; }
    template<class T> inline decltype(!std::declval<T>(), std::true_type{}) has_op_not(int) { return {}; }
    template<class T> inline decltype(~std::declval<T>(), std::true_type{}) has_op_cmp(int) { return {}; }
    
    template<class T, class U> inline decltype(std::declval<T>() *  std::declval<U>(), std::true_type{}) has_op_mul(int) { return {}; }
    template<class T, class U> inline decltype(std::declval<T>() /  std::declval<U>(), std::true_type{}) has_op_div(int) { return {}; }  
    template<class T, class U> inline decltype(std::declval<T>() %  std::declval<U>(), std::true_type{}) has_op_mod(int) { return {}; }
    template<class T, class U> inline decltype(std::declval<T>() +  std::declval<U>(), std::true_type{}) has_op_add(int) { return {}; }
    template<class T, class U> inline decltype(std::declval<T>() -  std::declval<U>(), std::true_type{}) has_op_sub(int) { return {}; }
    template<class T, class U> inline decltype(std::declval<T>() << std::declval<U>(), std::true_type{}) has_op_lsh(int) { return {}; }  
    template<class T, class U> inline decltype(std::declval<T>() >> std::declval<U>(), std::true_type{}) has_op_rsh(int) { return {}; }
    template<class T, class U> inline decltype(std::declval<T>() <  std::declval<U>(), std::true_type{}) has_op_lt (int) { return {}; }
    template<class T, class U> inline decltype(std::declval<T>() >  std::declval<U>(), std::true_type{}) has_op_gt (int) { return {}; }
    template<class T, class U> inline decltype(std::declval<T>() <= std::declval<U>(), std::true_type{}) has_op_le (int) { return {}; }  
    template<class T, class U> inline decltype(std::declval<T>() >= std::declval<U>(), std::true_type{}) has_op_ge (int) { return {}; }
    template<class T, class U> inline decltype(std::declval<T>() == std::declval<U>(), std::true_type{}) has_op_eq (int) { return {}; }
    template<class T, class U> inline decltype(std::declval<T>() != std::declval<U>(), std::true_type{}) has_op_ne (int) { return {}; }
    template<class T, class U> inline decltype(std::declval<T>() &  std::declval<U>(), std::true_type{}) has_op_int(int) { return {}; }
    template<class T, class U> inline decltype(std::declval<T>() ^  std::declval<U>(), std::true_type{}) has_op_xor(int) { return {}; }
    template<class T, class U> inline decltype(std::declval<T>() |  std::declval<U>(), std::true_type{}) has_op_un (int) { return {}; }  
    template<class T, class U> inline decltype(std::declval<T>() && std::declval<U>(), std::true_type{}) has_op_and(int) { return {}; }
    template<class T, class U> inline decltype(std::declval<T>() || std::declval<U>(), std::true_type{}) has_op_or (int) { return {}; }

    template<class T, class U> inline decltype(std::declval<T>() *=  std::declval<U>(), std::true_type{}) has_op_assign_mul(int) { return {}; }
    template<class T, class U> inline decltype(std::declval<T>() /=  std::declval<U>(), std::true_type{}) has_op_assign_div(int) { return {}; }  
    template<class T, class U> inline decltype(std::declval<T>() %=  std::declval<U>(), std::true_type{}) has_op_assign_mod(int) { return {}; }
    template<class T, class U> inline decltype(std::declval<T>() +=  std::declval<U>(), std::true_type{}) has_op_assign_add(int) { return {}; }
    template<class T, class U> inline decltype(std::declval<T>() -=  std::declval<U>(), std::true_type{}) has_op_assign_sub(int) { return {}; }
    template<class T, class U> inline decltype(std::declval<T>() <<= std::declval<U>(), std::true_type{}) has_op_assign_lsh(int) { return {}; }  
    template<class T, class U> inline decltype(std::declval<T>() >>= std::declval<U>(), std::true_type{}) has_op_assign_rsh(int) { return {}; }
    template<class T, class U> inline decltype(std::declval<T>() &=  std::declval<U>(), std::true_type{}) has_op_assign_int(int) { return {}; }
    template<class T, class U> inline decltype(std::declval<T>() ^=  std::declval<U>(), std::true_type{}) has_op_assign_xor(int) { return {}; }
    template<class T, class U> inline decltype(std::declval<T>() |=  std::declval<U>(), std::true_type{}) has_op_assign_un (int) { return {}; }

    // These function overloads always exist, but have lowest selection priority, and return false_type
    template<class T> inline std::false_type has_op_pos(...) { return {}; }
    template<class T> inline std::false_type has_op_neg(...) { return {}; }
    template<class T> inline std::false_type has_op_not(...) { return {}; }
    template<class T> inline std::false_type has_op_cmp(...) { return {}; }
    
    template<class T, class U> inline std::false_type has_op_mul(...) { return {}; }
    template<class T, class U> inline std::false_type has_op_div(...) { return {}; }  
    template<class T, class U> inline std::false_type has_op_mod(...) { return {}; }
    template<class T, class U> inline std::false_type has_op_add(...) { return {}; }
    template<class T, class U> inline std::false_type has_op_sub(...) { return {}; }
    template<class T, class U> inline std::false_type has_op_lsh(...) { return {}; }  
    template<class T, class U> inline std::false_type has_op_rsh(...) { return {}; }
    template<class T, class U> inline std::false_type has_op_lt (...) { return {}; }
    template<class T, class U> inline std::false_type has_op_gt (...) { return {}; }
    template<class T, class U> inline std::false_type has_op_le (...) { return {}; }  
    template<class T, class U> inline std::false_type has_op_ge (...) { return {}; }
    template<class T, class U> inline std::false_type has_op_eq (...) { return {}; }
    template<class T, class U> inline std::false_type has_op_ne (...) { return {}; }
    template<class T, class U> inline std::false_type has_op_int(...) { return {}; }
    template<class T, class U> inline std::false_type has_op_xor(...) { return {}; }
    template<class T, class U> inline std::false_type has_op_un (...) { return {}; }  
    template<class T, class U> inline std::false_type has_op_and(...) { return {}; }
    template<class T, class U> inline std::false_type has_op_or (...) { return {}; }

    template<class T, class U> inline std::false_type has_op_assign_mul(...) { return {}; }
    template<class T, class U> inline std::false_type has_op_assign_div(...) { return {}; }
    template<class T, class U> inline std::false_type has_op_assign_mod(...) { return {}; }
    template<class T, class U> inline std::false_type has_op_assign_add(...) { return {}; }
    template<class T, class U> inline std::false_type has_op_assign_sub(...) { return {}; }
    template<class T, class U> inline std::false_type has_op_assign_lsh(...) { return {}; }
    template<class T, class U> inline std::false_type has_op_assign_rsh(...) { return {}; }
    template<class T, class U> inline std::false_type has_op_assign_int(...) { return {}; }
    template<class T, class U> inline std::false_type has_op_assign_xor(...) { return {}; }
    template<class T, class U> inline std::false_type has_op_assign_un (...) { return {}; }
}
template<class T> struct has_op_pos : decltype(detail::has_op_pos<T>(0)) {};
template<class T> struct has_op_neg : decltype(detail::has_op_neg<T>(0)) {};
template<class T> struct has_op_not : decltype(detail::has_op_not<T>(0)) {};
template<class T> struct has_op_cmp : decltype(detail::has_op_cmp<T>(0)) {};

template<class T, class U> struct has_op_mul : decltype(detail::has_op_mul<T,U>(0)) {};
template<class T, class U> struct has_op_div : decltype(detail::has_op_div<T,U>(0)) {};
template<class T, class U> struct has_op_mod : decltype(detail::has_op_mod<T,U>(0)) {};
template<class T, class U> struct has_op_add : decltype(detail::has_op_add<T,U>(0)) {};
template<class T, class U> struct has_op_sub : decltype(detail::has_op_sub<T,U>(0)) {};
template<class T, class U> struct has_op_lsh : decltype(detail::has_op_lsh<T,U>(0)) {};
template<class T, class U> struct has_op_rsh : decltype(detail::has_op_rsh<T,U>(0)) {};
template<class T, class U> struct has_op_lt  : decltype(detail::has_op_lt <T,U>(0)) {};
template<class T, class U> struct has_op_gt  : decltype(detail::has_op_gt <T,U>(0)) {};
template<class T, class U> struct has_op_le  : decltype(detail::has_op_le <T,U>(0)) {};
template<class T, class U> struct has_op_ge  : decltype(detail::has_op_ge <T,U>(0)) {};
template<class T, class U> struct has_op_eq  : decltype(detail::has_op_eq <T,U>(0)) {};
template<class T, class U> struct has_op_ne  : decltype(detail::has_op_ne <T,U>(0)) {};
template<class T, class U> struct has_op_int : decltype(detail::has_op_int<T,U>(0)) {};
template<class T, class U> struct has_op_xor : decltype(detail::has_op_xor<T,U>(0)) {};
template<class T, class U> struct has_op_un  : decltype(detail::has_op_un <T,U>(0)) {};
template<class T, class U> struct has_op_and : decltype(detail::has_op_and<T,U>(0)) {};
template<class T, class U> struct has_op_or  : decltype(detail::has_op_or <T,U>(0)) {};

template<class T, class U> struct has_op_assign_mul: decltype(detail::has_op_assign_mul<T,U>(0)) {};
template<class T, class U> struct has_op_assign_div: decltype(detail::has_op_assign_div<T,U>(0)) {};
template<class T, class U> struct has_op_assign_mod: decltype(detail::has_op_assign_mod<T,U>(0)) {};
template<class T, class U> struct has_op_assign_add: decltype(detail::has_op_assign_add<T,U>(0)) {};
template<class T, class U> struct has_op_assign_sub: decltype(detail::has_op_assign_sub<T,U>(0)) {};
template<class T, class U> struct has_op_assign_lsh: decltype(detail::has_op_assign_lsh<T,U>(0)) {};
template<class T, class U> struct has_op_assign_rsh: decltype(detail::has_op_assign_rsh<T,U>(0)) {};
template<class T, class U> struct has_op_assign_int: decltype(detail::has_op_assign_int<T,U>(0)) {};
template<class T, class U> struct has_op_assign_xor: decltype(detail::has_op_assign_xor<T,U>(0)) {};
template<class T, class U> struct has_op_assign_un : decltype(detail::has_op_assign_un <T,U>(0)) {};

template<class T, int M> void check_vec_arithmetic_operators()
{
    // Can apply unary operations to vectors
    CHECK(has_op_pos<linalg::vec<T,M>>::value); // +a
    CHECK(has_op_neg<linalg::vec<T,M>>::value); // -a
    CHECK(has_op_not<linalg::vec<T,M>>::value); // !a
    
    // Can apply binary operations with vectors of same size
    CHECK(has_op_mul<linalg::vec<T,M>, linalg::vec<T,M>>::value); // a*b
    CHECK(has_op_div<linalg::vec<T,M>, linalg::vec<T,M>>::value); // a/b
    CHECK(has_op_add<linalg::vec<T,M>, linalg::vec<T,M>>::value); // a+b
    CHECK(has_op_sub<linalg::vec<T,M>, linalg::vec<T,M>>::value); // a-b

    // Can apply binary operations with scalars on the left or on the right
    CHECK(has_op_mul<T, linalg::vec<T,M>>::value); // a*b
    CHECK(has_op_div<T, linalg::vec<T,M>>::value); // a/b
    CHECK(has_op_add<T, linalg::vec<T,M>>::value); // a+b
    CHECK(has_op_sub<T, linalg::vec<T,M>>::value); // a-b

    CHECK(has_op_mul<linalg::vec<T,M>, T>::value); // a*b
    CHECK(has_op_div<linalg::vec<T,M>, T>::value); // a/b
    CHECK(has_op_add<linalg::vec<T,M>, T>::value); // a+b
    CHECK(has_op_sub<linalg::vec<T,M>, T>::value); // a-b

    // Cannot apply binary operations with vectors of different size
    CHECK_FALSE(has_op_mul<linalg::vec<T,M+1>, linalg::vec<T,M>>::value); // a*b
    CHECK_FALSE(has_op_div<linalg::vec<T,M+1>, linalg::vec<T,M>>::value); // a/b
    CHECK_FALSE(has_op_add<linalg::vec<T,M+1>, linalg::vec<T,M>>::value); // a+b
    CHECK_FALSE(has_op_sub<linalg::vec<T,M+1>, linalg::vec<T,M>>::value); // a-b

    CHECK_FALSE(has_op_mul<linalg::vec<T,M>, linalg::vec<T,M+1>>::value); // a*b
    CHECK_FALSE(has_op_div<linalg::vec<T,M>, linalg::vec<T,M+1>>::value); // a/b
    CHECK_FALSE(has_op_add<linalg::vec<T,M>, linalg::vec<T,M+1>>::value); // a+b
    CHECK_FALSE(has_op_sub<linalg::vec<T,M>, linalg::vec<T,M+1>>::value); // a-b

    // Can apply assignment operations with vectors of same size
    CHECK(has_op_assign_mul<linalg::vec<T,M> &, linalg::vec<T,M>>::value); // a*=b
    CHECK(has_op_assign_div<linalg::vec<T,M> &, linalg::vec<T,M>>::value); // a/=b
    CHECK(has_op_assign_add<linalg::vec<T,M> &, linalg::vec<T,M>>::value); // a+=b
    CHECK(has_op_assign_sub<linalg::vec<T,M> &, linalg::vec<T,M>>::value); // a-=b

    // Can apply assignment operations with scalars
    CHECK(has_op_assign_mul<linalg::vec<T,M> &, T>::value); // a*=b 
    CHECK(has_op_assign_div<linalg::vec<T,M> &, T>::value); // a/=b 
    CHECK(has_op_assign_add<linalg::vec<T,M> &, T>::value); // a+=b 
    CHECK(has_op_assign_sub<linalg::vec<T,M> &, T>::value); // a-=b
}

TEST_CASE_TEMPLATE("Overloaded operators for vectors with floating-point element type", T, float, double) 
{
    check_vec_arithmetic_operators<T,1>();
    check_vec_arithmetic_operators<T,2>();
    check_vec_arithmetic_operators<T,3>();
    check_vec_arithmetic_operators<T,4>();
}

template<class T, int M> void check_vec_integer_operators()
{
    check_vec_arithmetic_operators<T,M>();

    // Can apply unary operations to vectors
    CHECK(has_op_cmp<linalg::vec<T,M>>::value); // ~a
    
    // Can apply binary operations with vectors of same size
    CHECK(has_op_mod<linalg::vec<T,M>, linalg::vec<T,M>>::value); // a%b
    CHECK(has_op_lsh<linalg::vec<T,M>, linalg::vec<T,M>>::value); // a<<b
    CHECK(has_op_rsh<linalg::vec<T,M>, linalg::vec<T,M>>::value); // a>>b
    CHECK(has_op_int<linalg::vec<T,M>, linalg::vec<T,M>>::value); // a&b
    CHECK(has_op_xor<linalg::vec<T,M>, linalg::vec<T,M>>::value); // a^b
    CHECK(has_op_un <linalg::vec<T,M>, linalg::vec<T,M>>::value); // a|b

    // Can apply binary operations with scalars on the left or on the right
    CHECK(has_op_mod<T, linalg::vec<T,M>>::value); // a%b
    CHECK(has_op_lsh<T, linalg::vec<T,M>>::value); // a<<b
    CHECK(has_op_rsh<T, linalg::vec<T,M>>::value); // a>>b
    CHECK(has_op_int<T, linalg::vec<T,M>>::value); // a&b
    CHECK(has_op_xor<T, linalg::vec<T,M>>::value); // a^b
    CHECK(has_op_un <T, linalg::vec<T,M>>::value); // a|b

    CHECK(has_op_mod<linalg::vec<T,M>, T>::value); // a%b
    CHECK(has_op_lsh<linalg::vec<T,M>, T>::value); // a<<b
    CHECK(has_op_rsh<linalg::vec<T,M>, T>::value); // a>>b
    CHECK(has_op_int<linalg::vec<T,M>, T>::value); // a&b
    CHECK(has_op_xor<linalg::vec<T,M>, T>::value); // a^b
    CHECK(has_op_un <linalg::vec<T,M>, T>::value); // a|b

    // Cannot apply binary operations with vectors of different size
    CHECK_FALSE(has_op_mul<linalg::vec<T,M+1>, linalg::vec<T,M>>::value); // a*b
    CHECK_FALSE(has_op_lsh<linalg::vec<T,M+1>, linalg::vec<T,M>>::value); // a<<b
    CHECK_FALSE(has_op_rsh<linalg::vec<T,M+1>, linalg::vec<T,M>>::value); // a>>b
    CHECK_FALSE(has_op_int<linalg::vec<T,M+1>, linalg::vec<T,M>>::value); // a&b
    CHECK_FALSE(has_op_xor<linalg::vec<T,M+1>, linalg::vec<T,M>>::value); // a^b
    CHECK_FALSE(has_op_un <linalg::vec<T,M+1>, linalg::vec<T,M>>::value); // a|b

    CHECK_FALSE(has_op_mul<linalg::vec<T,M>, linalg::vec<T,M+1>>::value); // a*b
    CHECK_FALSE(has_op_lsh<linalg::vec<T,M>, linalg::vec<T,M+1>>::value); // a<<b
    CHECK_FALSE(has_op_rsh<linalg::vec<T,M>, linalg::vec<T,M+1>>::value); // a>>b
    CHECK_FALSE(has_op_int<linalg::vec<T,M>, linalg::vec<T,M+1>>::value); // a&b
    CHECK_FALSE(has_op_xor<linalg::vec<T,M>, linalg::vec<T,M+1>>::value); // a^b
    CHECK_FALSE(has_op_un <linalg::vec<T,M>, linalg::vec<T,M+1>>::value); // a|b

    // Can apply assignment operations with vectors of same size
    CHECK(has_op_assign_mod<linalg::vec<T,M> &, linalg::vec<T,M>>::value); // a%=b
    CHECK(has_op_assign_lsh<linalg::vec<T,M> &, linalg::vec<T,M>>::value); // a<<=b
    CHECK(has_op_assign_rsh<linalg::vec<T,M> &, linalg::vec<T,M>>::value); // a>>=b
    CHECK(has_op_assign_int<linalg::vec<T,M> &, linalg::vec<T,M>>::value); // a&=b
    CHECK(has_op_assign_xor<linalg::vec<T,M> &, linalg::vec<T,M>>::value); // a^=b
    CHECK(has_op_assign_un <linalg::vec<T,M> &, linalg::vec<T,M>>::value); // a|=b

    // Can apply assignment operations with scalars
    CHECK(has_op_assign_mul<linalg::vec<T,M> &, T>::value); // a*=b 
    CHECK(has_op_assign_div<linalg::vec<T,M> &, T>::value); // a/=b 
    CHECK(has_op_assign_mod<linalg::vec<T,M> &, T>::value); // a%=b 
    CHECK(has_op_assign_add<linalg::vec<T,M> &, T>::value); // a+=b 
    CHECK(has_op_assign_sub<linalg::vec<T,M> &, T>::value); // a-=b 
    CHECK(has_op_assign_lsh<linalg::vec<T,M> &, T>::value); // a<<=b 
    CHECK(has_op_assign_rsh<linalg::vec<T,M> &, T>::value); // a>>=b 
    CHECK(has_op_assign_int<linalg::vec<T,M> &, T>::value); // a&=b 
    CHECK(has_op_assign_xor<linalg::vec<T,M> &, T>::value); // a^=b 
    CHECK(has_op_assign_un <linalg::vec<T,M> &, T>::value); // a|=b 
}

TEST_CASE_TEMPLATE("Overloaded operators for vectors with integral element type", T, int, unsigned) 
{
    check_vec_integer_operators<T,1>();
    check_vec_integer_operators<T,2>();
    check_vec_integer_operators<T,3>();
    check_vec_integer_operators<T,4>();
}

template<class T, int M, int N> void check_matrix_operators()
{
    // Unary + and - are defined for matrices
    CHECK(has_op_pos<linalg::mat<T,M,N>>::value); // +a
    CHECK(has_op_neg<linalg::mat<T,M,N>>::value); // -a

    // Unary ! and ~ are NOT defined for matrices
    CHECK_FALSE(has_op_not<linalg::mat<T,M,N>>::value); // !a
    CHECK_FALSE(has_op_cmp<linalg::mat<T,M,N>>::value); // ~a
    
    // Binary + and - are defined for matrices of the same size
    CHECK(has_op_add<linalg::mat<T,M,N>, linalg::mat<T,M,N>>::value); // a+b
    CHECK(has_op_sub<linalg::mat<T,M,N>, linalg::mat<T,M,N>>::value); // a-b

    // Binary * is defined for matrices of compatible size
    CHECK(has_op_mul<linalg::mat<T,M,1>, linalg::mat<T,1,N>>::value);
    CHECK(has_op_mul<linalg::mat<T,M,2>, linalg::mat<T,2,N>>::value);
    CHECK(has_op_mul<linalg::mat<T,M,3>, linalg::mat<T,3,N>>::value);
    CHECK(has_op_mul<linalg::mat<T,M,4>, linalg::mat<T,4,N>>::value);

    // Binary * is not defined for matrices of incompatible size
    CHECK_FALSE(has_op_mul<linalg::mat<T,M,1>, linalg::mat<T,2,N>>::value);
    CHECK_FALSE(has_op_mul<linalg::mat<T,M,1>, linalg::mat<T,3,N>>::value);
    CHECK_FALSE(has_op_mul<linalg::mat<T,M,1>, linalg::mat<T,4,N>>::value);
    CHECK_FALSE(has_op_mul<linalg::mat<T,M,2>, linalg::mat<T,1,N>>::value);
    CHECK_FALSE(has_op_mul<linalg::mat<T,M,2>, linalg::mat<T,3,N>>::value);
    CHECK_FALSE(has_op_mul<linalg::mat<T,M,2>, linalg::mat<T,4,N>>::value);
    CHECK_FALSE(has_op_mul<linalg::mat<T,M,3>, linalg::mat<T,1,N>>::value);
    CHECK_FALSE(has_op_mul<linalg::mat<T,M,3>, linalg::mat<T,2,N>>::value);   
    CHECK_FALSE(has_op_mul<linalg::mat<T,M,3>, linalg::mat<T,4,N>>::value);
    CHECK_FALSE(has_op_mul<linalg::mat<T,M,4>, linalg::mat<T,1,N>>::value);
    CHECK_FALSE(has_op_mul<linalg::mat<T,M,4>, linalg::mat<T,2,N>>::value);
    CHECK_FALSE(has_op_mul<linalg::mat<T,M,4>, linalg::mat<T,3,N>>::value);
    
    // Binary /, %, <<, >>, &, ^, | are never defined for matrices
    CHECK_FALSE(has_op_div<linalg::mat<T,M,N>, linalg::mat<T,M,N>>::value); // a/b
    CHECK_FALSE(has_op_mod<linalg::mat<T,M,N>, linalg::mat<T,M,N>>::value); // a%b
    CHECK_FALSE(has_op_lsh<linalg::mat<T,M,N>, linalg::mat<T,M,N>>::value); // a<<b
    CHECK_FALSE(has_op_rsh<linalg::mat<T,M,N>, linalg::mat<T,M,N>>::value); // a>>b
    CHECK_FALSE(has_op_int<linalg::mat<T,M,N>, linalg::mat<T,M,N>>::value); // a&b
    CHECK_FALSE(has_op_xor<linalg::mat<T,M,N>, linalg::mat<T,M,N>>::value); // a^b
    CHECK_FALSE(has_op_un <linalg::mat<T,M,N>, linalg::mat<T,M,N>>::value); // a|b

    // Only operator * and / are defined with a scalar on the right
    CHECK(has_op_mul<linalg::mat<T,M,N>, T>::value); // a*b
    CHECK(has_op_div<linalg::mat<T,M,N>, T>::value); // a/b
    CHECK_FALSE(has_op_add<linalg::mat<T,M,N>, T>::value); // a+b
    CHECK_FALSE(has_op_sub<linalg::mat<T,M,N>, T>::value); // a-b
    CHECK_FALSE(has_op_mod<linalg::mat<T,M,N>, T>::value); // a%b
    CHECK_FALSE(has_op_lsh<linalg::mat<T,M,N>, T>::value); // a<<b
    CHECK_FALSE(has_op_rsh<linalg::mat<T,M,N>, T>::value); // a>>b
    CHECK_FALSE(has_op_int<linalg::mat<T,M,N>, T>::value); // a&b
    CHECK_FALSE(has_op_xor<linalg::mat<T,M,N>, T>::value); // a^b
    CHECK_FALSE(has_op_un <linalg::mat<T,M,N>, T>::value); // a|b

    // Only operator * is defined with a scalar on the left
    CHECK(has_op_mul<T, linalg::mat<T,M,N>>::value); // a*b
    CHECK_FALSE(has_op_div<T, linalg::mat<T,M,N>>::value); // a/b
    CHECK_FALSE(has_op_add<T, linalg::mat<T,M,N>>::value); // a+b
    CHECK_FALSE(has_op_sub<T, linalg::mat<T,M,N>>::value); // a-b
    CHECK_FALSE(has_op_mod<T, linalg::mat<T,M,N>>::value); // a%b
    CHECK_FALSE(has_op_lsh<T, linalg::mat<T,M,N>>::value); // a<<b
    CHECK_FALSE(has_op_rsh<T, linalg::mat<T,M,N>>::value); // a>>b
    CHECK_FALSE(has_op_int<T, linalg::mat<T,M,N>>::value); // a&b
    CHECK_FALSE(has_op_xor<T, linalg::mat<T,M,N>>::value); // a^b
    CHECK_FALSE(has_op_un <T, linalg::mat<T,M,N>>::value); // a|b
    
    // Binary +=, -= is defined for matrices of same size
    CHECK(has_op_assign_add<linalg::mat<T,M,N> &, linalg::mat<T,M,N>>::value); // a+=b
    CHECK(has_op_assign_sub<linalg::mat<T,M,N> &, linalg::mat<T,M,N>>::value); // a-=b

    // Binary *= is defined for matrices of compatible size
    CHECK(has_op_assign_mul<linalg::mat<T,M,N> &, linalg::mat<T,N,N>>::value);

    // Binary *= is not defined for matrices of incompatible size
    CHECK_FALSE(has_op_assign_mul<linalg::mat<T,M,1> &, linalg::mat<T,2,N>>::value);
    CHECK_FALSE(has_op_assign_mul<linalg::mat<T,M,1> &, linalg::mat<T,3,N>>::value);
    CHECK_FALSE(has_op_assign_mul<linalg::mat<T,M,1> &, linalg::mat<T,4,N>>::value);
    CHECK_FALSE(has_op_assign_mul<linalg::mat<T,M,2> &, linalg::mat<T,1,N>>::value);
    CHECK_FALSE(has_op_assign_mul<linalg::mat<T,M,2> &, linalg::mat<T,3,N>>::value);
    CHECK_FALSE(has_op_assign_mul<linalg::mat<T,M,2> &, linalg::mat<T,4,N>>::value);
    CHECK_FALSE(has_op_assign_mul<linalg::mat<T,M,3> &, linalg::mat<T,1,N>>::value);
    CHECK_FALSE(has_op_assign_mul<linalg::mat<T,M,3> &, linalg::mat<T,2,N>>::value);   
    CHECK_FALSE(has_op_assign_mul<linalg::mat<T,M,3> &, linalg::mat<T,4,N>>::value);
    CHECK_FALSE(has_op_assign_mul<linalg::mat<T,M,4> &, linalg::mat<T,1,N>>::value);
    CHECK_FALSE(has_op_assign_mul<linalg::mat<T,M,4> &, linalg::mat<T,2,N>>::value);
    CHECK_FALSE(has_op_assign_mul<linalg::mat<T,M,4> &, linalg::mat<T,3,N>>::value);

    // Binary /=, %=, <<=, >>=, &=, ^=, |= are never defined for matrices
    CHECK_FALSE(has_op_assign_div<linalg::mat<T,M,N> &, linalg::mat<T,M,N>>::value); // a/=b
    CHECK_FALSE(has_op_assign_mod<linalg::mat<T,M,N> &, linalg::mat<T,M,N>>::value); // a%=b
    CHECK_FALSE(has_op_assign_lsh<linalg::mat<T,M,N> &, linalg::mat<T,M,N>>::value); // a<<=b
    CHECK_FALSE(has_op_assign_rsh<linalg::mat<T,M,N> &, linalg::mat<T,M,N>>::value); // a>>=b
    CHECK_FALSE(has_op_assign_int<linalg::mat<T,M,N> &, linalg::mat<T,M,N>>::value); // a&=b
    CHECK_FALSE(has_op_assign_xor<linalg::mat<T,M,N> &, linalg::mat<T,M,N>>::value); // a^=b
    CHECK_FALSE(has_op_assign_un <linalg::mat<T,M,N> &, linalg::mat<T,M,N>>::value); // a|=b

    // Only operator *= and /= are defined for scalars
    CHECK(has_op_assign_mul<linalg::mat<T,M,N> &, T>::value); // a*=b 
    CHECK(has_op_assign_div<linalg::mat<T,M,N> &, T>::value); // a/=b 
    CHECK_FALSE(has_op_assign_mod<linalg::mat<T,M,N> &, T>::value); // a%=b 
    CHECK_FALSE(has_op_assign_add<linalg::mat<T,M,N> &, T>::value); // a+=b 
    CHECK_FALSE(has_op_assign_sub<linalg::mat<T,M,N> &, T>::value); // a-=b 
    CHECK_FALSE(has_op_assign_lsh<linalg::mat<T,M,N> &, T>::value); // a<<=b 
    CHECK_FALSE(has_op_assign_rsh<linalg::mat<T,M,N> &, T>::value); // a>>=b 
    CHECK_FALSE(has_op_assign_int<linalg::mat<T,M,N> &, T>::value); // a&=b 
    CHECK_FALSE(has_op_assign_xor<linalg::mat<T,M,N> &, T>::value); // a^=b 
    CHECK_FALSE(has_op_assign_un <linalg::mat<T,M,N> &, T>::value); // a|=b 
}

TEST_CASE_TEMPLATE("Overloaded operators for matrices", T, double, float, int, unsigned int) 
{
    check_matrix_operators<T,1,1>();
    check_matrix_operators<T,1,2>();
    check_matrix_operators<T,1,3>();
    check_matrix_operators<T,1,4>();
    check_matrix_operators<T,2,1>();
    check_matrix_operators<T,2,2>();
    check_matrix_operators<T,2,3>();
    check_matrix_operators<T,2,4>();
    check_matrix_operators<T,3,1>();
    check_matrix_operators<T,3,2>();
    check_matrix_operators<T,3,3>();
    check_matrix_operators<T,3,4>();
    check_matrix_operators<T,4,1>();
    check_matrix_operators<T,4,2>();
    check_matrix_operators<T,4,3>();
    check_matrix_operators<T,4,4>();
}

TEST_CASE_TEMPLATE("Overloaded operators for quaternions", T, double, float, int, unsigned int) 
{
    // Only operator + and - are defined for unary quaternions
    CHECK(has_op_pos<linalg::quat<T>>::value); // +a
    CHECK(has_op_neg<linalg::quat<T>>::value); // -a
    CHECK_FALSE(has_op_not<linalg::quat<T>>::value); // !a
    CHECK_FALSE(has_op_cmp<linalg::quat<T>>::value); // ~a
    
    // Only operator +, -, and * are defined for quat $ quat
    CHECK(has_op_add<linalg::quat<T>, linalg::quat<T>>::value); // a+b
    CHECK(has_op_sub<linalg::quat<T>, linalg::quat<T>>::value); // a-b
    CHECK(has_op_mul<linalg::quat<T>, linalg::quat<T>>::value); // a*b
    CHECK_FALSE(has_op_div<linalg::quat<T>, linalg::quat<T>>::value); // a/b
    CHECK_FALSE(has_op_mod<linalg::quat<T>, linalg::quat<T>>::value); // a%b
    CHECK_FALSE(has_op_lsh<linalg::quat<T>, linalg::quat<T>>::value); // a<<b
    CHECK_FALSE(has_op_rsh<linalg::quat<T>, linalg::quat<T>>::value); // a>>b
    CHECK_FALSE(has_op_int<linalg::quat<T>, linalg::quat<T>>::value); // a&b
    CHECK_FALSE(has_op_xor<linalg::quat<T>, linalg::quat<T>>::value); // a^b
    CHECK_FALSE(has_op_un <linalg::quat<T>, linalg::quat<T>>::value); // a|b

    // Only operator * and / are defined for quat $ scalar
    CHECK(has_op_mul<linalg::quat<T>, T>::value); // a*b
    CHECK(has_op_div<linalg::quat<T>, T>::value); // a/b
    CHECK_FALSE(has_op_add<linalg::quat<T>, T>::value); // a+b
    CHECK_FALSE(has_op_sub<linalg::quat<T>, T>::value); // a-b
    CHECK_FALSE(has_op_mod<linalg::quat<T>, T>::value); // a%b
    CHECK_FALSE(has_op_lsh<linalg::quat<T>, T>::value); // a<<b
    CHECK_FALSE(has_op_rsh<linalg::quat<T>, T>::value); // a>>b
    CHECK_FALSE(has_op_int<linalg::quat<T>, T>::value); // a&b
    CHECK_FALSE(has_op_xor<linalg::quat<T>, T>::value); // a^b
    CHECK_FALSE(has_op_un <linalg::quat<T>, T>::value); // a|b

    // Only operator * is defined for scalar $ quat
    CHECK(has_op_mul<T, linalg::quat<T>>::value); // a*b
    CHECK_FALSE(has_op_div<T, linalg::quat<T>>::value); // a/b
    CHECK_FALSE(has_op_add<T, linalg::quat<T>>::value); // a+b
    CHECK_FALSE(has_op_sub<T, linalg::quat<T>>::value); // a-b
    CHECK_FALSE(has_op_mod<T, linalg::quat<T>>::value); // a%b
    CHECK_FALSE(has_op_lsh<T, linalg::quat<T>>::value); // a<<b
    CHECK_FALSE(has_op_rsh<T, linalg::quat<T>>::value); // a>>b
    CHECK_FALSE(has_op_int<T, linalg::quat<T>>::value); // a&b
    CHECK_FALSE(has_op_xor<T, linalg::quat<T>>::value); // a^b
    CHECK_FALSE(has_op_un <T, linalg::quat<T>>::value); // a|b
    
    // Only operator +=, -=, *= are defined for quat $= quat
    CHECK(has_op_assign_add<linalg::quat<T> &, linalg::quat<T>>::value); // a+=b
    CHECK(has_op_assign_sub<linalg::quat<T> &, linalg::quat<T>>::value); // a-=b
    CHECK(has_op_assign_mul<linalg::quat<T> &, linalg::quat<T>>::value); // a*=b
    CHECK_FALSE(has_op_assign_div<linalg::quat<T> &, linalg::quat<T>>::value); // a/=b
    CHECK_FALSE(has_op_assign_mod<linalg::quat<T> &, linalg::quat<T>>::value); // a%=b
    CHECK_FALSE(has_op_assign_lsh<linalg::quat<T> &, linalg::quat<T>>::value); // a<<=b
    CHECK_FALSE(has_op_assign_rsh<linalg::quat<T> &, linalg::quat<T>>::value); // a>>=b
    CHECK_FALSE(has_op_assign_int<linalg::quat<T> &, linalg::quat<T>>::value); // a&=b
    CHECK_FALSE(has_op_assign_xor<linalg::quat<T> &, linalg::quat<T>>::value); // a^=b
    CHECK_FALSE(has_op_assign_un <linalg::quat<T> &, linalg::quat<T>>::value); // a|=b

    // Only operator *= and /= are defined for quat $= scalar
    CHECK(has_op_assign_mul<linalg::quat<T> &, T>::value); // a*=b 
    CHECK(has_op_assign_div<linalg::quat<T> &, T>::value); // a/=b 
    CHECK_FALSE(has_op_assign_mod<linalg::quat<T> &, T>::value); // a%=b 
    CHECK_FALSE(has_op_assign_add<linalg::quat<T> &, T>::value); // a+=b 
    CHECK_FALSE(has_op_assign_sub<linalg::quat<T> &, T>::value); // a-=b 
    CHECK_FALSE(has_op_assign_lsh<linalg::quat<T> &, T>::value); // a<<=b 
    CHECK_FALSE(has_op_assign_rsh<linalg::quat<T> &, T>::value); // a>>=b 
    CHECK_FALSE(has_op_assign_int<linalg::quat<T> &, T>::value); // a&=b 
    CHECK_FALSE(has_op_assign_xor<linalg::quat<T> &, T>::value); // a^=b 
    CHECK_FALSE(has_op_assign_un <linalg::quat<T> &, T>::value); // a|=b
}

// TODO: Check mixed-element expressions