#pragma once

#include <stdlib.h>
#include <math.h>
#include <values.h>
#include <vector>
#include <algorithm>
#include <array>

using Point2I = std::array<int,2>;
using Point3I = std::array<int,3>;
using Point2F = std::array<float,2>;
using Point3F = std::array<float,3>;
using Point4F = std::array<float,4>;

using Point2D = std::array<double,2>;
using Point3D = std::array<double,3>;
using Point4D = std::array<double,4>;

using Array3F = std::array<float,3>;
using Array4F = std::array<float,4>;

using Array2D = std::array<double,2>;
using Array3D = std::array<double,3>;
using Array4D = std::array<double,4>;

using Array2I = std::array<int,2> ;
using Array3I = std::array<int,3> ;
using Array4I = std::array<int,4> ;

using Array2F = std::array<float,2>;
using Array3F = std::array<float,3>;

using Vec2D = std::array<double,2>;
using Vec3D = std::array<double,3>;
using Vec4D = std::array<double,4>;

namespace JMath
{
template<class T>
inline T max_value( const T &a, const T &b, const T &c)
{
    return std::max(a,std::max(b, c));
}

template<class T>
inline T min_value( const T &a, const T &b, const T &c)
{
    return std::min(a,std::min(b, c));
}

template<class T>
inline T min_value( const T &a, const T &b, const T &c, const T &d)
{
    return std::min(d, std::min(a,std::min(b, c)));
}

template<class T>
T random_value(T minVal, T maxVal)
{
    return minVal + drand48()*(maxVal - minVal);
}

template<class T, size_t N>
inline T length(const std::array<T,N> &A, const std::array<T,N> &B)
{
    T  sum = 0;
    for( int i = 0; i < N; i++) {
	    T dl = A[i]-B[i];
	    sum  +=  dl*dl;
    }
    return sqrt(sum);
}

template<class T, size_t N>
inline T length2( const std::array<T,N> &A, const std::array<T,N> &B)
{
    T  sum = 0;
    for( int i = 0; i < N; i++) {
	    T dl  = A[i]-B[i];
	    sum  +=  dl*dl;
    }
    return sum;
}

template<class T, size_t N>
inline T magnitude( const std::array<T,N> &A )
{
    T  sum = 0;
    for( int i = 0; i < N; i++) sum += A[i]*A[i];
    return sqrt(sum);
}

template<class T, size_t N>
inline T dot_product( const std::array<T,N> &A, const std::array<T,N> &B)
{
    T  sum = 0;
    for( int i = 0; i < N; i++) sum += A[i]*B[i];
    return sqrt(sum);
}

template<class T>
inline std::array<T,3> cross_product( const std::array<T,3> &A, const std::array<T,3> &B)
{
    std::array<T,3> C;
    C[0] = A[1]*B[2] - A[2]*B[1];
    C[1] = A[2]*B[0] - A[0]*B[2];
    C[2] = A[0]*B[1] - A[1]*B[0];
    return C;
}

template<class T, size_t N>
inline std::array<T,3> make_vector( const std::array<T,N> &head, const std::array<T,N> &tail)
{
    std::array<T,N> xyz;

    for( int i = 0; i < N; i++) xyz[i] = head[i] - tail[i];
    return xyz;
}

///////////////////////////////////////////////////////////////////////////////
template<class T, size_t N>
inline std::array<T,N> unit_vector( const std::array<T,N> &vec)
{
    double dl  = magnitude(vec);
    std::array<T,N>  uvec;

    for( int i = 0; i < N; i++) uvec[i] = vec[i]/dl;

    return uvec;
}
///////////////////////////////////////////////////////////////////////////////

template<class T>
inline int unit_vector( const Point3D &head, const Point3D &tail)
{
    std::array<T,3> uvec;
    uvec = make_vector( head, tail);
    double dl  = magnitude(uvec);

    uvec[0]  /=  dl;
    uvec[1]  /=  dl;
    uvec[2]  /=  dl;
    return uvec;
}

template<class T>
inline T mean_value( const std::vector<T> &v)
{
    assert( !v.empty() );
    std::vector<T> tmp(v);
    std::sort( tmp.begin(), tmp.end() );
    return tmp[v.size()/2];
}

template<class T>
inline T average_value( const std::vector<T> &v)
{
    size_t nsize = v.size();

    T sum = 0.0;
    for( size_t i = 0; i < nsize; i++)
        sum += v[i];
    T avg = sum/(double)nsize;
    return avg;
}

template<class T>
inline T standard_deviation( const std::vector<T> &v)
{
    assert( !v.empty() );

    size_t nsize = v.size();

    T avg = average_value( v );

    T sum = 0.0;
    for( size_t i = 0; i < nsize; i++) {
        double vm = v[i] - avg;
        sum += vm*vm;
    }
    return sqrt(1.0/(double)(nsize-1)*sum );
}

///////////////////////////////////////////////////////////////////////////////
template<class T>
inline T angle( T x1, T y1, T x2, T y2)
{
    double theta1 = atan2( (double)y1, (double)x1);
    double theta2 = atan2( (double)y2, (double)x2);

    double dtheta = theta2-theta1;

    if( dtheta >  M_PI) dtheta -= 2.0*M_PI;
    if( dtheta < -M_PI) dtheta += 2.0*M_PI;

    return dtheta;
}
///////////////////////////////////////////////////////////////////////////////

template<class T, size_t N>
inline T angle( const std::array<T,N> &A, const std::array<T,N> &B)
{
    double AB = dot_product(A,B);
    double Am = magnitude(A);
    double Bm = magnitude(B);

    if( Am < 1.0E-15 || Bm < 1.0E-15) return 0.0;

    double x = AB/(Am*Bm);

    if( x > 1.0)  x = 1.0;
    if( x < -1.0) x = -1.0;

    return acos(x);
}

}
