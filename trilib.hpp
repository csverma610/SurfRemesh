#pragma once

#include "veclib.hpp"

#define ANGLE_IN_DEGREES  0
#define ANGLE_IN_RADIANS  1

using namespace JMath;

template<class T, size_t N>
inline T minlength( const std::array<T,N> &pa,
                    const std::array<T,N> &pb,
                    const std::array<T,N> &pc)
{
    T a  =  length( pb, pc );
    T b  =  length( pc, pa );
    T c  =  length( pa, pb );

    return min_value(a,b,c);
}

template<class T, size_t N>
inline T maxlength( const std::array<T,N> &pa,
                    const std::array<T,N> &pb,
                    const std::array<T,N> &pc)
{
    T a  =  length( pb, pc );
    T b  =  length( pc, pa );
    T c  =  length( pa, pb );

    return max_value(a,b,c);
}

template<class T, size_t N>
inline std::array<T,3> angles( const std::array<T,N> &pa,
                               const std::array<T,N> &pb,
                               const std::array<T,N> &pc, 
			       int measure = ANGLE_IN_DEGREES)
{
    std::array<T,3> angles = {0, 0, 0};

    double a2   =  length2( pb, pc );
    double b2   =  length2( pc, pa );
    double c2   =  length2( pa, pb );
    double cosA =  (b2 + c2 - a2)/(2*sqrt(b2*c2) );
    double cosB =  (a2 + c2 - b2)/(2*sqrt(a2*c2) );
    double cosC =  (a2 + b2 - c2)/(2*sqrt(a2*b2) );

    if( cosA >  1.0) cosA =  1.0;
    if( cosA < -1.0) cosA = -1.0;
    angles[0] = acos(cosA);

    if( cosB >  1.0) cosB =  1.0;
    if( cosB < -1.0) cosB = -1.0;
    angles[1] = acos(cosB);

    if( cosC >  1.0) cosC =  1.0;
    if( cosC < -1.0) cosC = -1.0;
    angles[2] = acos(cosC);

    if( measure ==  ANGLE_IN_DEGREES) {
        angles[0] *= 180/M_PI;
        angles[1] *= 180/M_PI;
        angles[2] *= 180/M_PI;
    }
    return angles;
}

///////////////////////////////////////////////////////////////////////////////

template<class T, size_t N>
inline T angleAt( const std::array<T,N> &pa,
                  const std::array<T,N> &pb,
                  const std::array<T,N> &pc, 
		  int measure = ANGLE_IN_DEGREES)
{
    T a2   =  length2( pb, pc );
    T b2   =  length2( pc, pa );
    T c2   =  length2( pa, pb );
    double cosA =  (b2 + c2 - a2)/(2*sqrt(b2*c2) );

    if( cosA >  1.0) cosA =  1.0;
    if( cosA < -1.0) cosA = -1.0;

    double angle;
    if( measure == ANGLE_IN_DEGREES)
        angle = 180*acos(cosA)/M_PI;
    else
        angle = acos(cosA);
    return angle;
}

////////////////////////////////////////////////////////////////////////////////

template<class T, size_t N>
std::pair<T,int> maxangle( const std::array<T,N> &pa,
                           const std::array<T,N> &pb,
                           const std::array<T,N> &pc, 
			   int measure = ANGLE_IN_DEGREES)
{
    T a2   =  length2( pb, pc );
    T b2   =  length2( pc, pa );
    T c2   =  length2( pa, pb );

    T maxlen = max_value(a2,b2,c2);

    std::pair<T,int> result;

    if( maxlen == a2) {
        double cosA =  (b2 + c2 - a2)/(2*sqrt(b2*c2) );
        if( cosA >  1.0) cosA =  1.0;
        if( cosA < -1.0) cosA = -1.0;
        double angle = acos(cosA);
        if( measure == ANGLE_IN_DEGREES)  angle *= 180.0/M_PI;
        result.first  = angle;
        result.second = 0;
    }

    if( maxlen == b2 ) {
        double cosB =  (a2 + c2 - b2)/(2*sqrt(a2*c2) );
        if( cosB >  1.0) cosB =  1.0;
        if( cosB < -1.0) cosB = -1.0;
        double angle  = acos(cosB);
        if( measure == ANGLE_IN_DEGREES)  angle *= 180.0/M_PI;
        result.first  = angle;
        result.second = 1;
    }

    if( maxlen == c2 ) {
        double cosC =  (a2 + b2 - c2)/(2*sqrt(a2*b2) );
        if( cosC >  1.0) cosC =  1.0;
        if( cosC < -1.0) cosC = -1.0;
        double angle = acos(cosC);
        if( measure == ANGLE_IN_DEGREES)  angle *= 180.0/M_PI;
        result.first  = angle;
        result.second = 2;
    }

    return result;
}

////////////////////////////////////////////////////////////////////////////////

template<class T, size_t N>
std::pair<T,int> minangle( const std::array<T,N> &pa,
                           const std::array<T,N> &pb,
                           const std::array<T,N> &pc, 
			   int measure = ANGLE_IN_DEGREES)
{
    double a2   =  length2( pb, pc );
    double b2   =  length2( pc, pa );
    double c2   =  length2( pa, pb );

    double minlen = min_value(a2,b2,c2);

    std::pair<T,int> result;

    if( minlen == a2) {
        double cosA =  (b2 + c2 - a2)/(2*sqrt(b2*c2) );
        if( cosA >  1.0) cosA =  1.0;
        if( cosA < -1.0) cosA = -1.0;
        double angle = acos(cosA);
        if( measure == ANGLE_IN_DEGREES)  angle *= 180.0/M_PI;
        result.first  = angle;
        result.second = 0;
    }

    if( minlen == b2 ) {
        double cosB =  (a2 + c2 - b2)/(2*sqrt(a2*c2) );
        if( cosB >  1.0) cosB =  1.0;
        if( cosB < -1.0) cosB = -1.0;
        double angle  = acos(cosB);
        if( measure == ANGLE_IN_DEGREES)  angle *= 180.0/M_PI;
        result.first  = angle;
        result.second = 1;
    }

    if( minlen == c2 ) {
        double cosC =  (a2 + b2 - c2)/(2*sqrt(a2*b2) );
        if( cosC >  1.0) cosC =  1.0;
        if( cosC < -1.0) cosC = -1.0;
        double angle = acos(cosC);
        if( measure == ANGLE_IN_DEGREES)  angle *= 180.0/M_PI;
        result.first  = angle;
        result.second = 2;
    }

    return result;
}
//
////////////////////////////////////////////////////////////////////////////////
//
template<class T, size_t N>
inline bool isObtuse( const std::array<T,N> &pa,
                      const std::array<T,N> &pb,
                      const std::array<T,N> &pc)
{
    auto result = maxangle(pa,pb,pc);
    if( result.first > 90.0) return 1;
    return 0;
}
////////////////////////////////////////////////////////////////////////////////

template<class T, size_t N>
inline bool isDegenerate( const std::array<T,N> &pa,
                          const std::array<T,N> &pb,
                          const std::array<T,N> &pc)
{
    auto result = maxangle(pa,pb,pc);
    if( result.first > 179.999) return 1;
    return 0;
}

////////////////////////////////////////////////////////////////////////////////

template<class T, size_t N>
inline bool isAcute( const std::array<T,N> &pa,
                     const std::array<T,N> &pb,
                     const std::array<T,N> &pc)
{
    auto result = maxangle(pa,pb,pc);
    if( result.first <= 90.0) return 1;
    return 0;
}
////////////////////////////////////////////////////////////////////////////////

template<class T>
inline std::array<T,3> normal( const std::array<T,3> &p0,
                               const std::array<T,3> &p1,
                               const std::array<T,3> &p2)
{
    auto p1p0   = make_vector( p1, p0);
    auto p2p0   = make_vector( p2, p0);
    auto normal = cross_product( p1p0, p2p0);

    double mag = magnitude( normal );
    normal[0] /= mag;
    normal[1] /= mag;
    normal[2] /= mag;

    return normal;
}

////////////////////////////////////////////////////////////////////////////////
template<class T, size_t N>
inline T area( const std::array<T,N> &pa,
               const std::array<T,N> &pb,
               const std::array<T,N> &pc)
{
    T a     =  length( pb, pc );
    T b     =  length( pc, pa );
    T c     =  length( pa, pb );
    T s     =  0.5*(a+b+c);
    T heron =  sqrt(s*(s-a)*(s-b)*(s-c));
    return heron;
}
////////////////////////////////////////////////////////////////////////////////
template<class T, size_t N>
inline std::array<T,N> centroid( const std::array<T,N> &pa,
                                 const std::array<T,N> &pb,
                                 const std::array<T,N> &pc)
{
    std::array<T,N> c;
    for( int i = 0; i < N; i++) c[i] = (pa[i] + pb[i] + pc[2])/3.0;
    return c;
}
////////////////////////////////////////////////////////////////////////////////
template<class T, size_t N>
inline std::array<T,3> barycoordinates( const std::array<T,N> &pa,
                                        const std::array<T,N> &pb,
                                        const std::array<T,N> &pc, 
					const std::array<T,N> &queryPoint)
{
    std::array<T,3> bcoords;

    T total_area = area(pa,pb,pc);

    bcoords[0] = area(pb,pc,queryPoint)/total_area;
    bcoords[1] = area(pc,pa,queryPoint)/total_area;
    bcoords[2] = area(pa,pb,queryPoint)/total_area;

    return bcoords;
}
////////////////////////////////////////////////////////////////////////////////
template<class T, size_t N>
inline std::array<T,N> circumcenter( const std::array<T,N> &pa,
                                     const std::array<T,N> &pb,
                                     const std::array<T,N> &pc)
{
   // Source : Wikipedia ...
    std::array<T,N> coords;
    T a   =  length( pb, pc );
    T b   =  length( pc, pa );
    T c   =  length( pa, pb );

    T u   =  a*a*(b*b + c*c - a*a);
    T v   =  b*b*(c*c + a*a - b*b);
    T w   =  c*c*(a*a + b*b - c*c);

    for( int i = 0; i < N; i++) 
         coords[i] = (u*pa[i] + v*pb[i] + w*pc[i] )/(u+v+w);

    return coords;
}

////////////////////////////////////////////////////////////////////////////////

template<class T, size_t N>
inline T circumradius( const std::array<T,N> &pa,
                       const std::array<T,N> &pb,
                       const std::array<T,N> &pc)
{
    T a  =  length( pb, pc );
    T b  =  length( pc, pa );
    T c  =  length( pa, pb );
    T s  =  0.5*(a+b+c);

    T r  = 0.25*a*b*c/sqrt(s*(s-a)*(s-b)*(s-c));
    return r;
}

////////////////////////////////////////////////////////////////////////////////

template<class T, size_t N>
inline std::array<T,N> incenter( const std::array<T,N> &pa,
                                 const std::array<T,N> &pb,
                                 const std::array<T,N> &pc)
{
    std::array<T,N> coords;
    T a  =  length( pb, pc );
    T b  =  length( pc, pa );
    T c  =  length( pa, pb );
    T t  =  a + b + c;

    for( int i = 0; i < N; i++) 
        coords[0] = (a*pa[i] + b*pb[i] + c*pc[i] )/t;

    return coords;
}
////////////////////////////////////////////////////////////////////////////////

template<class T, size_t N>
inline T inradius( const std::array<T,N> &pa,
                   const std::array<T,N> &pb,
                   const std::array<T,N> &pc)
{
    T a  =  length( pb, pc );
    T b  =  length( pc, pa );
    T c  =  length( pa, pb );
    T s  =  a+b+c;

    T r     = 0.5*sqrt((b+c-a)*(c+a-b)*(a+b-c)/s);
    return r;
}
////////////////////////////////////////////////////////////////////////////////
