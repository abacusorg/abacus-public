// Copyright 2012-2025 The Abacus Developers
// SPDX-License-Identifier: GPL-3.0-or-later

/*
ThreeVector class

Time-stamp: <threevector.cc on Saturday, 8 September, 2012 at 16:05:45 MST (philip)>

 */

#ifndef __THREEVECTOR_CC__
#define __THREEVECTOR_CC__
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include <cassert>

#include "promote_numeric.h"

enum IntervalType {
    Closed,
    HalfOpen
};

template <class T>
class ThreeVector {
public:
    T x, y, z;
    
    // constructors

    // Default is no initialization, just like built-in types
    inline ThreeVector() = default;

    template <class U>
    inline ThreeVector( const ThreeVector<U>& other )
    :  x(static_cast<T>(other.x)),
        y(static_cast<T>(other.y)),
        z(static_cast<T>(other.z))
    {}
    
    // want to be able to write: ThreeVector<double> x(2,2.5,4)
    template <class U, class V, class W>
    inline ThreeVector( const U rhsx, const V rhsy, W const rhsz ) 
    : x(static_cast<T>(rhsx)), 
        y(static_cast<T>(rhsy)), 
        z(static_cast<T>(rhsz))
    {}
    
    // broadcast a scalar
    // use SFINAE to avoid ambiguity with custom cast operators
    template <class U, typename = std::enable_if_t<std::is_convertible_v<U, T>>>
    inline explicit ThreeVector( const U& s)
    : x(static_cast<T>(s)),
        y(static_cast<T>(s)),
        z(static_cast<T>(s))
    {}

    // allow access as a three-vector 
    inline T operator [] ( const size_t i ) const {
        assert( i < 3 );
        return *(&x+i);
    }
    inline T& operator [] ( const size_t i ) {
        assert( i < 3 );
        return *(&x+i);
    }
    
    // assignment
    template <class U>
    inline ThreeVector<T>& operator = ( const ThreeVector<U>& other ) {
        x = static_cast<T>(other.x);
        y = static_cast<T>(other.y);
        z = static_cast<T>(other.z);
        return *this;
    }
    
    //unary operations (sign)
    inline const ThreeVector<T>& operator +() {
        return *this;
    }

    inline ThreeVector<T> operator -() {
        return ThreeVector<T>(-x, -y, -z);
    }

    template <class U>
    inline ThreeVector<T>& operator += ( const ThreeVector<U>& rhs ) {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        
        return *this;
    }
    
    template <class U>
    inline ThreeVector<T>& operator -= ( const ThreeVector<U>& rhs ) {
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
        
        return *this;
    }
    
    template <class U>
    inline ThreeVector<T>& operator *= ( const ThreeVector<U>& rhs ) {
        x *= rhs.x;
        y *= rhs.y;
        z *= rhs.z;
        
        return *this;
    }
    
    template <class U>
    inline ThreeVector<T>& operator *= ( const U rhs ) {
        x *= rhs;
        y *= rhs;
        z *= rhs;
        return *this;
    }

    template <class U>
    inline ThreeVector<T>& operator /= ( const U rhs ) {
        x /= rhs;
        y /= rhs;
        z /= rhs;
        return *this;
    }

    inline typename PromoteNumeric<T, float>::type norm() const {
        return sqrt( x * x + y * y + z * z );
    }

    inline typename PromoteNumeric<T, float>::type norm2() const {
        return ( x * x + y * y + z * z );
    }

    inline T maxcomponent() const {
        T maxxy = (x>y)?x:y;
        return maxxy>z?maxxy:z;
    }
    inline T mincomponent() const {
        T maxxy = (x<y)?x:y;
        return maxxy<z?maxxy:z;
    }
    inline T maxabscomponent() const {
        T maxxy = (fabs(x)>fabs(y))?fabs(x):fabs(y);
        return maxxy>fabs(z)?maxxy:fabs(z);
    }

    inline ThreeVector<T> round() const {
        return ThreeVector<T>(::round(x), ::round(y), ::round(z));
    }

    template<class U>
    inline typename PromoteNumeric<T, U>::type dot(const ThreeVector<U>& rhs) const {
        return x * rhs.x + y * rhs.y + z * rhs.z;
    }

    template <class U>
    inline ThreeVector<typename PromoteNumeric<T, U>::type>
    cross( const ThreeVector<U>& rhs ) const {
        return ThreeVector<typename PromoteNumeric<T, U>::type>(
                                                          y * rhs.z - z * rhs.y,
                                                          z * rhs.x - x * rhs.z,
                                                          x * rhs.y - y * rhs.x);
    }

    // abs. diff. between *this and rhs is less than epsilon in each component
    template <class U, class V>
    inline int absclose( const ThreeVector<U>& rhs, const V epsilon ) const {
        ThreeVector<typename PromoteNumeric<T, U>::type> diff;
        diff = *this - rhs;
        return 
            fabs(diff.x)  < epsilon  && 
            fabs(diff.y)  < epsilon  && 
            fabs(diff.z)  < epsilon;
    }

    // rel. diff. between *this and rhs is less than epsilon in each component
    template <class U, class V>
    inline int relclose( const ThreeVector<U>& rhs, const V epsilon ) const {
        ThreeVector<typename PromoteNumeric<T, U>::type> sum, diff;
        sum.x = fabs(x) + fabs(rhs.x);
        sum.y = fabs(y) + fabs(rhs.y);
        sum.z = fabs(z) + fabs(rhs.z);
        diff = *this - rhs;
        return 
            ( 2*fabs(diff.x) / sum.x ) < epsilon  && 
            ( 2*fabs(diff.y) / sum.y ) < epsilon  && 
            ( 2*fabs(diff.z) / sum.z ) < epsilon;
    }

    inline int is_finite() const {
        return isfinite(x) && isfinite(y) && isfinite(z);
    }

    // relational operators
    template <class U>
    inline bool operator == ( const ThreeVector<U>& rhs ) const {
        return ( x == rhs.x && y == rhs.y && z == rhs.z );
    }
    
    template <class U>
    inline bool operator != ( const ThreeVector<U>& rhs ) const {
        return ( x != rhs.x || y != rhs.y || z != rhs.z );
    }

    template <class U>
    inline bool operator < ( const ThreeVector<U>& rhs ) const {
        if( x < rhs.x && y < rhs.y && z < rhs.z ) return true;
        return false;
    }

    template <class U>
    inline bool operator <= ( const ThreeVector<U>& rhs ) const {
        if( x <= rhs.x && y <= rhs.y && z <= rhs.z ) return true;
        return false;
    }
    
    template <class U>
    inline bool operator > ( const ThreeVector<U>& rhs ) const {
        if( x > rhs.x && y > rhs.y && z > rhs.z ) return true;
        return false;
    }

    template <class U>
    inline bool operator >= ( const ThreeVector<U>& rhs ) const {
        if( x >= rhs.x && y >= rhs.y && z >= rhs.z ) return true;
        return false;
    }

    // stream operator: keep the format settings from being destroyed by the
    // non-numeric characters output
    inline friend std::ostream& operator <<( std::ostream& o, const ThreeVector<T>& v ) {
        std::streamsize tmpw = o.width();
        std::streamsize tmpp = o.precision();
        char tmps = o.fill();
        std::ios::fmtflags tmpf = o.flags();  // format flags like "scientific" and "left" and "showpoint"
        o << std::setw(1);
        o << "(";  
        o.flags(tmpf); o << std::setfill(tmps) << std::setprecision(tmpp) << std::setw(tmpw);
        o << v.x;
        o << ",";
        o.flags(tmpf); o << std::setfill(tmps) << std::setprecision(tmpp) << std::setw(tmpw);
        o << v.y;
        o << ",";
        o.flags(tmpf); o << std::setfill(tmps) << std::setprecision(tmpp) << std::setw(tmpw);
        o << v.z;
        o << ")";
        return o;
    }
    
    inline ThreeVector<T> zero() {
        x = 0;
        y = 0;
        z = 0;
        return *this;
    }

    // return true if all components are on [a,b] or [a,b)
    template <IntervalType interval = IntervalType::Closed>
    inline bool inrange(T low, T hi) {
        if constexpr (interval == HalfOpen) {
            return (x >= low && x < hi) && (y >= low && y < hi) && (z >= low && z < hi);
        } else {
            return (x >= low && x <= hi) && (y >= low && y <= hi) && (z >= low && z <= hi);
        }
    }

    template <IntervalType interval = IntervalType::Closed>
    inline bool inrange(ThreeVector<T> low, ThreeVector<T> hi) {
        if constexpr (interval == HalfOpen) {
            return (x >= low.x && x < hi.x) && (y >= low.y && y < hi.y) && (z >= low.z && z < hi.z);
        } else {
            return (x >= low.x && x <= hi.x) && (y >= low.y && y <= hi.y) && (z >= low.z && z <= hi.z);
        }
    }
};

// componentwise addition and subtraction
template <class T, class U>
inline ThreeVector<typename PromoteNumeric<T, U>::type>
operator + (const ThreeVector<T>& lhs, const ThreeVector<U>& rhs ) {
    return ThreeVector<typename PromoteNumeric<T, U>::type> (
                                                       lhs.x + rhs.x,
                                                       lhs.y + rhs.y,
                                                       lhs.z + rhs.z
                                                       );
}

template <class T, class U>
inline ThreeVector<typename PromoteNumeric<T, U>::type>
operator - (const ThreeVector<T>& lhs, const ThreeVector<U>& rhs ) {
    return ThreeVector<typename PromoteNumeric<T, U>::type> (
                                                       lhs.x - rhs.x,
                                                       lhs.y - rhs.y,
                                                       lhs.z - rhs.z
                                                       );
}

// left and right multiplication by a scalar
template <class T, class U>
inline ThreeVector<typename PromoteNumeric<T, U>::type>
operator * ( const ThreeVector<T>& lhs, const U rhs ) {
    return ThreeVector<typename PromoteNumeric<T, U>::type> (
                                                       lhs.x * rhs,
                                                       lhs.y * rhs,
                                                       lhs.z * rhs
                                                       );
}
template <class T, class U>
inline ThreeVector<typename PromoteNumeric<T, U>::type>
operator * ( const T lhs, const ThreeVector<U>& rhs ) {
    return ThreeVector<typename PromoteNumeric<T, U>::type> (
                                                       lhs * rhs.x,
                                                       lhs * rhs.y,
                                                       lhs * rhs.z
                                                       );
}

// element-wise multiplication
template <class T, class U>
inline ThreeVector<typename PromoteNumeric<T, U>::type>
operator * ( const ThreeVector<T>& lhs, const ThreeVector<U>& rhs ) {
    return ThreeVector<typename PromoteNumeric<T, U>::type> (
                                                       lhs.x * rhs.x,
                                                       lhs.y * rhs.y,
                                                       lhs.z * rhs.z
                                                       );
}

// right division by a scalar
template <class T, class U>
inline ThreeVector<typename PromoteNumeric<T, U>::type>
operator / ( const ThreeVector<T>& lhs, const U rhs ) {
    return ThreeVector<typename PromoteNumeric<T, U>::type> (
                                                       lhs.x / rhs,
                                                       lhs.y / rhs,
                                                       lhs.z / rhs
                                                       );
}

// three most common cases
using double3 = ThreeVector<double>;
using float3 = ThreeVector<float>;
using integer3 = ThreeVector<int>;

#ifdef _OPENMP
#pragma omp declare reduction(+:float3:omp_out += omp_in) initializer(omp_priv = float3(0.0f))
#pragma omp declare reduction(+:double3:omp_out += omp_in) initializer(omp_priv = double3(0.0))
#pragma omp declare reduction(+:integer3:omp_out += omp_in) initializer(omp_priv = integer3(0))
#endif

#endif // __THREEVECTOR_CC__
