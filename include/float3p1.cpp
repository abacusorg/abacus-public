// Copyright 2012-2025 The Abacus Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

// This is a class that has 4 floats, but omits the last element from
// scalar multiplies.  We use this for cases like tracking acceleration
// in the first three elements, and FoF link counts in the fourth.

// Because we don't have float3 and double3 templated, it's not 
// so easy to template this.

// Note that when adding a float3p1 to a float3, the ordering might
// matter.  float3p1+float3 will produce a float3p1, carrying the w along. 
// float3+float3p1 might be undefined, or might trigger the typecast
// and hence strip the w.  Not sure of the rules.

class FLOAT3p1 {
  public:
    FLOAT x, y, z, w;

    // Constructors and destructors
    // The null constructor doesn't initialize
    FLOAT3p1() = default; 
    FLOAT3p1(FLOAT _x) { x=y=z=w=_x; }   // Use with 0.0, normally
    FLOAT3p1(FLOAT _x, FLOAT _y, FLOAT _z, FLOAT _w) { 
        x=_x; y=_y; z=_z; ; w=_w; }
    FLOAT3p1(FLOAT3 _v, FLOAT _w) { 
        x = _v.x; y = _v.y; z = _v.z; w=_w; }
    // This constructor from a float3 will force w=0
    FLOAT3p1(FLOAT3 _v) { 
        x = _v.x; y = _v.y; z = _v.z; w=0.0; }
    // ~FLOAT3p1() { }

    // Add two float3p1's, inluding w.
    // We don't supply subtraction, because the behavior on w is ambiguous
    inline FLOAT3p1 operator + ( const FLOAT3p1 &rhs ) const {
        return FLOAT3p1(x+rhs.x, y+rhs.y, z+rhs.z, w+rhs.w); }
    inline FLOAT3p1& operator += ( const FLOAT3p1 &rhs ) {
        x+=rhs.x; y+=rhs.y; z+=rhs.z; w+=rhs.w; return *this; }

    // Add or subtract a float3, omitting w
    inline FLOAT3p1 operator + ( const FLOAT3 &rhs ) const {
        return FLOAT3p1(x+rhs.x, y+rhs.y, z+rhs.z, w); }
    inline FLOAT3p1 operator - ( const FLOAT3 &rhs ) const {
        return FLOAT3p1(x-rhs.x, y-rhs.y, z-rhs.z, w); }

    // Multiply and divide by scalars, omitting w
    inline FLOAT3p1 operator * ( const float &rhs ) const {
        return FLOAT3p1(x*rhs, y*rhs, z*rhs, w); }
    inline FLOAT3p1 operator *= ( const float &rhs ) {
    x*=rhs; y*=rhs; z*=rhs; return *this; }
    inline FLOAT3p1 operator * ( const double &rhs ) const {
        return FLOAT3p1(x*rhs, y*rhs, z*rhs, w); }
    inline FLOAT3p1 operator / ( const float &rhs ) const {
        float inv = 1/rhs; return FLOAT3p1(x*inv, y*inv, z*inv, w); }
    inline FLOAT3p1 operator / ( const double &rhs ) const {
        double inv = 1/rhs; return FLOAT3p1(x*inv, y*inv, z*inv, w); }

    // Supply a type cast to float3
    explicit operator FLOAT3() const { return FLOAT3(x,y,z); }

    // Provide a norm2() to act on the float3 portion
    FLOAT norm() const { return sqrt(x*x+y*y+z*z); }
    FLOAT norm2() const { return x*x+y*y+z*z; }

    // We don't supply a type cast to float4, but that could be easily
    // done with pointers, since the space is identical.
};
