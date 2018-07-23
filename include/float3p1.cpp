#pragma once

// This is a class that has 4 floats, but omits the last element from
// scalar multiplies.

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
    // There is no null constructor; user must initialize
    FLOAT3p1(FLOAT _x) { x=y=z=w=_x; }   // Use with 0.0, normally
    FLOAT3p1(FLOAT _x, FLOAT _y, FLOAT _z, FLOAT _w) { 
	x=_x; y=_y; z=_z; ; w=_w; }
    FLOAT3p1(FLOAT3 _v, FLOAT _w) { 
	x = _v.x; y = _v.y; z = _v.z; w=_w; }
    // This constructor from a float3 will force w=0
    FLOAT3p1(FLOAT3 _v) { 
	x = _v.x; y = _v.y; z = _v.z; w=0.0; }
    ~FLOAT3p1();

    // Add two float3p1's, inluding w.
    // We don't supply subtraction, because the behavior on w is ambiguous
    inline FLOAT3p1& operator += ( const FLOAT3p1 rhs ) {
        x += rhs.x; y += rhs.y; z += rhs.z; w += rhs.w; return *this; }

    // Add or subtract a float3, omitting w
    inline FLOAT3p1& operator += ( const FLOAT3 rhs ) {
        x += rhs.x; y += rhs.y; z += rhs.z; return *this; }
    inline FLOAT3p1& operator -= ( const FLOAT3 rhs ) {
        x -= rhs.x; y -= rhs.y; z -= rhs.z; return *this; }

    // Multiply and divide by scalars, omitting w
    inline FLOAT3p1& operator *= ( const float rhs ) {
        x *= rhs; y *= rhs; z *= rhs; return *this; }
    inline FLOAT3p1& operator *= ( const double rhs ) {
        x *= rhs; y *= rhs; z *= rhs; return *this; }
    inline FLOAT3p1& operator /= ( const float rhs ) {
        x /= rhs; y /= rhs; z /= rhs; return *this; }
    inline FLOAT3p1& operator /= ( const double rhs ) {
        x /= rhs; y /= rhs; z /= rhs; return *this; }

    // Supply a type cast to float3
    inline operator FLOAT3() { return FLOAT3(x,y,z); }

    // Provide a norm2() to act on the float3 portion
    FLOAT norm() { return sqrt(x*x+y*y+z*z); }
    FLOAT norm2() { return x*x+y*y+z*z; }

    // We don't supply a type cast to float4, but that could be easily
    // done with pointers, since the space is identical.
};

