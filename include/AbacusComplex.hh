/*

AbacusComplex.hh

A complex number class that is trivially constructible.

std::complex is technically not trivially constructible, so it can't be used
(safely) with memset. It also can lead to unexpected initializations when
making arrays of std::complex. This is especially bad when the initialization
is done in a single-threaded region and we want it done in parallel.

It's possible that the C++20 std::complex will be trivially constructible, but
until then, we can use this class.

*/

#ifndef __ABACUS_COMPLEX_HH__
#define __ABACUS_COMPLEX_HH__

#include <complex>
#include <fftw3.h>

template <typename T>
class AbacusComplex {
public:
    T real;
    T imag;

    // Default constructor
    AbacusComplex() = default;

    // Parameterized constructor
    AbacusComplex(T r, T i) : real(r), imag(i) {}

    // Copy constructor
    AbacusComplex(const AbacusComplex&) = default;

    // Move constructor
    AbacusComplex(AbacusComplex&&) = default;

    // Copy assignment operator
    AbacusComplex& operator=(const AbacusComplex&) = default;

    // Move assignment operator
    AbacusComplex& operator=(AbacusComplex&&) = default;

    // Destructor
    ~AbacusComplex() = default;

    // Addition operator
    AbacusComplex operator+(const AbacusComplex& other) const {
        return AbacusComplex(real + other.real, imag + other.imag);
    }

    // Addition operator with scalar
    AbacusComplex operator+(T scalar) const {
        return AbacusComplex(real + scalar, imag);
    }

    // In-place addition operator
    AbacusComplex& operator+=(const AbacusComplex& other) {
        real += other.real;
        imag += other.imag;
        return *this;
    }

    // Subtraction operator
    AbacusComplex operator-(const AbacusComplex& other) const {
        return AbacusComplex(real - other.real, imag - other.imag);
    }

    // Subtraction operator with scalar
    AbacusComplex operator-(T scalar) const {
        return AbacusComplex(real - scalar, imag);
    }

    // In-place subtraction operator
    AbacusComplex& operator-=(const AbacusComplex& other) {
        real -= other.real;
        imag -= other.imag;
        return *this;
    }

    // Multiplication operator
    AbacusComplex operator*(const AbacusComplex& other) const {
        return AbacusComplex(real * other.real - imag * other.imag, real * other.imag + imag * other.real);
    }

    // Multiplication operator with scalar
    AbacusComplex operator*(T scalar) const {
        return AbacusComplex(real * scalar, imag * scalar);
    }

    // In-place multiplication operator
    AbacusComplex& operator*=(const AbacusComplex& other) {
        T r = real * other.real - imag * other.imag;
        T i = real * other.imag + imag * other.real;
        real = r;
        imag = i;
        return *this;
    }

    // Division operator
    AbacusComplex operator/(const AbacusComplex& other) const {
        T denominator = other.real * other.real + other.imag * other.imag;
        return AbacusComplex((real * other.real + imag * other.imag) / denominator,
                             (imag * other.real - real * other.imag) / denominator);
    }

    // Division operator with scalar
    AbacusComplex operator/(T scalar) const {
        return AbacusComplex(real / scalar, imag / scalar);
    }

    // In-place division operator
    AbacusComplex& operator/=(const AbacusComplex& other) {
        T denominator = other.real * other.real + other.imag * other.imag;
        T r = (real * other.real + imag * other.imag) / denominator;
        T i = (imag * other.real - real * other.imag) / denominator;
        real = r;
        imag = i;
        return *this;
    }

    // Negation operator
    AbacusComplex operator-() const {
        return AbacusComplex(-real, -imag);
    }

    // Explicit cast for AbacusComplex<T> to AbacusComplex<U>,
    // e.g. for narrowing or widening conversions
    template <typename U>
    explicit operator AbacusComplex<U>() const {
        return AbacusComplex<U>(static_cast<U>(real), static_cast<U>(imag));
    }
};

// provide some std::complex overloads
template <typename T>
AbacusComplex<T> std::conj(const AbacusComplex<T>& z) {
    return AbacusComplex<T>(z.real, -z.imag);
}

template <typename T>
T std::real(const AbacusComplex<T>& z) {
    return z.real;
}

template <typename T>
T std::imag(const AbacusComplex<T>& z) {
    return z.imag;
}

// Attempt to check binary compatibility with std::complex and fftw_complex
static_assert(sizeof(AbacusComplex<double>) == sizeof(std::complex<double>), "AbacusComplex<double> is not binary compatible with std::complex<double>");
static_assert(sizeof(AbacusComplex<float>) == sizeof(std::complex<float>), "AbacusComplex<float> is not binary compatible with std::complex<float>");

// check alignment requirement compared to std::complex and fftw_complex
static_assert(alignof(AbacusComplex<double>) == alignof(std::complex<double>), "AbacusComplex<double> has different alignment requirement than std::complex<double>");
static_assert(alignof(AbacusComplex<float>) == alignof(std::complex<float>), "AbacusComplex<float> has different alignment requirement than std::complex<float>");

static_assert(sizeof(AbacusComplex<double>) == sizeof(fftw_complex), "AbacusComplex<double> is not binary compatible with fftw_complex");
static_assert(sizeof(AbacusComplex<float>) == sizeof(fftwf_complex), "AbacusComplex<float> is not binary compatible with fftwf_complex");

static_assert(alignof(AbacusComplex<double>) == alignof(fftw_complex), "AbacusComplex<double> has different alignment requirement than fftw_complex");
static_assert(alignof(AbacusComplex<float>) == alignof(fftwf_complex), "AbacusComplex<float> has different alignment requirement than fftwf_complex");

#endif // __ABACUS_COMPLEX_HH__
