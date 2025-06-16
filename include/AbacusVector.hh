// Copyright 2012-2025 The Abacus Developers
// SPDX-License-Identifier: GPL-3.0-or-later

/*

AbacusVector.hh

It's often convenient to have an array container that keeps track of its own
length and capacity.  Normally one would use std::vector for this, but
this has two drawbacks:
    
    1) std::vector is resizeable. In high-performance code, we want to be very sure
       that no hidden reallocations are occuring.
    2) std::vector does bounds checking on insertion.

We would like to be able to disable, or maybe toggle, these features.

*/

#include <utility>

template <typename T, bool bounds_check = false>
class alignas(CACHE_LINE_SIZE) AbacusVector {
    // At least resize() and clear() don't support non-POD for the moment
    static_assert(std::is_pod<T>::value, "T must be POD");

private:

    size_t size_ = 0;
    size_t capacity_;

    T *arr = nullptr;

public:
    AbacusVector() : capacity_(0) { }

    AbacusVector(size_t capacity){
        this->reserve(capacity);
    }

    ~AbacusVector() {
        delete[] this->arr;
    }

    AbacusVector(const AbacusVector& other) = delete;
    AbacusVector& operator=(const AbacusVector& other) = delete;
    AbacusVector(AbacusVector&& other) noexcept {
        this->arr = other.arr;
        this->size_ = other.size_;
        this->capacity_ = other.capacity_;

        other.arr = nullptr;
        other.size_ = 0;
        other.capacity_ = 0;
    };
    AbacusVector& operator=(AbacusVector&& other) noexcept = delete;

    void reserve(size_t capacity){
        assert(this->arr == nullptr);  // TODO

        this->arr = new T[capacity];

        this->capacity_ = capacity;
    }

    void push_back(const T& x){
        this->arr[this->size_++] = x;

        if constexpr (bounds_check) this->validate();
    }

    void push_back(T&& x){
        this->arr[this->size_++] = std::move(x);

        if constexpr (bounds_check) this->validate();
    }

    template<typename... Args>
    T& emplace_back(Args&&... args){
        new(&this->arr[this->size_++]) T{std::forward<Args>(args)...};

        if constexpr (bounds_check) this->validate();

        return this->arr[this->size_-1];
    }

    T& operator[](size_t i) {
        return arr[i];
    }

    const T& operator[](size_t i) const {
        return arr[i];
    }

    size_t size() const {
        return this->size_;
    }

    size_t capacity() const {
        return this->capacity_;
    }

    void resize(size_t newsize){
        this->size_ = newsize;
        if constexpr (bounds_check) this->validate();
    }

    void clear(){
        this->resize(0);
    }

    bool empty() const {
        return this->size_ == 0;
    }

    // TODO: do these need different properties for faster sorting...?
    T* begin() const {
        return this->arr;
    }

    T* end() const {
        return this->arr + this->size_;
    }

    void validate() const {
        assert(this->size_ <= this->capacity_);
    }
};
