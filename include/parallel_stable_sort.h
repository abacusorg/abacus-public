#ifndef __PARALLEL_STABLE_SORT_H
#define __PARALLEL_STABLE_SORT_H

/*
  Copyright (C) 2014 Intel Corporation
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions
  are met:

  * Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.
  * Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in
    the documentation and/or other materials provided with the
    distribution.
  * Neither the name of Intel Corporation nor the names of its
    contributors may be used to endorse or promote products derived
    from this software without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
  OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
  AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY
  WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
  POSSIBILITY OF SUCH DAMAGE.
*/
#include <iterator>
#include <algorithm>
#include <tbb/task.h>

/*
  Copyright (C) 2014 Intel Corporation
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions
  are met:

  * Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.
  * Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in
    the documentation and/or other materials provided with the
    distribution.
  * Neither the name of Intel Corporation nor the names of its
    contributors may be used to endorse or promote products derived
    from this software without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
  OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
  AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY
  WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
  POSSIBILITY OF SUCH DAMAGE.
*/
#include <utility>
#include <iterator>

namespace pss {

namespace internal {

//! Destroy sequence [xs,xe)
template<class RandomAccessIterator>
void serial_destroy( RandomAccessIterator zs, RandomAccessIterator ze ) {
    typedef typename std::iterator_traits<RandomAccessIterator>::value_type T;
    while( zs!=ze ) {
        --ze;
        (*ze).~T();
    }
}

//! Merge sequences [xs,xe) and [ys,ye) to output sequence [zs,(xe-xs)+(ye-ys)), using std::move
template<class RandomAccessIterator1, class RandomAccessIterator2, class RandomAccessIterator3, class Compare>
void serial_move_merge( RandomAccessIterator1 xs, RandomAccessIterator1 xe, RandomAccessIterator2 ys, RandomAccessIterator2 ye, RandomAccessIterator3 zs, Compare comp ) {
    if( xs!=xe ) {
        if( ys!=ye )
            for(;;)
                if( comp(*ys,*xs) ) {
                    *zs = std::move(*ys);
                    ++zs;
                    if( ++ys==ye ) break;
                } else {
                    *zs = std::move(*xs);
                    ++zs;
                    if( ++xs==xe ) goto movey;
                }
        ys = xs;
        ye = xe;
    }
movey:
    std::move( ys, ye, zs );
}

template<typename RandomAccessIterator1, typename RandomAccessIterator2, typename Compare>
void stable_sort_base_case(  RandomAccessIterator1 xs, RandomAccessIterator1 xe, RandomAccessIterator2 zs, int inplace, Compare comp) {
    std::stable_sort( xs, xe, comp );
    if( inplace!=2 ) {
        RandomAccessIterator2 ze = zs + (xe-xs);
        typedef typename std::iterator_traits<RandomAccessIterator2>::value_type T;
        if( inplace )
            // Initialize the temporary buffer
            for( ; zs<ze; ++zs )
                new(&*zs) T;
        else
            // Initialize the temporary buffer and move keys to it.
            for( ; zs<ze; ++xs, ++zs )
                new(&*zs) T(std::move(*xs));
    }
}

//! Raw memory buffer with automatic cleanup.
class raw_buffer {
    void* ptr;
public:
    //! Try to obtain buffer of given size.
    raw_buffer( size_t bytes ) : ptr( operator new(bytes,std::nothrow) ) {}
    //! True if buffer was successfully obtained, zero otherwise.
    operator bool() const {return ptr;}
    //! Return pointer to buffer, or  NULL if buffer could not be obtained.
    void* get() const {return ptr;}
    //! Destroy buffer
    ~raw_buffer() {operator delete(ptr);}
};

} // namespace internal

} // namespace pss


namespace pss {

namespace internal {

template<typename RandomAccessIterator1, typename RandomAccessIterator2, typename RandomAccessIterator3, typename Compare>
class merge_task: public tbb::task {
    tbb::task* execute();
    RandomAccessIterator1 xs, xe;
    RandomAccessIterator2 ys, ye;
    RandomAccessIterator3 zs;
    Compare comp;
    bool destroy;
public:
    merge_task( RandomAccessIterator1 xs_, RandomAccessIterator1 xe_, RandomAccessIterator2 ys_, RandomAccessIterator2 ye_, RandomAccessIterator3 zs_, bool destroy_, Compare comp_ ) :
        xs(xs_), xe(xe_), ys(ys_), ye(ye_), zs(zs_), destroy(destroy_), comp(comp_)
    {}
};

template<typename RandomAccessIterator1, typename RandomAccessIterator2, typename RandomAccessIterator3, typename Compare>
tbb::task* merge_task<RandomAccessIterator1,RandomAccessIterator2,RandomAccessIterator3,Compare>::execute() {
    const size_t MERGE_CUT_OFF = 2000;
    auto n = (xe-xs) + (ye-ys);
    if( n <= MERGE_CUT_OFF ) {
        serial_move_merge( xs, xe, ys, ye, zs, comp );
        if( destroy ) {
            serial_destroy(xs,xe);
            serial_destroy(ys,ye);
        }
        return NULL;
    } else {
        RandomAccessIterator1 xm;
        RandomAccessIterator2 ym;
        if( xe-xs < ye-ys  ) {
            ym = ys+(ye-ys)/2;
            xm = std::upper_bound(xs,xe,*ym,comp);
        } else {
            xm = xs+(xe-xs)/2;
            ym = std::lower_bound(ys,ye,*xm,comp);
        }
        RandomAccessIterator3 zm = zs + ((xm-xs) + (ym-ys));
        tbb::task* right = new( allocate_additional_child_of(*parent()) ) merge_task( xm, xe, ym, ye, zm, destroy, comp );
        spawn(*right);
        recycle_as_continuation();
        xe = xm;
        ye = ym;
        return this;
    }
}

template<typename RandomAccessIterator1, typename RandomAccessIterator2, typename Compare>
class stable_sort_task: public tbb::task {
    tbb::task* execute();
    RandomAccessIterator1 xs, xe;
    RandomAccessIterator2 zs;
    Compare comp;
    signed char inplace;
public:
    stable_sort_task(RandomAccessIterator1 xs_, RandomAccessIterator1 xe_, RandomAccessIterator2 zs_, int inplace_, Compare comp_ ) : 
        xs(xs_), xe(xe_), zs(zs_), inplace(inplace_), comp(comp_) 
    {}
};

template<typename RandomAccessIterator1, typename RandomAccessIterator2, typename Compare>
tbb::task* stable_sort_task<RandomAccessIterator1, RandomAccessIterator2, Compare>::execute() {
    const size_t SORT_CUT_OFF = 500;
    if (xe - xs <= SORT_CUT_OFF) {
        stable_sort_base_case(xs, xe, zs, inplace, comp); 
        return NULL;
    } else {
        RandomAccessIterator1 xm = xs + (xe - xs) / 2;
        RandomAccessIterator2 zm = zs + (xm - xs);
        RandomAccessIterator2 ze = zs + (xe - xs);
        task* m;
        if (inplace)
            m = new (allocate_continuation()) merge_task<RandomAccessIterator2,RandomAccessIterator2,RandomAccessIterator1,Compare>(zs, zm, zm, ze, xs, inplace==2, comp);
        else
            m = new (allocate_continuation()) merge_task<RandomAccessIterator1,RandomAccessIterator1,RandomAccessIterator2,Compare>(xs, xm, xm, xe, zs, false, comp);
        m->set_ref_count(2);
        task* right = new(m->allocate_child()) stable_sort_task(xm,xe,zm,!inplace, comp);
        spawn(*right);
        recycle_as_child_of(*m); 
        xe=xm;
        inplace=!inplace;
        return this;
    }
}

} // namespace internal 

template<typename RandomAccessIterator, typename Compare>
void parallel_stable_sort( RandomAccessIterator xs, RandomAccessIterator xe, Compare comp ) {
    typedef typename std::iterator_traits<RandomAccessIterator>::value_type T;
    if( internal::raw_buffer z = internal::raw_buffer( sizeof(T)*(xe-xs) ) ) {
        using tbb::task;
        typedef typename std::iterator_traits<RandomAccessIterator>::value_type T;
        internal::raw_buffer buf( sizeof(T)*(xe-xs) );
        task::spawn_root_and_wait(*new( task::allocate_root() ) internal::stable_sort_task<RandomAccessIterator,T*,Compare>( xs, xe, (T*)buf.get(), 2, comp ));
    } else
        // Not enough memory available - fall back on serial sort
        std::stable_sort( xs, xe, comp );
}
    
//! Wrapper for sorting with default comparator.
template<class RandomAccessIterator>
void parallel_stable_sort( RandomAccessIterator xs, RandomAccessIterator xe ) {
    typedef typename std::iterator_traits<RandomAccessIterator>::value_type T;
    parallel_stable_sort( xs, xe, std::less<T>() );
}

} // namespace pss

#endif