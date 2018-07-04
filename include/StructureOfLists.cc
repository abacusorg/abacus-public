// Easy to use 3 component float structure of lists

#ifndef __STRUCTUREOFLISTS_CC
#define __STRUCTUREOFLISTS_CC

#include <stdint.h>
#include <stdlib.h>

// TODO: we could probably make template specializations of this instead of an owndata field
template<typename T> class List3
{
    public:
        // This constructor simply makes an empty wrapper object
        // that will not try to free memory on destruction
        List3(){
            N = 0;
            X = NULL, Y = NULL, Z = NULL;
            owndata = false;
        }
        
        // This constructor allocates memory and frees it on destruction
        List3(uint64_t count){
            N = count;
            size_t size = sizeof(T)*N;
            
            assert(posix_memalign((void **) &X, 4096, size) == 0);
            assert(posix_memalign((void **) &Y, 4096, size) == 0);
            assert(posix_memalign((void **) &Z, 4096, size) == 0);
            owndata = true;
        }
        ~List3(){
            if(!owndata)
                return;
            free(X);
            free(Y);
            free(Z);
        }

        T x(uint64_t n){
            return X[n];
        }

        T y(uint64_t n)
        {
            return Y[n];
        }

        T z(uint64_t n){
            return Z[n];
        }

        void SetZero(uint64_t start, uint64_t count)
        {
            memset(&(X[start]), 0, sizeof(T) * count);
            memset(&(Y[start]), 0, sizeof(T) * count);
            memset(&(Z[start]), 0, sizeof(T) * count);
        }
 
        uint64_t N;
        T *X, *Y, *Z;
        bool owndata;
};

#endif
