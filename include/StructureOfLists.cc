// Easy to use 3 component float structure of lists

#ifndef __STRUCTUREOFLISTS_CC
#define __STRUCTUREOFLISTS_CC

template<typename T> class List3
{
    public:
        List3(){ }
        List3(uint64 count){
            N = count;
            size_t size = sizeof(T)*N;
            
            assert(posix_memalign((void **) &X, 4096, size) == 0);
            assert(posix_memalign((void **) &Y, 4096, size) == 0);
            assert(posix_memalign((void **) &Z, 4096, size) == 0);
        }
        ~List3(){
            if(N == 0)
                return;
            free(X);
            free(Y);
            free(Z);
        }

        T x(uint64 n){
            return X[n];
        }

        T y(uint64 n)
        {
            return Y[n];
        }

        T z(uint64 n){
            return Z[n];
        }

        void SetZero(uint64 start, uint64 count)
        {
            memset(&(X[start]), 0, sizeof(T) * count);
            memset(&(Y[start]), 0, sizeof(T) * count);
            memset(&(Z[start]), 0, sizeof(T) * count);
        }
 
        uint64 N = 0;
        T* X = NULL;
        T* Y = NULL;
        T* Z = NULL;
};

#endif
