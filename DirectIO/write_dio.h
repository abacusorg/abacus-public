#ifndef __WRITE_DIO_H__
#define __WRITE_DIO_H__

class WriteDirect {
public:
    WriteDirect(int no_directio, size_t buffersize) {
        alignedbytes = 0;
        alignedbuffer = NULL;
        allow_directio = !no_directio;

        if(allow_directio){
            alignedbytes = buffersize;
            int rv = posix_memalign( (void **) (&alignedbuffer), 4096, buffersize);
            assert(rv==0);
            assert(alignedbuffer!=NULL);
        }
    }
    ~WriteDirect(void) { if(alignedbuffer!=NULL) free(alignedbuffer); }

    void BlockingAppend( const fs::path &fn, char *x, size_t length );
    void BlockingAppend( const fs::path &fn, char *x, size_t length, int no_dio);
    void BlockingAppend( FILE *f, char *x, size_t length);

private:

    size_t  alignedbytes;
    char *alignedbuffer;

    int wropenflags(void);
    int wropenfd(const fs::path &fn, int writeflags);

    void BlockingWrite_Append_Aligned(const fs::path &fn, char *x, size_t length);
    size_t min(size_t a, size_t b) { if(a<b) return a; return b; }

    void BlockingAppendDirect( const fs::path &fn, char *x, size_t length);
    void BlockingAppendfwrite( const fs::path &fn, char *x, size_t length);
    void BlockingAppendPointer( FILE *f, char *x, size_t length);

    int allow_directio;
};

#endif
