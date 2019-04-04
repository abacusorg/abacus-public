class WriteDirect { 
public:
    WriteDirect(int isramdisk, size_t buffersize) {
        alignedbytes = 0;
        alignedbuffer = NULL;
        ramdiskflag = isramdisk;

        if(!ramdiskflag){
            alignedbytes = buffersize;
            int rv = posix_memalign( (void **) (&alignedbuffer), 4096, buffersize);
            assert(rv==0);
            assert(alignedbuffer!=NULL);
        }
    }
    ~WriteDirect(void) { if(alignedbuffer!=NULL) free(alignedbuffer); }

    void BlockingAppend( char *fn, char *x, size_t length );
    void BlockingAppend( char *fn, char *x, size_t length, int no_dio);

private:

    size_t  alignedbytes;
    char *alignedbuffer;

    int wropenflags(void);
    int wropenfd(char *fn, int writeflags);

    void BlockingWrite_Append_Aligned(char *fn, char *x, size_t length);
    size_t min(size_t a, size_t b) { if(a<b) return a; return b; }

    void BlockingAppendDirect( char *fn, char *x, size_t length);
    void BlockingAppendfwrite( char *fn, char *x, size_t length);

    int ramdiskflag;
};
