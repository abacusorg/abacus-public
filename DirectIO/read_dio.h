class ReadDirect { 
public:
    ReadDirect(int no_directio, size_t buffersize) {
        alignedbytes = 0;
        alignedbuffer = NULL;
        allow_directio = !no_directio;
        
        if(allow_directio){
            alignedbytes = buffersize;
            int rv = posix_memalign((void **) (&alignedbuffer), 4096, buffersize);
            assert(rv==0);
            assert(alignedbuffer!=NULL);
        }
    }
    ~ReadDirect(void) { if(alignedbuffer!=NULL) free(alignedbuffer); }

    // Use blockingfread if ramdisk -- can't use O_DIRECT 
    // direct algorithm is :
    //       firstly align file offset by freading until offset is aligned 
    //       if bytesleft>=4096 use direct read 
    //       finally fread if anything else remains 
    void BlockingRead(const fs::path &fn, char *x, size_t length, size_t fileoffsetbytes);

    void BlockingRead(const fs::path &fn, char *x, size_t length, size_t fileoffsetbytes, int no_dio);

private:

    size_t  alignedbytes;
    char   *alignedbuffer;

    int rdopenflags(void);
    int rdopenfd(const fs::path &fn, int readflags);

    void BlockingDirectReadAligned(const fs::path &fn, char *x, size_t length, size_t fileoffsetbytes);
    size_t min(size_t a, size_t b) { if(a<b) return a; return b; }

    void BlockingReadDirect(const fs::path &fn, char *x, size_t length, size_t fileoffsetbytes);
    void Blockingfread(const fs::path &fn, char *x, size_t length, size_t fileoffsetbytes);

    int allow_directio;
};
