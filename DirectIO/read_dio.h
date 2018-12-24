class ReadDirect { 
public:
    ReadDirect(int isramdisk, size_t buffersize) {
        alignedbytes = 0;
        alignedbuffer = NULL;
        ramdiskflag = isramdisk;
        if(!isramdisk) {
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
    //       if bytesleft>4069 use direct read 
    //       finally fread if anything else remains 
    void BlockingRead(char *fn, char *x, size_t length, off_t fileoffsetbytes);

    // Passing the ramdisk flag explicitly will override the global "ramdiskflag"
    void BlockingRead(char *fn, char *x, size_t length, off_t fileoffsetbytes, int ramdisk);

private:

    size_t  alignedbytes;
    char   *alignedbuffer;

    int rdopenflags(void);
    int rdopenfd(char *fn, int readflags);

    void BlockingDirectReadAligned(char *fn, char *x, size_t length, off_t fileoffsetbytes);
    size_t min(size_t a, size_t b) { if(a<b) return a; return b; }

    void BlockingReadDirect(char *fn, char *x, size_t length, off_t fileoffsetbytes);
    void Blockingfread(char *fn, char *x, size_t length, off_t fileoffsetbytes);

    int ramdiskflag;
};
