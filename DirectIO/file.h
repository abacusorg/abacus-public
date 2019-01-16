#define FPL         __FILE__ , __PRETTY_FUNCTION__ , __LINE__
#define DFPL        do { fprintf(stderr, "  %s::%s::%d  ", FPL); } while(0)
#define DERRNO      do {fprintf(stderr, "  errno(%d)::%s  ", errno, strerror(errno) ); } while(0)
#define FAILERRNO   do { if(errno) { DFPL; DERRNO; assert(errno==0); } } while(0)

void ExpandPathName(char *foo);
int DirectoryExists(const char *fn);
void CheckDirectoryExists(const char *fn);
int FileExists(const char *fn);
void CheckFileExists(const char *fn);
off_t fsize(const char *filename);

int is_path_on_ramdisk(std::string path);
int is_path_on_ramdisk(const char* path);
