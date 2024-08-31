#ifndef INCLUDE_FILE
#define INCLUDE_FILE

#include <sys/stat.h>
#include <dirent.h>
#include <libgen.h>

#include <string.h>
#include <limits.h>     /* PATH_MAX */
#include <errno.h>

#include "config.h"

#define FPL         __FILE__ , __PRETTY_FUNCTION__ , __LINE__
#define DFPL        do { fmt::print(stderr, "  {:s}::{:s}::{:d}  ", FPL); } while(0)
#define DERRNO      do {fmt::print(stderr, "  errno({:d})::{:s}  ", errno, strerror(errno) ); } while(0)
#define FAILERRNO   do { if(errno) { DFPL; DERRNO; assert(errno==0); } } while(0)

#ifndef assertf
#define assertf(_mytest,...) do { \
    if (!(_mytest)) { \
        fprintf(stderr,"Failed Assertion: %s\n", #_mytest); \
        fprintf(stderr, __VA_ARGS__); \
        assert(0==99); \
    }} while(0)
#endif


// Check if two paths point to the same file/directory
// This is semantically similar to Python's os.path.samefile()
int samefile(const fs::path &path1, const fs::path &path2) {
    std::error_code ec;  // ignore errors
    if (fs::equivalent(path1, path2, ec)) return 1;

    // Maybe they don't exist but are the same string?
    return path1 == path2;
}

int is_path_on_ramdisk(const fs::path &path){
    // This Ramdisk detection via path name is not very elegant
    // But if there is a programmatic way to determine it, I haven't found it
    return strncmp(fs::absolute(path).c_str(), RAMDISK_PATH, strlen(RAMDISK_PATH)) == 0;
}

int IsTrueLocalDirectory(const fs::path &path);

#endif // INCLUDE_FILE
