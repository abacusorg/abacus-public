#ifndef INCLUDE_FILE
#define INCLUDE_FILE

#include <sys/stat.h>
#include <dirent.h>
#include <libgen.h>

#include <string.h>
#include <limits.h>     /* PATH_MAX */
#include <errno.h>

#include "config.h"

// TODO: This routine is too heavy on asserts.  Would be better
// to return an answer to the calling program and let it decide how
// to react.

// TODO: We might want some routines to concatenate path names.

#define FPL         __FILE__ , __PRETTY_FUNCTION__ , __LINE__
#define DFPL        do { fprintf(stderr, "  %s::%s::%d  ", FPL); } while(0)
#define DERRNO      do {fprintf(stderr, "  errno(%d)::%s  ", errno, strerror(errno) ); } while(0)
#define FAILERRNO   do { if(errno) { DFPL; DERRNO; assert(errno==0); } } while(0)

#ifndef assertf
#define assertf(_mytest,...) do { \
    if (!(_mytest)) { \
        fprintf(stderr,"Failed Assertion: %s\n", #_mytest); \
        fprintf(stderr, __VA_ARGS__); \
        assert(0==99); \
    }} while(0)
#endif


// This expands to a real path, but the given name must exist!
void ExpandPathName(char *foo) {
    char str[1024];
    sprintf(str,"%s",foo);
    errno = 0;
    char *retval = realpath(str,foo);
    if(errno)
        fprintf(stderr, "realpath error code %d %s\n", errno, strerror(errno));
    // TODO: This STDLOG appears needed to avoid a buffer overflow in -O3 in g++
    // on ted.
    assertf(retval!=NULL, "realpath failed on path \"%s\"\n", foo);
}

void CheckDirectoryExists(const char *fn) {
    if ( access( fn, 0 ) == 0 ) {
        struct stat status;
        stat( fn, &status );

        if (!( status.st_mode & S_IFDIR ))  {
            fprintf(stderr,"%s is a file\n",fn);
            assert(1==0);
        }
    }
    else {
        fprintf(stderr,"%s doesn't even exist\n",fn);
        assert(1==0);
    }
}

// A recursive mkdir function.
// This is semantically similar to Python's os.makedirs()
// from https://gist.github.com/JonathonReinhart/8c0d90191c38af2dcadb102c4e202950
int CreateDirectories(const char *path){
    const size_t len = strlen(path);
    char _path[PATH_MAX];
    char *p; 

    errno = 0;

    /* Copy string so its mutable */
    if (len > sizeof(_path)-1) {
        errno = ENAMETOOLONG;
        return -1; 
    }   
    strcpy(_path, path);

    /* Iterate the string */
    for (p = _path + 1; *p; p++) {
        if (*p == '/') {
            /* Temporarily truncate */
            *p = '\0';

            if (mkdir(_path, S_IRWXU) != 0) {
                if (errno != EEXIST)
                    return -1; 
            }

            *p = '/';
        }
    }   

    if (mkdir(_path, S_IRWXU) != 0) {
        if (errno != EEXIST)
            return -1;
    }   

    return 0;
}

// Recursively remove directories
// Semantically similar to Python's shutil.rmtree()
// from https://stackoverflow.com/questions/2256945/removing-a-non-empty-directory-programmatically-in-c-or-c
int RemoveDirectories(const char *path){
   DIR *d = opendir(path);
   size_t path_len = strlen(path);
   int r = -1;

   if (d)
   {
      struct dirent *p;

      r = 0;

      while (!r && (p=readdir(d)))
      {
          int r2 = -1;
          char *buf;
          size_t len;

          /* Skip the names "." and ".." as we don't want to recurse on them. */
          if (!strcmp(p->d_name, ".") || !strcmp(p->d_name, ".."))
          {
             continue;
          }

          len = path_len + strlen(p->d_name) + 2; 
          buf = (char *) malloc(len);
          assert(buf != NULL);

         struct stat statbuf;

         snprintf(buf, len, "%s/%s", path, p->d_name);

         if (!lstat(buf, &statbuf))  // stat or lstat?
         {
            if (S_ISDIR(statbuf.st_mode))
            {
               r2 = RemoveDirectories(buf);
            }
            else
            {
               r2 = unlink(buf);
            }
         }

         free(buf);

          r = r2;
      }

      closedir(d);
   }

   if (!r)
   {
      r = rmdir(path);
   }

   return r;
}

// Check if two paths point to the same file/directory
// This is semantically similar to Python's os.path.samefile()
int samefile(const char *path1, const char *path2) {
    struct stat s1, s2;
    /*assertf(stat(path1, &s1) == 0, "Stat on \"%s\" failed\n", path1);
    assertf(stat(path2, &s2) == 0, "Stat on \"%s\" failed\n", path2);*/

    int res1 = stat(path1, &s1);
    int res2 = stat(path2, &s2);

    // If one exists but the other doesn't, then they're not the same file!
    if(res1 != res2)
        return 0;

    // Neither exists... do we compare paths now?
    if(res1 != 0 && res2 != 0){
        //assertf(0, "stat failed on both \"%s\" and \"%s\"", path1, path2);
        return strcmp(path1, path2) == 0;
    }

    return (s1.st_ino == s2.st_ino) && (s1.st_dev == s2.st_dev);
}

int CreateSubDirectory(const char *path, const char *subdir) {
    // This should check whether the subdirectory exists and if not make it.
    // Return 0 if all well.

    char fn[1100];
    sprintf(fn, "%s/%s", path, subdir);

    // It's okay for mkdir to fail! Another node could have beaten it to the punch
    int ret = mkdir(fn, 0775);
    int reason = errno;
    assertf(ret == 0 || reason == EEXIST, "mkdir(\"%s\") failed for reason %s\n", fn, strerror(reason));
    assert(access(fn,0)==0);

    // Or the parent directory might not exist, flag that for debugging
    CheckDirectoryExists(path);

    return 0;
}

int CreateSymlink(const char* target, const char* link){
    int ret = symlink(target, link);
    int reason = errno;
    assertf(ret == 0 || errno == EEXIST, "symlink(\"%s\", \"%s\") failed for reason %s", target, link, strerror(reason));

    return 0;
}

int FileExists(const char *fn) {
    return (access(fn,0) == 0);
}

// Returns 0 if `fn` exists
// Returns 1 if `fn` does not exist
// Returns 2 if `fn` exists but is a directory
int CheckFileExists(const char *fn) {
    if ( access( fn, 0 ) == 0 ) {
        struct stat status;
        stat( fn, &status );

        if (status.st_mode & S_IFDIR )  {
            return 2;
        }
    }
    else {
        return 1;
    }

    return 0;
}

off_t fsize(const char *filename) {
    struct stat st;
    assertf(CheckFileExists(filename) == 0, "fsize failed to get size of file \"%s\" because it does not exist or is a dir\n", filename);

    if (stat(filename, &st) == 0)
        return st.st_size;

    return -1;
}

// Places the name of the immediate parent directory of `filename` in `dir`
// We use this for recording the IO performance in different directories/filesystems
// Special behavior: if the containing directory is "Step*", then the parent is returned.
void containing_dirname(const char *filename, char dir[1024]){
    // Believe it or not, dirname modifies its argument
    char buffer[1024];
    strncpy(buffer,filename,1024);
    char *_dir = basename(dirname(buffer));
        if(strncmp(_dir, "Step", 4) == 0){
                strncpy(buffer,filename,1024);
                _dir = basename(dirname(dirname(buffer)));
        }
    strncpy(dir, _dir, 1024);
    
    // now append a trailing slash
    int len = strlen(dir);
    assertf(len <= 1022, "Directory \"%s\" name too long!", dir);
    dir[len] = '/'; dir[len+1] = '\0';
}


// Extract the dir and name components of the path
void split_path(const char path[1024], char dir[1024], char name[1024]){
    char buf[1024];
    strncpy(buf, path, 1024);
    char *dn = dirname(buf);
    strncpy(dir, dn, 1024);

    strncpy(buf, path, 1024);
    char *bn = basename(buf);
    strncpy(name, bn, 1024);
}


int is_path_on_ramdisk(std::string path){
    const char *c_str = path.c_str();
    int res = is_path_on_ramdisk(c_str);
    return res;
}

int is_path_on_ramdisk(const char* path){
    // This Ramdisk detection via path name is not very elegant
    // But if there is a programmatic way to determine it, I haven't found it

    char str[1024];
    strncpy(str, path, 1024);
    // We've been handed a file name, but that file may not yet exist.
    // We'd prefer to only pass in the directory path, not the file itself.
    // Clip off the file name
    for (int j=strlen(str)-1; j>=0; j--)
        if (str[j]=='/') break; else str[j]='\0';

    ExpandPathName(str);
    assert(strlen(str) < 1024);
    return strncmp(str, RAMDISK_PATH, strlen(RAMDISK_PATH)) == 0;
}

int IsTrueLocalDirectory(const char*);

#endif // INCLUDE_FILE
