#ifndef INCLUDE_FILE
#define INCLUDE_FILE

#include <sys/stat.h>
#include <libgen.h>

// TODO: This routine is too heavy on asserts.  Would be better
// to return an answer to the calling program and let it decide how
// to react.

// TODO: We might want some routines to concatenate path names.

#define FPL         __FILE__ , __PRETTY_FUNCTION__ , __LINE__
#define DFPL        do { fprintf(stderr, "  %s::%s::%d  ", FPL); } while(0)
#define DERRNO      do {fprintf(stderr, "  errno(%d)::%s  ", errno, strerror(errno) ); } while(0)
#define FAILERRNO   do { if(errno) { DFPL; DERRNO; assert(errno==0); } } while(0)

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
    assert(retval!=NULL);
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

int CreateSubDirectory(const char *path, const char *subdir) {
    // This should check whether the subdirectory exists and if not make it.
    // Return 0 if all well.
    CheckDirectoryExists(path);		// Does the parent exist?
    char fn[1100];
    sprintf(fn, "%s/%s", path, subdir);
    if (access(fn,0)==0) {
	// Subdir name already exists
        struct stat status;
	stat(fn, &status);
	if (!(status.st_mode & S_IFDIR)) {
	    fprintf(stderr,"%s exists but is a file\n", fn);
	    assert(1==0);
	} else {
	    // Subdirectory already exists
	    return 0;
	}
    } else {
        // Subdir doesn't exist
	mkdir(fn, 0775);
	assert(access(fn,0)==0);
	return 0;
    }
}

int FileExists(const char *fn) {
    return (access(fn,0) == 0);
}

void CheckFileExists(const char *fn) {
    if ( access( fn, 0 ) == 0 ) {
        struct stat status;
        stat( fn, &status );

        if (status.st_mode & S_IFDIR )  {
            fprintf(stderr,"%s is a directory\n",fn);
            assert(1==0);
        }
    }
    else {
        fprintf(stderr,"%s doesn't even exist\n",fn);
        assert(1==0);
    }
}

off_t fsize(const char *filename) {
    struct stat st;
    CheckFileExists(filename);

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
    if (len >= 1022){
        fprintf(stderr,"Directory \"%s\" name too long!", dir);
    }
    dir[len] = '/'; dir[len+1] = '\0';
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

#endif // INCLUDE_FILE
