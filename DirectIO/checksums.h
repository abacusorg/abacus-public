/*

This is the Abacus object-oriented interface to the checksum
functions defined in crc32_fast.

*/

#include "crc32_fast.cpp"

class CRC32 {
private:
    uint32 partial_crc;

public:
	size_t size;
	
    CRC32(){
        partial_crc = CRC32_FAST_SEED;
        size = 0;
    }

    void ingest(const void *data, size_t len){
    	partial_crc = crc32_fast_partial(data, len, partial_crc);
    	size += len;
    }

    uint32 finalize(){
    	return crc32_fast_finalize(size, partial_crc);
    }
};
