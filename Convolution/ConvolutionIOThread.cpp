#include "tbb/concurrent_queue.h"

class ConvIOThread {
private:
    // Queue to send blocks to the compute code
    tbb::concurrent_bounded_queue<Block *> read_queue;

    // Queue to receive blocks from the compute code
    tbb::concurrent_bounded_queue<Block *> write_queue;

    // Queue with pointers to free blocks that can be read into
    tbb::concurrent_bounded_queue<Block *> free_queue;
    
    ConvolutionParameters CP;
    
    size_t taylor_bytes_written = 0, multipole_bytes_read = 0, derivative_bytes_read = 0;
    double taylor_write_time = 0, multipole_read_time = 0, derivative_read_time = 0;
    int nblocks = 0;
    int thread_num = -1;

    // This is the main thread loop
    void ConvolutionIOThread(){
        {
        int io_core = CP.io_cores[thread_num];
        if(io_core >= 0){
            set_core_affinity(io_core);
            STDLOG(0, "IO thread %d started on core %d\n", thread_num, io_core);
        } else {
            STDLOG(0, "IO thread %d started; not bound to core\n", thread_num);
        }
        }
        
        int n_blocks_read = 0;
        int n_blocks_written = 0;
        int cpd = CP.runtime_cpd;
        int zwidth = CP.zwidth;

        assert(!free_queue.empty());
        while(n_blocks_written < nblocks){
            Block *write_buffer = NULL, *read_buffer = NULL;
            // Write anything on the write queue.
            if(!write_queue.try_pop(write_buffer))
                // If nothing on the write queue and no free blocks for read, wait for a write so read will have something to read into
                if(free_queue.empty())
                    write_queue.pop(write_buffer);
                
            if(write_buffer != NULL) {
                int zstart = n_blocks_written*CP.zwidth;
                int this_zwidth = min(zwidth, (cpd+1)/2-zstart);
                write_buffer->write(n_blocks_written*CP.zwidth, this_zwidth, thread_num);
                free_queue.push(write_buffer);
                write_buffer = NULL;
                n_blocks_written++;
            }

            // Read into any free blocks
            if(n_blocks_read < nblocks && free_queue.try_pop(read_buffer)){
                int zstart = n_blocks_read*CP.zwidth;
                int this_zwidth = min(zwidth, (cpd+1)/2-zstart);
                
                read_buffer->read_derivs(zstart, this_zwidth, thread_num);
                read_buffer->read(zstart, this_zwidth, thread_num);
                
                n_blocks_read++;
                read_queue.push(read_buffer);
            }
        }

        // Clean up and exit
        assert(n_blocks_read == n_blocks_written);
        assert(read_queue.empty());
        assert(write_queue.empty());
        assert(free_queue.size() == read_ahead);
    }

    static void *start_thread(void *iothread_obj){
        ((ConvIOThread *) iothread_obj)->ConvolutionIOThread();
        return NULL;
    }

    pthread_t conv_io_thread;
    bool thread_joined = false;
public:
    // Read N blocks ahead of the last block written
    int read_ahead;

    ConvIOThread(ConvolutionParameters &_CP, int _read_ahead, Block **blocks, int _thread_num){
        CP = _CP;
        thread_num = _thread_num;
        read_ahead = _read_ahead;
        nblocks = (int) ceil((CP.runtime_cpd+1)/2./CP.zwidth);
        STDLOG(0, "Starting IO thread %d\n", thread_num);
        
        for(int i = 0; i < read_ahead; i++)
            free_queue.push(blocks[i]);
        
        assert(pthread_create(&conv_io_thread, NULL, start_thread, this) == 0);
    }

    ~ConvIOThread(){
        if(!thread_joined)
            join();
    }

    void push(Block *&block_to_write){
        write_queue.push(block_to_write);
    }
    
    void pop(Block *&block_to_read){
        read_queue.pop(block_to_read);
    }
    
    // The IO thread will automatically exit after writing all blocks.
    // This just ensures that it did in fact finish
    void join(){
        STDLOG(0, "Stopping IO thread\n");
        assert(pthread_join(conv_io_thread, NULL) == 0);
        
        Block *block;
        while(free_queue.try_pop(block)){
            taylor_bytes_written += block->WriteTaylorBytes;
            multipole_bytes_read += block->ReadMultipoleBytes;
            derivative_bytes_read += block->ReadDerivativeBytes;
            taylor_write_time += block->WriteTaylor.Elapsed();
            multipole_read_time += block->ReadMultipoles.Elapsed();
            derivative_read_time += block->ReadDerivatives.Elapsed();
        }
        
        thread_joined = true;
    }
    
    size_t get_multipole_bytes_read(){
        assert(thread_joined);
        return multipole_bytes_read;
    }
    
    size_t get_derivative_bytes_read(){
        assert(thread_joined);
        return derivative_bytes_read;
    }
    
    size_t get_taylor_bytes_written(){
        assert(thread_joined);
        return taylor_bytes_written;
    }
    
    double get_taylor_write_time(){
        assert(thread_joined);
        return taylor_write_time;
    }
    
    double get_multipole_read_time(){
        assert(thread_joined);
        return multipole_read_time;
    }
    
    double get_deriv_read_time(){
        assert(thread_joined);
        return derivative_read_time;
    }
};
