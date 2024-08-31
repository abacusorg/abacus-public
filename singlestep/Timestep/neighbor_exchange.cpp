/* file: neighbor_exchange.cpp
 *
 * The "neighbor exchange" is a 2D parallel concept where ghost
 * particles are updated across the z-split of a slab.
 */

#ifdef PARALLEL
// Partition function
inline bool is_in_range(ilstruct *particle, int slab, int left, int len){
    // returns true if particle is in [left,left+len)
    // assumes left is in [0,CPD]
    if(particle->newslab != slab) return 0;

    int cellz = particle->global_cellz();

    int dist = cellz - left;
    if (dist < 0) dist += P.cpd;

    return dist < len;
}

/* NeighborExchange happens when a slab and both adjacent slabs have Drifted,
   and does the following for each neighbor rank:
 * 1. Partition the IL to find particles the neighbor needs (for its primary or ghost)
 * 2. Copy them to a buffer
 * 3. Initiate an MPI_Isend
 * 4. Partition the tail of the IL, this time for particles that are *only* in the neighbor's primary
 * 5. "Delete" them
 * 
 * And then later...
 * 6. Check if we have an incoming MPI message
 * 7. Allocate space to receive the message
 * 8. Start the MPI_Irecv
 * 
 * And then later...
 * 9. Check if the receive is done
 * 10. Push the buffer onto the insert list
 * 11. Free the receive buffer and release the MPI handle
 *
 * And sometime after the send...
 * 12. Check if the send is complete
 * 13. Free the send buffer and release the MPI handle
 */

class NeighborExchanger {
/* Represents a z exchange with one neighbor. */
private:
    int slab;
    int zneigh;  // neighbor z rank
    int right;   // sending the left or right side of our domain

    int send_status = 0;  // 0 = not sent; 1 = sending; 2 = done
    ilstruct *sendbuf = NULL;  // staging the MPI_Isend
    size_t sendelem = 0;    // number of ilstructs

    int recv_status = 0;    // 0 = waiting; 1 = receiving; 2 = done
    ilstruct *recvbuf = NULL;  // staging the MPI_Irecv
    size_t recvelem = 0;

    MPI_Request send_handle = MPI_REQUEST_NULL;
    MPI_Request recv_handle = MPI_REQUEST_NULL;

    // If MPI_size_z=2, then the right and left neighbor will be the same.
    // So we'll need separate tag "namespaces".
    static const int RIGHT_TAG_OFFSET = 100000;
    static const int LEFT_TAG_OFFSET = 200000;

public:
    NeighborExchanger(int slab, int zneigh, int right) :
        slab(slab), zneigh(zneigh), right(right) { }

    
    void send(){
        // 1. partition
        ilstruct *start = partition(&sendelem);

        // 2. copy
        copy_to_sendbuf(start, sendelem);

        // 3. MPI_Isend
        launch_mpi_send();

        // 4. tail partition
        size_t removed = tail_partition(start, sendelem);

        // 5. delete
        IL->ShrinkMAL(IL->length - removed);
    }

    int try_free_sendbuf(){
        int ret = 0;
        if(send_status == 1){
            int done = 0;
            // 12. done send?
            MPI_Test(&send_handle, &done, MPI_STATUS_IGNORE);
            if (done) {
                STDLOG(1, "Done ilstruct send to zrank {:d} for slab {:d}\n", zneigh, slab);

                // 13. free
                free(sendbuf);
                sendbuf = NULL;
                send_status = 2;
                ret = 1;
            }
        }
        return ret;
    }

    int done_send(){
        return send_status == 2;
    }

private:
    ilstruct *partition(size_t *nhigh){
        IL->CollectGaps();
        
        uint64 mid;
        if(right){
            // A "right" ghost must be closer to our right primary domain edge than the left edge.
            // In the case of ties, the right ghost wins.
            int len = (P.cpd - node_z_size + 1)/2;

            // [0..mid) are not in slab, [mid..length) are in slab
            mid = ParallelPartition(IL->list, IL->length,
                is_in_range,
                slab,
                node_z_start + node_z_size - MERGE_GHOST_RADIUS,
                len + MERGE_GHOST_RADIUS
                );
        }
        else {
            int len = (P.cpd - node_z_size)/2;

            mid = ParallelPartition(IL->list, IL->length,
                is_in_range,
                slab,
                CP->WrapSlab(node_z_start - len),
                len + MERGE_GHOST_RADIUS
                );
        }
        
        *nhigh = (size_t) (IL->length - mid);
        return IL->list + mid;
    }


    size_t tail_partition(ilstruct *start, size_t nelem){
        /* Normal steps will not need to execute this.
           This only does anything when the IL has particles that are so
           far away from our domain that we will not need them as ghosts.
           Since we can only drift 1 cell under normal conditions, and
           MERGE_GHOST_RADIUS is usually 2, this usually cannot happen.
           Common exceptions are ingestion of ICs (which have some
           disorder), and the displacement-flipping of 2LPT.
        */

        if(WriteState.FullStepNumber != 0 &&
            WriteState.LPTStepNumber == 0 &&
            MERGE_GHOST_RADIUS > 0)
                return 0;

        // No need to collect gaps; still contiguous from partition()

        uint64 mid;
        if(right){
            int len = (P.cpd - node_z_size + 1)/2;

            // [0..mid) are not in slab, [mid..length) are in slab
            mid = ParallelPartition(start, nelem,
                is_in_range,
                slab,
                CP->WrapSlab(node_z_start + node_z_size + MERGE_GHOST_RADIUS),
                len - MERGE_GHOST_RADIUS
                );
        }
        else {
            int len = (P.cpd - node_z_size)/2;

            mid = ParallelPartition(start, nelem,
                is_in_range,
                slab,
                CP->WrapSlab(node_z_start - len),
                len - MERGE_GHOST_RADIUS
                );
        }

        size_t nhigh = (size_t) nelem - mid;
        STDLOG(1,"Removing {:d} particles that are strictly neighbor's primary\n", nhigh);
        return nhigh;
    }

    void copy_to_sendbuf(ilstruct *start, size_t nelem){
        size_t sz = sizeof(ilstruct)*nelem;
        posix_memalign((void **) &(this->sendbuf), PAGE_SIZE, sz);
        assertf(this->sendbuf != NULL,
            "Failed neighbor exchange sendbuf alloc of {:d} bytes\n", sz);

        #pragma omp parallel for schedule(static)
        for(size_t i = 0; i < nelem; i++){
            this->sendbuf[i] = start[i];
        }
    }

    void launch_mpi_send(){
        assertf(send_status == 0,
            "send() called but exchange already sent to zrank {:d}, slab {:d}?\n", zneigh, slab);

        assertf(sendelem < INT32_MAX,
            "2D neighbor exchanger trying to send {:d} particles, which overflows 32-bit int\n",
            sendelem);

        STDLOG(1, "Sending {:d} ilstructs ({:.4g} MB) to zrank {:d} for slab {:d}\n", sendelem, sendelem*sizeof(ilstruct)/1e6, zneigh, slab);
        // each pair of ranks has only one send and one recv per slab
        int tag = right ? (slab + RIGHT_TAG_OFFSET) : (slab + LEFT_TAG_OFFSET);
        MPI_Isend((const void *) sendbuf, (int) sendelem, MPI_ilstruct, zneigh, tag, comm_1d_z, &send_handle);

        send_status = 1;  // sending
    }

    void push_buffer_to_IL(){
        // Bypass IL->Push, which builds new ilstructs.
        // We just need to do a raw copy!
        IL->CollectGaps();
        uint64 oldlen = IL->length;
        IL->GrowMAL(oldlen + recvelem);

        #pragma omp parallel for schedule(static)
        for(uint64 i = 0; i < recvelem; i++){
            // Adjust sorting key to be relative to the node-local ghost_z_start.
            // TODO: could probably emplace, copying all fields except zlocal
            ilstruct p = recvbuf[i];
            p.setzlocal( CP->WrapSlab(p.local_cellz() - (node_z_start - all_node_z_start[zneigh])) );
            IL->list[oldlen + i] = p;
        }
    }

public:
    int try_receive(){
        int ret = 0;
        if (send_status > 0      // started send
            && recv_status == 1  // receiving
            ){
            // Not allowed to import recvbuf until we send! Otherwise we will send things we received

            int done = 0;
            // 9. done receive?
            MPI_Test(&recv_handle, &done, MPI_STATUS_IGNORE);
            if (done) {
                STDLOG(1, "Done ilstruct receive from zrank {:d} for slab {:d}\n", zneigh, slab);
                // 10. push to IL
                push_buffer_to_IL();

                // 11. free
                free(recvbuf);
                recvbuf = NULL;
                recv_status = 2;
                ret = 1;
            }
            return ret;
        }

        if (recv_status == 2)
            return ret;  // all done!

        if (recv_status == 0){
            // 6. check for incoming
            MPI_Status status;
            // if we sent to the right, receive from the left
            int tag = right ? (slab + LEFT_TAG_OFFSET) : (slab + RIGHT_TAG_OFFSET);
            int ready = 0;
            MPI_Iprobe(zneigh, tag, comm_1d_z, &ready, &status);
            if(ready){
                // 7. alloc recvbuf
                int _recvelem32;
                MPI_Get_count(&status, MPI_ilstruct, &_recvelem32);
                recvelem = (size_t) _recvelem32;
                size_t sz = sizeof(ilstruct)*recvelem;
                posix_memalign((void **) &recvbuf, PAGE_SIZE, sz);
                assertf(recvbuf != NULL, "Failed neighbor exchange recvbuf alloc of {:d} bytes\n", sz);
                
                // 8. start recv
                STDLOG(1, "Receiving {:d} ilstructs ({:.4g} MB) from zrank {:d} for slab {:d}\n",
                    recvelem, recvelem*sizeof(ilstruct)/1e6, zneigh, slab);
                MPI_Irecv(recvbuf, (int) recvelem, MPI_ilstruct, zneigh, tag, comm_1d_z, &recv_handle);
                recv_status = 1;
                ret = 1;
            }
            return ret;
        }

        assertf(recv_status == 1, "Unknown recv_status {:d}\n", recv_status);
        return ret;
    }

    int done_receive(){
        return recv_status == 2;
    }

    void force_done(){
        assertf(recv_status == 0, "Tried to force done an active Neighbor Exchange?");
        assertf(send_status == 0, "Tried to force done an active Neighbor Exchange?");

        recv_status = 2;
        send_status = 2;
    }
};


NeighborExchanger **left_exchanger;
NeighborExchanger **right_exchanger;
int neighbor_exchange_is_noop;

// Lightweight setup of data structures
void SetupNeighborExchange(int first, int nslab){
    if(MPI_size_z == 1){
        neighbor_exchange_is_noop = 1;
        STDLOG(1,"Skipping Neighbor Exchange setup\n");
        return;
    }

    STDLOG(1,"Setting up Neighbor Exchange to do exchanges on slabs [{:d},{:d})\n",
        first, CP->WrapSlab(first+nslab));

    left_exchanger = new NeighborExchanger*[P.cpd]();
    right_exchanger = new NeighborExchanger*[P.cpd]();

    int leftz = MPI_rank_z - 1;
    if(leftz < 0) leftz += MPI_size_z;
    
    int rightz = MPI_rank_z + 1;
    if(rightz >= MPI_size_z) rightz -= MPI_size_z;

    for(int i = first; i < first+nslab; i++){
        int iw = CP->WrapSlab(i);
        left_exchanger[iw] = new NeighborExchanger(iw, leftz, 0);
        right_exchanger[iw] = new NeighborExchanger(iw, rightz, 1);
    }
    // slabs not on this node have NULL exchangers
}


// Called before Finish, after incoming Drifts
void DoNeighborSend(int slab){
    if(neighbor_exchange_is_noop)
        return;

    STDLOG(1,"Triggering Neighbor Exchange sends on slab {:d}\n", slab);

    assertf(left_exchanger[slab] != NULL,
        "Slab {:d} not set up to do neighbor send?\n", slab);
    assertf(right_exchanger[slab] != NULL,
        "Slab {:d} not set up to do neighbor send?\n", slab);

    left_exchanger[slab]->send();
    right_exchanger[slab]->send();
}

// Called repeatedly in the timestep loop
int AttemptNeighborReceive(int first, int receive_ahead){
    // This will typically be called with
    // first = (Finish.last_slab_executed - FINISH_WAIT_RADIUS), receive_ahead ~ 3
    // It's also safe to call with (0,cpd) if one always wants to receive everything
    int ret = 0;

    if(neighbor_exchange_is_noop){
        return ret;  // safely no-op
    }

    // Try to receive data from one or more slabs
    for(int i = first; i < first + receive_ahead; i++){
        int iw = CP->WrapSlab(i);
        // try_receive() will run twice for every exchange: once to launch the receive,
        // and again to install the data
        if(left_exchanger[iw] != NULL) ret |= left_exchanger[iw]->try_receive();
        if(right_exchanger[iw] != NULL) ret |= right_exchanger[iw]->try_receive();
    }

    // Do a lightweight release of the send buffer and MPI handle
    for(int i = 0; i < P.cpd; i++){
        if(left_exchanger[i] != NULL){
            ret |= left_exchanger[i]->try_free_sendbuf();
        }
        if(right_exchanger[i] != NULL){
            ret |= right_exchanger[i]->try_free_sendbuf();
        }
    }

    return ret;
}

// Checked in the Finish precondition
int IsNeighborReceiveDone(int slab){
    if(neighbor_exchange_is_noop){
        return 1;  // safely no-op
    }

    int sw = CP->WrapSlab(slab);

    assertf(left_exchanger[sw] != NULL,
        "Checking neighbor receive status on slab {:d} not set up to receive?\n", sw);
    assertf(right_exchanger[sw] != NULL,
        "Checking neighbor receive status on slab {:d} not set up to receive?\n", sw);

    return left_exchanger[sw]->done_receive() && right_exchanger[sw]->done_receive();
}


// Called in the Epilogue
void TeardownNeighborExchange(){
    if(neighbor_exchange_is_noop){
        return;  // safely no-op
    }

    STDLOG(3,"Tearing down Neighbor Exchange\n");
    
    for(int i = 0; i < P.cpd; i++){
        if(left_exchanger[i] != NULL){
            assertf(left_exchanger[i]->done_send(), "left exchanger for slab {:d} not done send?\n", i);
            assertf(left_exchanger[i]->done_receive(), "left exchanger for slab {:d} not done receive?\n", i);
            delete left_exchanger[i];
        }
        if(right_exchanger[i] != NULL){
            assertf(right_exchanger[i]->done_send(), "right exchanger for slab {:d} not done send?\n", i);
            assertf(right_exchanger[i]->done_receive(), "right exchanger for slab {:d} not done receive?\n", i);
            delete right_exchanger[i];
        }
    }
    delete[] left_exchanger;
    delete[] right_exchanger;
}

#else // PARALLEL

void SetupNeighborExchange(int first, int nslab) { }
int AttemptNeighborReceive(int first, int receive_ahead){ return 0; }
void DoNeighborSend(int slab){ }
int IsNeighborReceiveDone(int slab){ return 1; }
void TeardownNeighborExchange(){ }

#endif
