// Copyright 2012-2025 The Abacus Developers
// SPDX-License-Identifier: GPL-3.0-or-later

/*

SlabTaylorMPI is an implementation of the virtual SlabTaylor class
that does the Taylors in the 2D code.  Compared to SlabTaylorLocal,
SlabTaylorMPI does the Y-FFT first, and then executes an MPI transpose
to gather Z and execute the Z-FFT.

*/

class SlabTaylorMPI : public SlabTaylor {
public:
    SlabTaylorMPI(int order, int cpd);
    ~SlabTaylorMPI(void);

    void EvaluateSlabTaylor(int x, FLOAT3 *FA, const FLOAT3 *spos,
                        const int *count, const int *offset, const int *ghost_offsets,
                        const FLOAT3 *cc, const MTCOMPLEX *TaylorCoefficients);
    
    void ComputeIFFTZAndMPI(int x, MTCOMPLEX *outslab);
    
    int CheckAnyMPIDone();
    int IsMPIDone(int slab);

private:
    int node_ky_size;  // size of ky on this node, after transpose with z
    mpi_complex_t *transposetmp;
    mpi_complex_t *ytmp;

    mpi_complex_t **sendbuf;
    mpi_complex_t **recvbuf;

    MPI_Request *handle;
    int *mpi_status;  // 0 = waiting; 1 = send/recv; 2 = done

    int *sendcounts;  // num elements to send to rank i
    int *senddispls;  // offset of chunk to send to rank i

    int *recvcounts;  // num elements to send to rank i
    int *recvdispls;  // offset of chunk to send to rank i


    void InverseFFTZ(mpi_complex_t *out, const MTCOMPLEX *in);
    void InverseFFTY(double *out, const mpi_complex_t *in);

    void MakeSendRecvBufs(mpi_complex_t **sbuf, mpi_complex_t **rbuf, const mpi_complex_t *in);
    void DoMPIAllToAll(int slab, MPI_Request *handle, const mpi_complex_t *sbuf, mpi_complex_t *rbuf);

};


SlabTaylorMPI::~SlabTaylorMPI(){
    free(transposetmp);
    free(ytmp);

    delete[] sendbuf;
    delete[] recvbuf;
    delete[] handle;
    delete[] sendcounts;
    delete[] senddispls;
    delete[] recvcounts;
    delete[] recvdispls;
}


SlabTaylorMPI::SlabTaylorMPI(int order, int cpd)
    : SlabTaylor(order, cpd) {

    assert(MPI_size_z > 1);
    
    node_ky_size = node_cpdp1half;

    assert(posix_memalign((void **) &transposetmp, PAGE_SIZE, sizeof(mpi_complex_t)*node_ky_size*rml*cpd) == 0);
    assert(posix_memalign((void **) &ytmp, PAGE_SIZE, sizeof(mpi_complex_t)*node_z_size*rml*cpdp1half) == 0);

    sendbuf = new mpi_complex_t*[cpd]();  // sendbuf[x] will be allocated on-demand
    recvbuf = new mpi_complex_t*[cpd]();

    handle = new MPI_Request[cpd];
    mpi_status = new int[cpd];
    for(int i = 0; i < cpd; i++){
        handle[i] = MPI_REQUEST_NULL;
        mpi_status[i] = 0;
    }

    sendcounts = new int[MPI_size_z];
    senddispls = new int[MPI_size_z];
    for(int i = 0; i < MPI_size_z; i++){
        sendcounts[i] = all_node_z_size[i]*node_ky_size*rml;
        senddispls[i] = (i>0) ? (senddispls[i-1] + sendcounts[i-1]) : 0;
    }

    recvcounts = new int[MPI_size_z];
    recvdispls = new int[MPI_size_z];
    for(int i = 0; i < MPI_size_z; i++){
        recvcounts[i] = all_node_ky_size[i]*node_z_size*rml;
        recvdispls[i] = (i>0) ? (recvdispls[i-1] + recvcounts[i-1]) : 0;
    }
}


void SlabTaylorMPI::ComputeIFFTZAndMPI(int x, MTCOMPLEX *tslab){
    // tslab: [cpd, rml, node_ky_size]

    // Just finished the ParallelConvolve MPI, so we have all kz for some ky.
    // Do the z-iFFT.
    InverseFFTZ(transposetmp, tslab);

    MakeSendRecvBufs(&sendbuf[x], &recvbuf[x], transposetmp);
    DoMPIAllToAll(x, &handle[x], sendbuf[x], recvbuf[x]);
}


void SlabTaylorMPI::DoMPIAllToAll(int slab, MPI_Request *handle, const mpi_complex_t *sbuf, mpi_complex_t *rbuf){
    AllToAll.Start();
    MPI_Ialltoallv((const void *) sbuf, sendcounts,
        senddispls, mpi_dtype,
        (void *) rbuf, recvcounts,
        recvdispls, mpi_dtype, comm_taylors_z,
        handle);
    
    mpi_status[slab] = 1;
    AllToAll.Stop();
}


int SlabTaylorMPI::CheckAnyMPIDone(){
    CheckMPI.Start();
    int ret = 0;
    for(int i = 0; i < cpd; i++){
        if(mpi_status[i] == 1){
            int done = 0;
            MPI_Test(&handle[i], &done, MPI_STATUS_IGNORE);
            if(done){
                STDLOG(2, "Taylors y-z MPI transpose done on slab {:d}\n", i);
                mpi_status[i] = 2;
                free(sendbuf[i]);
                ret = 1;
            }
        }
    }
    CheckMPI.Stop();
    return ret;
}


int SlabTaylorMPI::IsMPIDone(int slab){
    return mpi_status[CP->WrapSlab(slab)] == 2;
}


void SlabTaylorMPI::MakeSendRecvBufs(mpi_complex_t **sbuf, mpi_complex_t **rbuf, const mpi_complex_t *in){
    // in: [node_ky_size, rml, cpd]
    // sbuf: [cpd, rml, node_ky_size]

    FillMPIBufs.Start();
    assert(posix_memalign((void **) sbuf, PAGE_SIZE,
        sizeof(mpi_complex_t) * node_ky_size * cpd * rml) == 0);  // all z, some ky
    assert(posix_memalign((void **) rbuf, PAGE_SIZE,
        sizeof(Complex) * cpdp1half * node_z_size * rml) == 0);  // all ky, some z
        // N.B. rbuf is intentionally overallocated to type Complex so we can reuse it as tbuf

    #pragma omp parallel for schedule(static)
    for(int64_t z = 0; z < cpd; z++){
        for(int64_t m = 0; m < rml; m++){
            for(int64_t ky = 0; ky < node_ky_size; ky++){
                (*sbuf)[ z*rml*node_ky_size + m*node_ky_size + ky ] = 
                    in[ky*rml*cpd + m*cpd + z];
            }
        }
    }
    FillMPIBufs.Stop();
}


void SlabTaylorMPI::InverseFFTZ(mpi_complex_t *out, const MTCOMPLEX *in){
    // in: [cpd, rml, node_ky_size]
    // out: [node_ky_size, rml, cpd]
    
    FFTZTaylor.Start();
    #pragma omp parallel for schedule(static)
    for(int64_t y = 0; y < node_ky_size; y++){
        for(int64_t m = 0; m < rml; m++){
            int g = omp_get_thread_num();
            
            for(int64_t kz = 0; kz < cpd; kz++)
                in_1d[g][kz] = (Complex) in[kz*rml*node_ky_size + m*node_ky_size + y];
            
            fftw_execute( plan_backward_c2c_1d[g] );
            
            for(int64_t z = 0; z < cpd; z++)
                out[y*rml*cpd + m*cpd + z] = (mpi_complex_t) out_1d[g][z];
        }
    }
    FFTZTaylor.Stop();
}

void SlabTaylorMPI::InverseFFTY(double *out, const mpi_complex_t *in){
    // in: [node_z_size, rml, cpdp1half]
    // out: [node_z_size, rml, cpd]

    #pragma omp parallel for schedule(static) collapse(2)
    for(int64_t z = 0; z < node_z_size; z++){
        for(int64_t m = 0; m < rml; m++) {
            int g = omp_get_thread_num();
            for(int64_t ky = 0; ky < cpdp1half; ky++)
                in_c2r[g][ky] = static_cast<Complex>(in[z*rml*cpdp1half + m*cpdp1half + ky]);
            fftw_execute( plan_backward_c2r_1d[g] );
            for(int64_t y = 0; y < cpd; y++)
                out[z*rml*cpd + m*cpd + y] = out_c2r[g][y];
        }
    }
}


void SlabTaylorMPI::EvaluateSlabTaylor(int x, FLOAT3 *FA, const FLOAT3 *spos,
                                        const int *count, const int *offset, const int *ghost_offset,
                                        const FLOAT3 *cc, const MTCOMPLEX *_taylors [[maybe_unused]]){
    // FA: particle accelerations

    UnpackRecvBuf.Start();
    mpi_complex_t *rbuf = recvbuf[x];

    // unpack the MPI rbuf into ytmp
    // rbuf: [MPI_size_z, node_z_size, rml, all_node_ky_size[node]]
    // ytmp: [node_z_size, rml, cpdp1half]
    #pragma omp parallel for schedule(static) collapse(2)
    for(int64_t z = 0; z < node_z_size; z++){
        for(int64_t m = 0; m < rml; m++){
            for(int64_t node = 0; node < MPI_size_z; node++){
                int64_t kysize = all_node_ky_size[node];
                int64_t kystart = all_node_ky_start[node];
                for(int64_t ky = 0; ky < kysize; ky++){
                    ytmp[ z*rml*cpdp1half + m*cpdp1half + (kystart + ky) ] = 
                        rbuf[ kystart*node_z_size*rml + z*rml*kysize + m*kysize + ky ];
                }
            }
        }
    }
    UnpackRecvBuf.Stop();

    FFTTaylor.Start();
    double *tbuf = (double *) rbuf;  // reuse recvbuf, was overallocated
    InverseFFTY(tbuf, ytmp);
    FFTTaylor.Stop();

    STimer wc;
    PTimer _r2c, _tkernel;
    wc.Start();

    NUMA_FOR(y,0,cpd, NO_CLAUSE, FALLBACK_DYNAMIC){
        int gh_off = ghost_offset[y];
        int g = omp_get_thread_num();
        for(int64_t z = 0; z < node_z_size; z++) {
            int64_t i = y*node_z_size + z;
            //_r2c.Start();
            CellTaylorFromPencil(y, cartesian[g], &tbuf[z*rml*cpd]);
            //_r2c.Stop();
            
            //_tkernel.Start();
            FLOAT3 *aa = &(FA[offset[i]]);
            memset(aa, 0, sizeof(FLOAT3)*count[i]);

            EvaluateTaylor( cartesian[g], 
                               cc[i], count[i], (float3*) &spos[offset[i] + gh_off], aa);
            //_tkernel.Stop();

        }
    }
    NUMA_FOR_END;

    free(rbuf);
    
    wc.Stop();
    
    /*double seq = _r2c.Elapsed() + _tkernel.Elapsed();
    double f_r2c = _r2c.Elapsed()/seq;
    double f_kernel = _tkernel.Elapsed()/seq;

    struct timespec  seq_r2c = scale_timer(f_r2c, wc.get_timer() );
    struct timespec  seq_tkernel = scale_timer(f_kernel, wc.get_timer() );

    TaylorR2C.increment( seq_r2c  );
    TaylorKernel.increment( seq_tkernel );*/
    TaylorKernel.increment( wc.get_timer() );
}
