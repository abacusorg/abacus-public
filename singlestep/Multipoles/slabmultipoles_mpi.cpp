
class SlabMultipolesMPI : public SlabMultipoles { 
public:
    
    SlabMultipolesMPI(int order, int cpd); 
    ~SlabMultipolesMPI(void);

    void ComputeMultipoleFFT( int x,  FLOAT3 *spos, int *count, int *offset,  FLOAT3 *cc, MTCOMPLEX *tmp);
    void ComputeFFTZ(int x, MTCOMPLEX *outslab);
    int CheckAnyMPIDone();
    int IsMPIDone(int slab);

private:
    int node_ky_size;  // size of ky on this node, after transpose with z

    double *reducedtmp;
    Complex *ztmp;
    
    Complex **sendbuf;
    Complex **recvbuf;

    MPI_Request *handle;
    int *mpi_status;  // 0 = waiting; 1 = send/recv; 2 = done

    int *sendcounts;  // num elements to send to rank i
    int *senddispls;  // offset of chunk to send to rank i

    int *recvcounts;  // num elements to send to rank i
    int *recvdispls;  // offset of chunk to send to rank i

    void FFTY(Complex *out, const double *in);
    void FFTZ(MTCOMPLEX *out, const Complex *in);
    
    void MakeSendRecvBufs(Complex **sbuf, Complex **rbuf, const Complex *in);
    void DoMPIAllToAll(int slab, MPI_Request *handle, const Complex *sbuf, Complex *rbuf);

    const MPI_Datatype mpi_dtype = MPI_C_DOUBLE_COMPLEX;
};

SlabMultipolesMPI::~SlabMultipolesMPI(void) {
    free(reducedtmp);
    free(transposetmp);
    free(ztmp);

    delete[] sendbuf;
    delete[] recvbuf;
    delete[] handle;
    delete[] sendcounts;
    delete[] senddispls;
    delete[] recvcounts;
    delete[] recvdispls;
}

SlabMultipolesMPI::SlabMultipolesMPI(int order, int cpd)
    : SlabMultipoles(order, cpd)  {

    assert(MPI_size_z > 1);

    node_ky_size = node_cpdp1half;

    sendbuf = new Complex*[cpd]();  // sendbuf[x] will be allocated on-demand
    recvbuf = new Complex*[cpd]();

    handle = new MPI_Request[cpd];
    mpi_status = new int[cpd];
    for(int i = 0; i < cpd; i++){
        handle[i] = MPI_REQUEST_NULL;
        mpi_status[i] = 0;
    }

    sendcounts = new int[MPI_size_z];
    senddispls = new int[MPI_size_z];
    for(int i = 0; i < MPI_size_z; i++){
        sendcounts[i] = all_node_ky_size[i]*node_z_size*rml;
        senddispls[i] = (i>0) ? (senddispls[i-1] + sendcounts[i-1]) : 0;
    }

    recvcounts = new int[MPI_size_z];
    recvdispls = new int[MPI_size_z];
    for(int i = 0; i < MPI_size_z; i++){
        recvcounts[i] = node_ky_size*all_node_z_size[i]*rml;
        recvdispls[i] = (i>0) ? (recvdispls[i-1] + recvcounts[i-1]) : 0;
    }

    assert(posix_memalign((void **) &reducedtmp, PAGE_SIZE, sizeof(double) * cpd * node_z_size * rml) == 0);
    // TODO: could halve MPI data if we use MTCOMPLEX
    assert(posix_memalign((void **) &transposetmp, PAGE_SIZE, sizeof(Complex) * cpdp1half * rml * node_z_size) == 0);
    
    assert(posix_memalign((void **) &ztmp, PAGE_SIZE,
        sizeof(Complex) * node_ky_size * cpd * rml) == 0);  // all z, some ky
}

void SlabMultipolesMPI::FFTY(Complex *out, const double *in) {
    // in:  [cpd, node_z_size, rml]
    // out: [node_z_size, rml, (cpd+1)/2]

    FFTMultipole.Start();
    #pragma omp parallel for schedule(static)
    for(int64_t z = 0; z < node_z_size; z++){
        int64_t g = omp_get_thread_num();
        for(int64_t m=0;m<rml;m++) {
            for(int64_t y=0;y<cpd;y++)
                in_r2c[g][y] = in[y*node_z_size*rml + z*rml + m];
            fftw_execute(plan_forward_r2c_1d[g]);

            for(int64_t y=0; y < cpdp1half; y++)
                out[ z*cpdp1half*rml + m*cpdp1half + y]  = out_r2c[g][y];
        }
    }
    FFTMultipole.Stop();
}

void SlabMultipolesMPI::MakeSendRecvBufs(Complex **sbuf, Complex **rbuf, const Complex *in){
    // in: [node_z_size, rml, (cpd+1)/2]
    // out: [(cpd+1)/2, rml, node_z_size]

    FillMPIBufs.Start();

    // TODO: could make these arenas, would get good reuse
    assert(posix_memalign((void **) sbuf, PAGE_SIZE,
        sizeof(Complex) * cpdp1half * node_z_size * rml) == 0);  // all ky, some z
    assert(posix_memalign((void **) rbuf, PAGE_SIZE,
        sizeof(Complex) * node_ky_size * cpd * rml) == 0);  // all z, some ky

    // fill sendbuf
    #pragma omp parallel for schedule(static)
    for(int64_t y = 0; y < cpdp1half; y++){
        for(int64_t m = 0; m < rml; m++){
            for(int64_t z = 0; z < node_z_size; z++){
                (*sbuf)[y*node_z_size*rml + m*node_z_size + z] = in[z*rml*cpdp1half + m*cpdp1half + y];
            }
        }
    }

    FillMPIBufs.Stop();
}

void SlabMultipolesMPI::DoMPIAllToAll(int slab, MPI_Request *handle, const Complex *sbuf, Complex *rbuf){
    // FUTURE: MPI-4 supports MPI_Ialltoallv_c, with 64-bit counts.
    // But all HPC MPI-3 implementations seem to support > 2 GB data, as long as the counts are 32-bit.
    AllToAll.Start();

    MPI_Ialltoallv((const void *) sbuf, sendcounts,
        senddispls, mpi_dtype,
        (void *) rbuf, recvcounts,
        recvdispls, mpi_dtype, comm_multipoles_z,
        handle);
    
    mpi_status[slab] = 1;

    AllToAll.Stop();
}

int SlabMultipolesMPI::CheckAnyMPIDone(){
    CheckMPI.Start();
    int ret = 0;
    for(int i = 0; i < cpd; i++){
        if(mpi_status[i] == 1){
            int done = 0;
            MPI_Test(&handle[i], &done, MPI_STATUS_IGNORE);
            if(done){
                STDLOG(2, "Multipoles y-z MPI transpose done on slab %d\n", i);
                mpi_status[i] = 2;
                free(sendbuf[i]);
                ret = 1;
            }
        }
    }
    CheckMPI.Stop();
    return ret;
}

int SlabMultipolesMPI::IsMPIDone(int slab){
    return mpi_status[CP->WrapSlab(slab)] == 2;
}

void SlabMultipolesMPI::ComputeMultipoleFFT( int x, FLOAT3 *spos, 
                     int *count, int *offset, FLOAT3 *cc, MTCOMPLEX *_out_unused) {
    STimer wc;
    PTimer _kernel, _c2r;
    pdouble localMassSlabX[nprocs];
    padded<double3> localdipole[nprocs];
    double *localMassSlabZ = new double[nprocs*node_z_size];
    
    // Init with OpenMP to preserve thread locality
    #pragma omp parallel for schedule(static)
    for(int64_t g = 0; g < nprocs; g++){
        localMassSlabX[g] = 0.;
        localdipole[g] = double3(0.);
        for(int64_t z = 0; z < node_z_size; z++)
            localMassSlabZ[g*node_z_size + z] = 0.;
    }
    
    wc.Start();
    // compute the cell-by-cell multipoles for this node's subslab
    NUMA_FOR(y,0,cpd)
        int64_t g = omp_get_thread_num();
        for(int64_t z = 0; z < node_z_size; z++) {
            _kernel.Start();
            int64_t i = y*node_z_size + z;
            EvaluateCartesianMultipoles( &(spos[offset[i]]),
                                    count[i], cc[i], cartesian[g] );
            _kernel.Stop();
            
            _c2r.Start();
            double *reducedcell = &reducedtmp[y*node_z_size*rml + z*rml];
            DispatchCartesian2Reduced(order, cartesian[g], reducedcell);
            
            double Mxyz = reducedcell[0];
            localMassSlabX[g] += Mxyz;
            localMassSlabZ[g*node_z_size + z] += Mxyz;
            MassSlabY[y] += Mxyz;  // thread dimension, no race TODO: still false sharing
            
            localdipole[g] += double3(reducedcell[rmap(1,0,0) ],
                                reducedcell[rmap(0,1,0) ],
                                reducedcell[rmap(0,0,1) ] );
            _c2r.Stop();
        }
    }
    wc.Stop();
    
    // do thread reductions
    for(int64_t g = 0; g < nprocs; g++){
        MassSlabX[x] += localMassSlabX[g];
        globaldipole += localdipole[g];
        for(int64_t z = 0; z < node_z_size; z++){
            // MassSlabZ uses "global" z indices
            MassSlabZ[z + node_z_start] += localMassSlabZ[g*node_z_size + z];
        }
    }
    delete[] localMassSlabZ;


    FFTY(transposetmp, reducedtmp);
    
    MakeSendRecvBufs(&sendbuf[x], &recvbuf[x], transposetmp);
    DoMPIAllToAll(x, &handle[x], sendbuf[x], recvbuf[x]);

    double seq = _kernel.Elapsed() + _c2r.Elapsed();
    
    struct timespec seq_kernel = scale_timer(_kernel.Elapsed()/seq, wc.get_timer() );
    struct timespec seq_c2r = scale_timer(_c2r.Elapsed()/seq, wc.get_timer() );

    MultipoleKernel.increment(seq_kernel);
    MultipoleC2R.increment(seq_c2r);
}


void SlabMultipolesMPI::ComputeFFTZ(int x, MTCOMPLEX *outslab){
    // out: [cpd, rml, node_ky_size]

    assertf(mpi_status[x] == 2, "ComputeFFTZ() called before MPI receive done on slab %d?\n", x);

    // unpack recvbuf into ztmp
    // recvbuf holds each node's chunk, one after the other

    UnpackRecvBuf.Start();

    Complex *rbuf = recvbuf[x];

    // recvbuf: [MPI_size_z, node_ky_size, rml, all_node_z_size[node]]
    // ztmp: [node_ky_size, rml, cpd]
    #pragma omp parallel for schedule(static)
    for(int64_t y = 0; y < node_ky_size; y++){
        for(int64_t m = 0; m < rml; m++){
            for(int64_t node = 0; node < MPI_size_z; node++){
                int64_t zstart = all_node_z_start[node];
                int64_t zsize = all_node_z_size[node];
                for(int64_t zoff = 0; zoff < zsize; zoff++){
                    ztmp[y*rml*cpd + m*cpd + (zstart + zoff)] =
                        rbuf[zstart*node_ky_size*rml + y*rml*zsize + m*zsize + zoff];
                }
            }
        }
    }

    UnpackRecvBuf.Stop();

    // TODO: decide if we're going to do MPI in complex<double> or not
    MTCOMPLEX *rbuf32 = (MTCOMPLEX *) recvbuf[x];

    // Now FFT from ztmp back into recvbuf[x]
    FFTZ(rbuf32, ztmp);

    FFTZTranspose.Start();
    
    // and finally transpose from recvbuf[x] into outslab
    // recvbuf: [node_ky_size, rml, cpd]
    // outslab: [cpd, rml, node_ky_size]
    
    #pragma omp parallel for schedule(static)
    for(int64_t kz = 0; kz < cpd; kz++){
        for(int64_t m = 0; m < rml; m++){
            for(int64_t ky = 0; ky < node_ky_size; ky++){
                outslab[kz*rml*node_ky_size + m*node_ky_size + ky] =
                    rbuf32[ky*rml*cpd + m*cpd + kz];
            }
        }
    }

    free(rbuf);

    FFTZTranspose.Stop();
}


void SlabMultipolesMPI::FFTZ(MTCOMPLEX *out, const Complex *in) {
    // out: shape [node_ky_size, rml, cpd]
    // in: shape [node_ky_size, rml, cpd]

    FFTZMultipole.Start();

    // collapse(2): thread over the combined y*m dimension
    #pragma omp parallel for schedule(static) collapse(2)
    for(int64_t y = 0; y < node_ky_size; y++){
        for(int64_t m = 0; m < rml; m++) {
            int g = omp_get_thread_num();
            for(int64_t z = 0; z < cpd; z++)
                in_1d[g][z] = in[y*rml*cpd + m*cpd + z];
            
            fftw_execute( plan_forward_c2c_1d[g] );

            for(int64_t kz = 0; kz < cpd; kz++)
                out[y*rml*cpd + m*cpd + kz] = out_1d[g][kz];
        }
    }

    FFTZMultipole.Stop();
}
