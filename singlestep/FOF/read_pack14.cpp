// This is a particle header for pack14 inputs


#include "cell_header.h"
#include "pack14_storage.cpp"
#include <sys/stat.h>

#include "particle_subsample.cpp"

// Unpack a FieldTimeSlice pack14 slab into pos,vel,aux,cellinfo slabs
// This is used by our standalone_fof pipeline
int64_t unpack_slab_pack14(int slab, double taggable_frac) {
    void *rawslab = SB->GetSlabPtr(FieldTimeSlice, slab);
    uint64 rawsize = SB->SlabSizeBytes(FieldTimeSlice, slab);
    FILE *buffer_file = fmemopen(rawslab, rawsize, "rb");
    assert(buffer_file != NULL);

    // Fast-forward past header
    uint64 header_size = HeaderStream::SkipHeaderFromFP(buffer_file);
    STDLOG(3,"Skipped a %d byte header\n", header_size);

    // This includes cell headers, so it's actually an upper limit
    // We will allocate SB slabs of this size, then shrink later
    // We will return the number of particles actually read.
    uint64 datasize = rawsize - header_size;
    assert(datasize % 14 == 0);
    uint64 maxnp = datasize/14;
    STDLOG(3,"Allocating space for %d particles\n", maxnp);

    posstruct *pos = (posstruct *) SB->AllocateSpecificSize(PosSlab,slab,maxnp*sizeof(posstruct), RAMDISK_NO);
    velstruct *vel = (velstruct *) SB->AllocateSpecificSize(VelSlab,slab,maxnp*sizeof(velstruct), RAMDISK_NO);
    auxstruct *aux = (auxstruct *) SB->AllocateSpecificSize(AuxSlab,slab,maxnp*sizeof(auxstruct), RAMDISK_NO);
    cellinfo *ci = (cellinfo *) SB->AllocateSpecificSize(CellInfoSlab,slab,maxnp*sizeof(cellinfo), RAMDISK_NO);
        // this one is grossly overallocated! but L0 might have more than cpd^2 cells...

    // Initialize to empty cells
    int cpd = CP->cpd;
    for (int j=0; j<maxnp; j++)
        ci[j].makenull();

    int64_t nump = 0;
    int thiscell = -1;

    double3 posd, veld;
    uint64_t id;
    cell_header cellhead;
    pack14 particle;
    while (fread(&particle, sizeof(pack14), 1, buffer_file)==1) {
        if (particle.iscell()) {
            // Starting a new cell, so finish the old one.
            if(thiscell>=0)
                ci[thiscell].count = ci[thiscell].active = nump - ci[thiscell].startindex;
            cellhead = particle.unpack_cell();
            integer3 ijk = cellhead.cellid();
            assert(ijk.maxcomponent() < cpd);
            thiscell = ijk.y*cpd + ijk.z;    // The ID of this cell.
            assert(thiscell < cpd*cpd);

            // Make sure we haven't seen this cell before.
            // This is a valid pack14 pattern, but in Abacus we expect all particles in be gathered in cell-order
            // TODO: to support L0 outputs, we need to push everything to the insert list (or index the L0 cells for fast lookup and merging)
            assert(ci[thiscell].startindex == 0);

            ci[thiscell].startindex = nump;  // The starting point
        } else {
            assert(cellhead.islegal());
            particle.unpack(posd,veld,id,cellhead);
            pos[nump] = posd - CP->CellCenter(cellhead.cellid());
                // Move back to cell-centered positions
            vel[nump] = veld*ReadState.VelZSpace_to_Canonical;
            aux[nump].aux = id;
            // Note: pack14 gives us the PID, but not the full aux field. So we have to restore any other aux fields.
            if (is_subsample_particle((int64_t) id, taggable_frac))
                aux[nump].set_taggable();
            nump++;
        }
    }
    // Finish the last cell
    if(thiscell>=0)
        ci[thiscell].count = ci[thiscell].active = nump - ci[thiscell].startindex;
    fclose(buffer_file);

    // Shrink to the actual read size
    SB->ResizeSlab(PosSlab, slab, nump*sizeof(posstruct));
    SB->ResizeSlab(VelSlab, slab, nump*sizeof(velstruct));
    SB->ResizeSlab(AuxSlab, slab, nump*sizeof(auxstruct));
    // TODO: figure out how much we can shrink CellInfoSlab

    return nump;
}


#ifdef STANDALONE_FOF
class CellParticles {
public:
    int cpd;
    int cpdhalf;
    int cpd3;
    int np;
    float invcpd; 

    // These will all be arrays [0,cpd] that point to the buffers for each slab.
    posstruct **posslab;
    velstruct **velslab;
    auxstruct **auxslab;
    accstruct **accslab;
    int **cellslab;
    int **nslab;

    CellParticles(int _cpd) {
        cpd = _cpd;
        invcpd = 1.0/cpd;
        cpdhalf = cpd/2;
        cpd3 = cpd*cpd*cpd;
        np = 0;

        posslab = new posstruct *[cpd];
        velslab = new velstruct *[cpd];
        auxslab = new auxstruct *[cpd];
        accslab = new accstruct *[cpd];
        cellslab = new int *[cpd];
        nslab    = new int *[cpd];
        for (int j=0; j<cpd; j++) {
            posslab[j] = NULL;
            velslab[j] = NULL;
            auxslab[j] = NULL;
            accslab[j] = NULL;
            cellslab[j] = NULL;
            nslab[j] = NULL;
        }
        return;
    }

    void free_slab(int slab) {
        if (posslab[slab] !=NULL) free(posslab[slab]);  posslab[slab] = NULL;
        if (velslab[slab] !=NULL) free(velslab[slab]);  velslab[slab] = NULL;
        if (auxslab[slab] !=NULL) free(auxslab[slab]);  auxslab[slab] = NULL;
        if (accslab[slab] !=NULL) free(accslab[slab]);  accslab[slab] = NULL;
        if (cellslab[slab]!=NULL) free(cellslab[slab]); cellslab[slab] = NULL;
        if (nslab[slab]!=NULL) free(nslab[slab]); nslab[slab] = NULL;
        return;
    }

    ~CellParticles() {
        for (int j=0; j<cpd; j++) free_slab(j);
        delete[] posslab;
        delete[] velslab;
        delete[] auxslab;
        delete[] accslab;
        delete[] cellslab;
        delete[] nslab;
        return;
    }

    int read_slab_pack14(const char fname[], int slab, double taggable_frac) {
        // We are going to allocate space for the particles and cells, so 
        // pass a pointer to this slab's pointers.
        // Return the number of particles actually read.
        double3 posd, veld;
        uint64_t id;
        cell_header cellhead;
        pack14 particle;

        FILE *fp = fopen(fname,"rb");
        assert(fp!=NULL);

        struct stat st;
        assert(stat(fname,&st)==0);
        int max = st.st_size/14;     // This includes the header, but we don't care if
            // we overallocate a bit.

        assert(posslab[slab]==NULL);   // We shouldn't overwrite something

        posstruct *pos = posslab[slab] = (posstruct *)malloc(sizeof(posstruct)*max);
        velstruct *vel = velslab[slab] = (velstruct *)malloc(sizeof(velstruct)*max);
        auxstruct *aux = auxslab[slab] = (auxstruct *)malloc(sizeof(auxstruct)*max);
        accstruct *acc = accslab[slab] = (accstruct *)malloc(sizeof(accstruct)*max);

        int *cell = cellslab[slab] = (int *) malloc(sizeof(int)*cpd*cpd);
        int *n = nslab[slab] = (int *) malloc(sizeof(int)*cpd*cpd);
        // Initialize to empty cells
        for (int j=0; j<cpd*cpd; j++) cell[j] = 0;
        for (int j=0; j<cpd*cpd; j++) n[j] = 0;

        int nump = 0;
        int thiscell = -1;

        // Skip the file head
        int c, clast;
        int numhead =0;
        if ((clast=getc(fp))==EOF) return 0;
        while ((c=getc(fp))!=EOF) {
            numhead++;
            if (clast==''&&c=='\n') break;
            clast = c;
        }
        // printf("Skipped a %d byte header\n", numhead);

        while (fread(&particle, sizeof(pack14), 1, fp)==1) {
            if (particle.iscell()) {
                // Starting a new cell, so finish the old one.
                if(thiscell>=0) n[thiscell] = nump-cell[thiscell];
                cellhead = particle.unpack_cell();
                integer3 ijk = cellhead.cellid();
                assert(ijk.maxcomponent() < cpd);
                thiscell = ijk.y*cpd+ijk.z;    // The ID of this cell.
                cell[thiscell] = nump;        // The starting point
            } else {
                assert(cellhead.islegal());  
                particle.unpack(posd,veld,id,cellhead);
                pos[nump] = posd+double3(0.5,0.5,0.5)-(double3(0.5,0.5,0.5)+cellhead.cellid())*invcpd;
                    // Move back to cell-centered positions
                vel[nump] = veld;
                aux[nump].aux = id;
                if (is_subsample_particle((int64_t) id, taggable_frac)) aux[nump].set_taggable();
                acc[nump] = accstruct(0.0);
                nump++;
            }
        }
        // Finish the last cell
        if(thiscell>=0) n[thiscell] = nump-cell[thiscell];
        fclose(fp);
        np += nump;
        return nump;
    }

    Cell GetCell(int i, int j, int k) {
        int id = j*cpd+k;
        Cell c;
        c.pos = posslab[i]+cellslab[i][id];
        c.vel = velslab[i]+cellslab[i][id];
        c.aux = auxslab[i]+cellslab[i][id];
        c.acc = accslab[i]+cellslab[i][id];
        cellinfo *ci = new cellinfo;
        ci->makenull();
        c.ci = ci;
        c.ci->count = nslab[i][id];
        return c;
    }
    Cell GetCell(integer3 ijk) {
        return GetCell(ijk.x, ijk.y, ijk.z);
    }
    float3 CellCenter(int i, int j, int k) {
        float3 v;
        v.x = (i-cpdhalf)*invcpd;;
        v.y = (j-cpdhalf)*invcpd;;
        v.z = (k-cpdhalf)*invcpd;;
        return v;
    }
};

Cell::~Cell(){
    delete cellinfo;
}

CellParticles *CP;
#endif
