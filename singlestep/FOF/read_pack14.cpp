// This is a particle header for pack14 inputs


#include "cell_header.h"
#include "pack14_storage.cpp"
#include <sys/stat.h>

class Particles {
public:
    int cpd;
    int cpd3;
    int np;
    float invcpd; 

    // These will all be arrays [0,cpd] that point to the buffers for each slab.
    posstruct **posslab;
    velstruct **velslab;
    auxstruct **auxslab;
    int **cellslab;
    int **nslab;

    Particles(int _cpd) {
	cpd = _cpd;
	invcpd = 1.0/cpd;
	cpd3 = cpd*cpd*cpd;
	np = 0;

	posslab = new posstruct *[cpd];
	velslab = new velstruct *[cpd];
	auxslab = new auxstruct *[cpd];
	cellslab = new int *[cpd];
	nslab    = new int *[cpd];
	for (int j=0; j<cpd; j++) {
	    posslab[j] = NULL;
	    velslab[j] = NULL;
	    auxslab[j] = NULL;
	    cellslab[j] = NULL;
	    nslab[j] = NULL;
	}
	return;
    }

    void free_slab(int slab) {
	if (posslab[slab] !=NULL) free(posslab[slab]);  posslab[slab] = NULL;
	if (velslab[slab] !=NULL) free(velslab[slab]);  velslab[slab] = NULL;
	if (auxslab[slab] !=NULL) free(auxslab[slab]);  auxslab[slab] = NULL;
	if (cellslab[slab]!=NULL) free(cellslab[slab]); cellslab[slab] = NULL;
	if (nslab[slab]!=NULL) free(nslab[slab]); nslab[slab] = NULL;
	return;
    }

    ~Particles() {
	for (int j=0; j<cpd; j++) free_slab(j);
	delete[] posslab;
	delete[] velslab;
	delete[] auxslab;
	delete[] cellslab;
	delete[] nslab;
	return;
    }

    int read_slab_pack14(char fname[], int slab) {
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

	int *cell = cellslab[slab] = (int *) malloc(sizeof(int)*cpd*cpd);
	int *n = nslab[slab] = (int *) malloc(sizeof(int)*cpd*cpd);
	// Initialize to empty cells
	for (int j=0; j<cpd*cpd; j++) cell[j] = 0;
	for (int j=0; j<cpd*cpd; j++) n[j] = 0;

	int nump = 0;
	int thiscell = -1;

	// Skip the file head
	int c, clast;
	if ((clast=getc(fp))==EOF) return 0;
	while ((c=getc(fp))!=EOF) {
	    if (clast=='^B'&&c=='\n') break;
	    clast = c;
	}

	while (fread(&particle, sizeof(pack14), 1, fp)!=EOF) {
	    if (particle.iscell()) {
		// Starting a new cell, so finish the old one.
		if(thiscell>=0) n[thiscell] = nump-cell[thiscell];
		cellhead = particle.unpack_cell();
		integer3 ijk = cellhead.cellid();
		thiscell = ijk.y*cpd+ijk.z;    // The ID of this cell.
		cell[thiscell] = nump;	// The starting point
	    } else {
		assert(cellhead.islegal());  
		particle.unpack(posd,veld,id,cellhead);
		pos[nump] = posd+double3(0.5,0.5,0.5)-(double3(0.5,0.5,0.5)+cellhead.cellid())*invcpd;
		    // Move back to cell-centered positions
		vel[nump] = veld;
		aux[nump].val = id;
		nump++;
	    }
	}
	fclose(fp);
        np += nump;
	return np;
    }

    Cell GetCell(int i, int j, int k) {
        int id = j*cpd+k;
        Cell c;
        c.pos = posslab[i]+cellslab[i][id];
        c.vel = velslab[i]+cellslab[i][id];
        c.aux = auxslab[i]+cellslab[i][id];
        c.n = nslab[i][id];
        return c;
    }
    Cell GetCell(integer3 ijk) {
        return GetCell(ijk.x, ijk.y, ijk.z);
    }

};

Particles *PP;
