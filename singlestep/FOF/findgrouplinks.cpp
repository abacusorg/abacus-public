/* This is the code that searches for links between CellGroups in neighboring
cells.  When a link is found, it is placed in the GroupLinkList GLL.
Aside from those GroupLinks, there is no other output from this code.

The code creates a great deal of temporary objects to speed this
search.  Each cell has a face of width linking_length that contains
the PseudoParticles that are on the face.  A PseudoParticle can
either be a normal particle, with its index in the cell, or it can
be a FaceGroup, which is the set of particles from a CellGroup that
are on the face as well as the radius of a bounding sphere for those
particles.

The code creates the 5 relevant faces for all cells in the slab, plus
one face for the neighboring slab.  Then it searches all faces in order
to find all of the links that involve the given slab.

We make no effort to restrict down to edges or corners; we expect that
the faces are already small enough that the bookkeeping overheads 
would overwhelm the value.
*/

class FaceGroup {
    public:
    int start;
    int n; 
        // FaceGroupParticles are in locations [start, start+n)
    int cellgroupID;
    FaceGroup(int _start, int _n, int _id) {
        start = _start;
	n = _n;
	cellgroupID = _id;
	return;
    }
};

typedef FOFparticle FacePseudoParticle;
typedef FOFparticle FaceGroupParticle;

class CellFaceSlab {
    public:
    SlabAccum<FacePseudoParticle> pseudoParticles;
	// We create a set of searchable FacePseudoParticles for each face.
	// index = cell group number for singlets, -FaceGroup number for multiplets
    SlabAccum<float> pseudoRadii;
	// Every pseudo particle has a bounding radius; 0 for singlets.

    SlabAccum<FaceGroupParticle> faceParticles;
	// The actual face particles in the multiplets, index undefined
    SlabAccum<FaceGroup> faceGroups;
        // The lookup information for each multiplet, plus its cell group number.

    int edgebit;	// The edge we're going to select
    int slab; 		// The slab we're computing (wrapped)
    int slab_prewrap; 	// The slab we're told to compute (unwrapped)

    CellFaceSlab(int _slab, int _edgebit, int cpd, int maxsize) {
	pseudoParticles.setup(cpd, maxsize); // ETC
	pseudoRadii.setup(cpd, maxsize);
	faceParticles.setup(cpd, maxsize);
	faceGroups.setup(cpd, maxsize);
	edgebit = _edgebit;
	slab_prewrap = _slab;
	slab = GFC->WrapSlab(_slab);
    }

    void CreateSlabFace();

    ~CellFaceSlab() {
	pseudoParticles.destroy(); 
	pseudoRadii.destroy();
	faceParticles.destroy();
	faceGroups.destroy();
	return;
    }
};

class FaceSet {
    public:
    // This is a simple accumulation, so that we can make other code more concise
    PencilAccum<FacePseudoParticle> *pP;
    PencilAccum<float>              *pR;
    PencilAccum<FaceGroupParticle>  *fP;
    PencilAccum<FaceGroup>          *fG;
    // We have to guarantee these float4 to be SSE aligned.  
    // No cache lines either
    /// float *radius; 
    /// FacePseudoParticle *midpos; 

    CellFaceSlab *face;
    FaceSet(CellFaceSlab *_face) { 
    	face = _face; 
	int ret;
	/// ret = posix_memalign((void **)&radius, 64, 4*sizeof(float)); assert(ret==0);
	/// ret = posix_memalign((void **)&midpos, 64, sizeof(FacePseudoParticle)); assert(ret==0);
	pP = NULL;
	pR = NULL;
	fP = NULL;
	fG = NULL;
    }
    ~FaceSet() {
	/// free(radius);
	/// free(midpos);
    }

    void StartPencil(int j) {
	pP = face->pseudoParticles.StartPencil(j);
	pR = face->pseudoRadii.StartPencil(j);
	fP = face->faceParticles.StartPencil(j);
	fG = face->faceGroups.StartPencil(j);
    }
    void FinishCell() {
	pP->FinishCell();
	pR->FinishCell();
	fP->FinishCell();
	fG->FinishCell();
    }
    void FinishPencil() {
	pP->FinishPencil();
	pR->FinishPencil();
	fP->FinishPencil();
	fG->FinishPencil();
    }

    void ProcessCell(int j, int k) {
	// Search this cell for face particles
	Cell c = PP->GetCell(face->slab, j, k);
	integer3 wc = GFC->WrapCell(face->slab, j, k);
	CellPtr<CellGroup> cg = GFC->cellgroups[wc.x][wc.y][wc.z];
	const FOFloat boundary = GFC->boundary;   // The boundaries, cell-centered

	// Store a SSE version of 0.5 for later
	float half1 = 0.5;
	__m128 half = _mm_load_ps1(&half1);
	int facegroupnumstart = fG->get_pencil_size();
		// We're going to start the facegroup numbering from the start of the cell
		// So store the pencil index here and subtract it later

	int fPcellstart = fP->get_pencil_size();

	// printf("Cell %d %d %d for bit %02x\n", face->slab, j, k, face->edgebit >> 25 );
	
	// Take a look at the CellGroups, which include singlet boundary particles
	for (int g=0; g<cg.size(); g++) {
	    // printf("CG %d: %d %d %02x\n", g, cg[g].start, cg[g].size(), cg[g].n>>25);
	    if (cg[g].test(face->edgebit)!=0) {
		// Group is in the face
		// But now we want to find the particles in the face
		int fPstart = fP->get_pencil_size()-fPcellstart;
		int fGsize = 0;
		FOFparticle BBmin, BBmax;

		// Look at the particles in this group
		posstruct *ppos = c.pos+cg[g].start;
		for (int p=0; p<cg[g].size(); p++, ppos++)
		    if (test_particle(*ppos,face->edgebit,boundary)) {
			// Particle is in the face; save it
			// This is done regardless of whether this will
			// end up unused, as a single particle pseudoGroup.
			FaceGroupParticle pp = FaceGroupParticle(*ppos,0);
			// printf("FGPart: %f %f %f -> %f %f %f %f\n",
				// ppos->x, ppos->y, ppos->z,
				// pp.x/FOF_RESCALE, pp.y/FOF_RESCALE, pp.z/FOF_RESCALE, pp.n);
			fP->append(pp);
			/// __m128 newp = _mm_loadu_ps((float *)&pp);
			// Note that the 4th element here is garbage
			if (fGsize>0) {
			    BBmin.min(pp);
			    BBmax.max(pp);
			    /// BBmin = _mm_min_ps(BBmin, newp);
			    /// BBmax = _mm_max_ps(BBmax, newp);
			} else {
			    BBmin = pp;
			    BBmax = pp;
			    /// BBmin = newp;
			    /// BBmax = newp;
			}
			fGsize++;
		    }
		// Done looking at the particles

		if (fGsize==1) {
		    // We only found one particle, so it just goes on 
		    // the PseudoParticle list with its cellgroup number.
		    // No need to create a FaceGroup
		    /// _mm_store_ps((float *)midpos, BBmin);
		    /// midpos->n = g;
		    BBmin.n = g;
		    // printf("PsPart: %f %f %f, n=%f r=%f\n",
				// midpos->x/FOF_RESCALE, midpos->y/FOF_RESCALE, midpos->z/FOF_RESCALE, midpos->n, 0.0);
		    pP->append(BBmin);
		    pR->append(0.0);
		} else if (fGsize>1) {
		    // We found multiple particles, so need to compute
		    // the pseudoparticle.  

		    // Add this to the FaceGroup list
		    int facegroupnum = fG->get_pencil_size() - facegroupnumstart;
		    	// This numbering is from the start of the cell
		    fG->append(FaceGroup(fPstart, fGsize, g));
		        // fPstart is relative to the start of the cell
		    assert(fP->get_pencil_size()==fPcellstart+fPstart+fGsize);

		    // Use the Min/Max to create the pseudo particle
		    FOFloat radius = BBmax.diff2(&BBmin)*0.25;;
		    BBmin.x = 0.5*(BBmin.x+BBmax.x);
		    BBmin.y = 0.5*(BBmin.y+BBmax.y);
		    BBmin.z = 0.5*(BBmin.z+BBmax.z);
		    BBmin.n = 0.5*(BBmin.n+BBmax.n);
		    /// __m128 mid = _mm_mul_ps(_mm_add_ps(BBmin, BBmax),half);
		    /// _mm_store_ps((float *)(midpos), mid);
		    /// __m128 diff = _mm_sub_ps(BBmax, BBmin);
		    /// diff = _mm_mul_ps(diff,diff);  // Square diameter
		    /// _mm_store_ps(radius, diff);
		    /// radius[0] += radius[1]+radius[2];
		    /// radius[0] *= 0.25;   // Now this is the square radius

		    BBmin.n = -1-facegroupnum;  // Overwrite the index
		    // printf("PsPart: %f %f %f, n=%d r=%f\n",
				// midpos->x/FOF_RESCALE, midpos->y/FOF_RESCALE, midpos->z/FOF_RESCALE, midpos->index(), sqrt(radius[0])/FOF_RESCALE);
		    pP->append(BBmin);
		    pR->append(sqrt(radius));
		} else {
		    assert(fGsize>0);
		}
	    }
	}
	FinishCell();
	// Done with this cell
    }

};  // End class FaceSet


void CreateFaces( CellFaceSlab &xm, CellFaceSlab &xp,
	CellFaceSlab &ym, CellFaceSlab &yp, CellFaceSlab &zm, CellFaceSlab &zp) {
    // We're going to create all the faces, looping over the j,k cells for the slab.
    // This routine does all 5 at once for each cell, on the hopes that
    // we only load the particles into cache once.

    #pragma omp parallel for schedule(static)
    for (int j=0; j<GFC->cpd; j++) {
	// The FaceSet holds the PencilAccum for this face.
	FaceSet fxm(&xm);
	FaceSet fxp(&xp);
	FaceSet fym(&ym);
	FaceSet fyp(&yp);
	FaceSet fzm(&zm);
	FaceSet fzp(&zp);

	fxm.StartPencil(j);
	fxp.StartPencil(j);
	fym.StartPencil(j);
	fyp.StartPencil(j);
	fzm.StartPencil(j);
	fzp.StartPencil(j);

	// xp is the one oddball, because it is in the neighboring slab.
	// So we do all of that pencil first, so that it benefits from 
	// whatever memory look-ahead might be there.
	for (int k=0; k<GFC->cpd; k++) {
	    fxp.ProcessCell(j,k);
	}
	// The other 5 proceed together, so that we re-use the particles
	// that have been loaded into cache.
	for (int k=0; k<GFC->cpd; k++) {
	    fxm.ProcessCell(j,k);
	    fym.ProcessCell(j,k);
	    fyp.ProcessCell(j,k);
	    fzm.ProcessCell(j,k);
	    fzp.ProcessCell(j,k);
	}
	
	// Done with this pencil
	fxm.FinishPencil();
	fxp.FinishPencil();
	fym.FinishPencil();
	fyp.FinishPencil();
	fzm.FinishPencil();
	fzp.FinishPencil();
    }
    return;
}

// ===================== Faces: Searching for Links =================================

void SearchPair(CellFaceSlab &c1, int j1, int k1, 
	    CellFaceSlab &c2, int j2, int k2) {
    // Given two CellFaceSlabs and the desired cells, search for Links.
    // This starts by comparing the pseudoParticles.
    // If they are close enough, then consider whether they are compound.
    // If either or both is compound, search all those pairs.
    // Whenever one finds a link, mark the pair and move on to the next one.

    int i1 = c1.slab;
    int i2 = c2.slab;

    // We also need to account for the cell-centered positions
    FOFparticle *offset;
    int ret = posix_memalign((void **)&offset, 64, sizeof(FOFparticle));   assert(ret==0); // Need this aligned
    int del_i = c2.slab_prewrap - c1.slab_prewrap;   // Need this info before wrapping
    offset->x = GFC->invcpd*FOF_RESCALE*(del_i);
    offset->y = GFC->invcpd*FOF_RESCALE*(j2-j1);
    offset->z = GFC->invcpd*FOF_RESCALE*(k2-k1);
    offset->n = 0.0;
    // This is what we should add to the 2nd positions

    // Now wrap them
    j1 = GFC->WrapSlab(j1);
    k1 = GFC->WrapSlab(k1);
    j2 = GFC->WrapSlab(j2);
    k2 = GFC->WrapSlab(k2);
    // printf("%d %d %d to %d %d %d\n", i1, j1, k1, i2, j2, k2);


    CellPtr<FacePseudoParticle> pP1 = c1.pseudoParticles[j1][k1];
    CellPtr<FacePseudoParticle> pP2 = c2.pseudoParticles[j2][k2];
    CellPtr<float> pR1 = c1.pseudoRadii[j1][k1];
    CellPtr<float> pR2 = c2.pseudoRadii[j2][k2];
    CellPtr<FaceGroup> fG1 = c1.faceGroups[j1][k1];
    CellPtr<FaceGroup> fG2 = c2.faceGroups[j2][k2];
    CellPtr<FaceGroupParticle> fP1 = c1.faceParticles[j1][k1];
    CellPtr<FaceGroupParticle> fP2 = c2.faceParticles[j2][k2];

    FOFloat link = GFC->linking_length*FOF_RESCALE;
    FOFloat bsq = link*link;

    for (int p=0; p<pP1.size(); p++) {
	FOFloat b1 = link+pR1[p];
        for (int q=0; q<pP2.size(); q++) {
	    // We need the square distance compared to (link+r1+r2)**2
	    FOFloat b2 = b1+pR2[q];
	    // printf("%d %f to %d %f = %f vs %f\n",
	    	// p, pR1[p]/FOF_RESCALE, q, pR2[q]/FOF_RESCALE, 
		// pP1[p].diff2(pP2.ptr(q),offset)/FOF_RESCALE/FOF_RESCALE, b2*b2/FOF_RESCALE/FOF_RESCALE);
	    if (pP1[p].diff2(pP2.ptr(q),offset)<b2*b2) {
	        // We've found a pair of pseudogroups
		// if (pR1[p]==0) 
		// printf("Candidate: %d %d %d %d to %d %d %d %d\n",
			// i1, j1, k1, pP1[p].index(),
			// i2, j2, k2, pP2[q].index());
		// Any time index()<0, we get valgrind errors
		if (pP1[p].index()>=0) {
		    if (pP2[q].index()>=0) {
			// A pair of singlets.  Just link them.
			// printf("Link: %d %d %d %d to %d %d %d %d\n",
				// i1, j1, k1, pP1[p].index(),
			        // i2, j2, k2, pP2[q].index());
			GFC->GLL->DoublePush(LinkID(i1, j1, k1, pP1[p].index()),
					   LinkID(i2, j2, k2, pP2[q].index()));
		    } else {
			// Group 2 was a multiplet.  Need to search it.
			FaceGroup *fg = fG2.ptr(-1-pP2[q].index());
			FaceGroupParticle *fp = fP2.ptr(fg->start);
			for (int qq=0; qq<fg->n; qq++, fp++) {
			    FOFloat d2 = pP1[p].diff2(fp,offset);
			    if (d2<bsq) {
				// We've found a pair
				// printf("Link: %d %d %d %d to %d %d %d %d %d\n",
				    // i1, j1, k1, pP1[p].index(),
				    // i2, j2, k2, pP2[q].index(), fg->cellgroupID);
				GFC->GLL->DoublePush(LinkID(i1, j1, k1, pP1[p].index()),
						   LinkID(i2, j2, k2, fg->cellgroupID));
				break;  // No need to continue
			    }
			}
		    }
		} else {
		    if (pP2[q].index()>=0) {
			// Group 1 was a multiplet.  Need to search it.
			FaceGroup *fg = fG1.ptr(-1-pP1[p].index());
			FaceGroupParticle *fp = fP1.ptr(fg->start);
			for (int pp=0; pp<fg->n; pp++, fp++) {
			    FOFparticle *p2 = pP2.ptr(q);
			    FOFloat d2 = fp->diff2(p2,offset);
			    if (d2<bsq) {
			    // if (pP2[q].diff2(fp)<bsq) 
				// We've found a pair
				// printf("Link: %d %d %d %d %d to %d %d %d %d\n",
				    // i1, j1, k1, pP1[p].index(),fg->cellgroupID,
				    // i2, j2, k2, pP2[q].index()); 
				GFC->GLL->DoublePush(LinkID(i1, j1, k1, fg->cellgroupID),
						   LinkID(i2, j2, k2, pP2[q].index()));
				break;  // No need to continue
			    }
			}
		    } else {
			// Both groups are multiplets
			FaceGroup *fg1 = fG1.ptr(-1-pP1[p].index());
			FaceGroup *fg2 = fG2.ptr(-1-pP2[q].index());
			FaceGroupParticle *fp1 = fP1.ptr(fg1->start);
			for (int pp=0; pp<fg1->n; pp++, fp1++) {
			    FaceGroupParticle *fp2 = fP2.ptr(fg2->start);
			    for (int qq=0; qq<fg2->n; qq++, fp2++) {
				if (fp1->diff2(fp2,offset)<bsq) {
				    // We've found a pair
				    // printf("Link: %d %d %d %d %d to %d %d %d %d %d\n",
					// i1, j1, k1, pP1[p].index(),fg1->cellgroupID,
					// i2, j2, k2, pP2[q].index(), fg2->cellgroupID); 
				    GFC->GLL->DoublePush(LinkID(i1, j1, k1, fg1->cellgroupID),
						       LinkID(i2, j2, k2, fg2->cellgroupID));
				    goto endsearch;
				}
			    }
			}
		    }
		}
	        endsearch: ;
	    }  // Done with this pseudoParticle pair.
	}
    }
    free(offset);
    return; 
}

void FindGroupLinks(int slab) {
    // Find the CellGroup borders for all cells in slab, including
    // those to slab-1.
    // Then search them, adding to the List of Links.
    GFC->CreateFaceTime.Start();
    slab = GFC->WrapSlab(slab);
    uint64 Ltot_start = GFC->GLL->length;

    // TODO: This maxsize is likely inefficient
    CellFaceSlab xp(slab-1, XP_BIT, GFC->cpd, 1024);
    CellFaceSlab xm(slab,   XM_BIT, GFC->cpd, 1024);
    CellFaceSlab yp(slab,   YP_BIT, GFC->cpd, 1024);
    CellFaceSlab ym(slab,   YM_BIT, GFC->cpd, 1024);
    CellFaceSlab zp(slab,   ZP_BIT, GFC->cpd, 1024);
    CellFaceSlab zm(slab,   ZM_BIT, GFC->cpd, 1024);

    // Now load up these slabs
    // We do all at once, in order to get repeated access to individual cells.
    CreateFaces(xm, xp, ym, yp, zm, zp);

    STDLOG(2,"XP for slab %d: %u pseudoParticles, %u faceParticles, %u faceGroups\n",
	slab, xp.pseudoParticles.get_slab_size(), xp.faceParticles.get_slab_size(), xp.faceGroups.get_slab_size());
    STDLOG(2,"XM for slab %d: %u pseudoParticles, %u faceParticles, %u faceGroups\n",
	slab, xm.pseudoParticles.get_slab_size(), xm.faceParticles.get_slab_size(), xm.faceGroups.get_slab_size());
    STDLOG(2,"YP for slab %d: %u pseudoParticles, %u faceParticles, %u faceGroups\n",
	slab, yp.pseudoParticles.get_slab_size(), yp.faceParticles.get_slab_size(), yp.faceGroups.get_slab_size());
    STDLOG(2,"YM for slab %d: %u pseudoParticles, %u faceParticles, %u faceGroups\n",
	slab, ym.pseudoParticles.get_slab_size(), ym.faceParticles.get_slab_size(), ym.faceGroups.get_slab_size());
    STDLOG(2,"ZP for slab %d: %u pseudoParticles, %u faceParticles, %u faceGroups\n",
	slab, zp.pseudoParticles.get_slab_size(), zp.faceParticles.get_slab_size(), zp.faceGroups.get_slab_size());
    STDLOG(2,"ZM for slab %d: %u pseudoParticles, %u faceParticles, %u faceGroups\n",
	slab, zm.pseudoParticles.get_slab_size(), zm.faceParticles.get_slab_size(), zm.faceGroups.get_slab_size());

    GFC->pPtot += 
    	xp.pseudoParticles.get_slab_size() + xm.pseudoParticles.get_slab_size() +
    	yp.pseudoParticles.get_slab_size() + ym.pseudoParticles.get_slab_size() +
    	zp.pseudoParticles.get_slab_size() + zm.pseudoParticles.get_slab_size();

    GFC->fPtot += 
    	xp.faceParticles.get_slab_size() + xm.faceParticles.get_slab_size() +
    	yp.faceParticles.get_slab_size() + ym.faceParticles.get_slab_size() +
    	zp.faceParticles.get_slab_size() + zm.faceParticles.get_slab_size();

    GFC->fGtot += 
    	xp.faceGroups.get_slab_size() + xm.faceGroups.get_slab_size() +
    	yp.faceGroups.get_slab_size() + ym.faceGroups.get_slab_size() +
    	zp.faceGroups.get_slab_size() + zm.faceGroups.get_slab_size();

    /* 
    // These files are consistent, even with different threading.
    char fname[200];
    sprintf(fname, "Z/pP.xp.%02d", slab); xp.pseudoParticles.dump_to_file(fname);
    sprintf(fname, "Z/pP.xm.%02d", slab); xm.pseudoParticles.dump_to_file(fname);
    sprintf(fname, "Z/pP.yp.%02d", slab); yp.pseudoParticles.dump_to_file(fname);
    sprintf(fname, "Z/pP.ym.%02d", slab); ym.pseudoParticles.dump_to_file(fname);
    sprintf(fname, "Z/pP.zp.%02d", slab); zp.pseudoParticles.dump_to_file(fname);
    sprintf(fname, "Z/pP.zm.%02d", slab); zm.pseudoParticles.dump_to_file(fname);

    sprintf(fname, "Z/fP.xp.%02d", slab); xp.faceParticles.dump_to_file(fname);
    sprintf(fname, "Z/fP.xm.%02d", slab); xm.faceParticles.dump_to_file(fname);
    sprintf(fname, "Z/fP.yp.%02d", slab); yp.faceParticles.dump_to_file(fname);
    sprintf(fname, "Z/fP.ym.%02d", slab); ym.faceParticles.dump_to_file(fname);
    sprintf(fname, "Z/fP.zp.%02d", slab); zp.faceParticles.dump_to_file(fname);
    sprintf(fname, "Z/fP.zm.%02d", slab); zm.faceParticles.dump_to_file(fname);

    sprintf(fname, "Z/fG.xp.%02d", slab); xp.faceGroups.dump_to_file(fname);
    sprintf(fname, "Z/fG.xm.%02d", slab); xm.faceGroups.dump_to_file(fname);
    sprintf(fname, "Z/fG.yp.%02d", slab); yp.faceGroups.dump_to_file(fname);
    sprintf(fname, "Z/fG.ym.%02d", slab); ym.faceGroups.dump_to_file(fname);
    sprintf(fname, "Z/fG.zp.%02d", slab); zp.faceGroups.dump_to_file(fname);
    sprintf(fname, "Z/fG.zm.%02d", slab); zm.faceGroups.dump_to_file(fname);

    sprintf(fname, "Z/pR.xp.%02d", slab); xp.pseudoRadii.dump_to_file(fname);
    sprintf(fname, "Z/pR.xm.%02d", slab); xm.pseudoRadii.dump_to_file(fname);
    sprintf(fname, "Z/pR.yp.%02d", slab); yp.pseudoRadii.dump_to_file(fname);
    sprintf(fname, "Z/pR.ym.%02d", slab); ym.pseudoRadii.dump_to_file(fname);
    sprintf(fname, "Z/pR.zp.%02d", slab); zp.pseudoRadii.dump_to_file(fname);
    sprintf(fname, "Z/pR.zm.%02d", slab); zm.pseudoRadii.dump_to_file(fname);
    */

    /* 
    // These also matched between threadings:
    char fname[200];
    sprintf(fname, "Z/pP.xp.%02d", slab); xp.pseudoParticles.dump_cells_to_file(fname);
    sprintf(fname, "Z/pP.xm.%02d", slab); xm.pseudoParticles.dump_cells_to_file(fname);
    sprintf(fname, "Z/pP.yp.%02d", slab); yp.pseudoParticles.dump_cells_to_file(fname);
    sprintf(fname, "Z/pP.ym.%02d", slab); ym.pseudoParticles.dump_cells_to_file(fname);
    sprintf(fname, "Z/pP.zp.%02d", slab); zp.pseudoParticles.dump_cells_to_file(fname);
    sprintf(fname, "Z/pP.zm.%02d", slab); zm.pseudoParticles.dump_cells_to_file(fname);

    sprintf(fname, "Z/fP.xp.%02d", slab); xp.faceParticles.dump_cells_to_file(fname);
    sprintf(fname, "Z/fP.xm.%02d", slab); xm.faceParticles.dump_cells_to_file(fname);
    sprintf(fname, "Z/fP.yp.%02d", slab); yp.faceParticles.dump_cells_to_file(fname);
    sprintf(fname, "Z/fP.ym.%02d", slab); ym.faceParticles.dump_cells_to_file(fname);
    sprintf(fname, "Z/fP.zp.%02d", slab); zp.faceParticles.dump_cells_to_file(fname);
    sprintf(fname, "Z/fP.zm.%02d", slab); zm.faceParticles.dump_cells_to_file(fname);

    sprintf(fname, "Z/fG.xp.%02d", slab); xp.faceGroups.dump_cells_to_file(fname);
    sprintf(fname, "Z/fG.xm.%02d", slab); xm.faceGroups.dump_cells_to_file(fname);
    sprintf(fname, "Z/fG.yp.%02d", slab); yp.faceGroups.dump_cells_to_file(fname);
    sprintf(fname, "Z/fG.ym.%02d", slab); ym.faceGroups.dump_cells_to_file(fname);
    sprintf(fname, "Z/fG.zp.%02d", slab); zp.faceGroups.dump_cells_to_file(fname);
    sprintf(fname, "Z/fG.zm.%02d", slab); zm.faceGroups.dump_cells_to_file(fname);

    sprintf(fname, "Z/pR.xp.%02d", slab); xp.pseudoRadii.dump_cells_to_file(fname);
    sprintf(fname, "Z/pR.xm.%02d", slab); xm.pseudoRadii.dump_cells_to_file(fname);
    sprintf(fname, "Z/pR.yp.%02d", slab); yp.pseudoRadii.dump_cells_to_file(fname);
    sprintf(fname, "Z/pR.ym.%02d", slab); ym.pseudoRadii.dump_cells_to_file(fname);
    sprintf(fname, "Z/pR.zp.%02d", slab); zp.pseudoRadii.dump_cells_to_file(fname);
    sprintf(fname, "Z/pR.zm.%02d", slab); zm.pseudoRadii.dump_cells_to_file(fname);
    */

    GFC->CreateFaceTime.Stop();
    GFC->FindLinkTime.Start();
    
    // Search each cell, putting the GroupLinks onto that insertlist.
    #pragma omp parallel for schedule(dynamic,1)
    for (int j=0; j<GFC->cpd; j++) {
        for (int k=0; k<GFC->cpd; k++) {
	    SearchPair(xp, j  , k-1, xm, j, k);
	    SearchPair(xp, j  , k  , xm, j, k);
	    SearchPair(xp, j  , k+1, xm, j, k);
	    SearchPair(xp, j+1, k-1, xm, j, k);
	    SearchPair(xp, j+1, k  , xm, j, k);
	    SearchPair(xp, j+1, k+1, xm, j, k);
	    SearchPair(xp, j-1, k-1, xm, j, k);
	    SearchPair(xp, j-1, k  , xm, j, k);
	    SearchPair(xp, j-1, k+1, xm, j, k);
	    SearchPair(ym, j+1, k-1, yp, j, k);
	    SearchPair(ym, j+1, k  , yp, j, k);
	    SearchPair(ym, j+1, k+1, yp, j, k);
	    SearchPair(zm, j, k+1  , zp, j, k);
	}
    }
    GFC->GLL->CollectGaps();   // We've now found all of the links from this slab
    GFC->Ltot += GFC->GLL->length - Ltot_start;

    GFC->FindLinkTime.Stop();
    return;
}
