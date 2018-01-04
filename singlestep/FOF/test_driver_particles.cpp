// This is the code to input or create particles for a simple global test.

class SortStruct {
  // We need this do sort particles
  public:
      auxstruct aux;
      posstruct pos;
      velstruct vel;
      accstruct acc;
      bool operator< (const SortStruct& c) const { return (aux.val<c.aux.val); }
};

class Particles {
    public:
    int cpd;
    int cpd3;
    int cpdhalf;
    int np;
    float invcpd;

    posstruct *pos;
    velstruct *vel;
    auxstruct *aux;
    accstruct *acc;
    int *cell;   // Store the starting point of each cell
    int *n;      // Store the size of each cell

    Particles() {
    	pos = NULL; vel = NULL; aux = NULL; acc = NULL; cell = NULL; n = NULL;
    }
    ~Particles() {
	if (pos!=NULL) free(pos); pos=NULL;
	if (vel!=NULL) free(vel); vel=NULL;
	if (aux!=NULL) free(aux); aux=NULL;
	if (acc!=NULL) free(acc); acc=NULL;
	if (cell!=NULL) free(cell); cell=NULL;
	if (n!=NULL) free(n); n=NULL;
    }

    void MakeParticles(int _cpd, int _np) {
        cpd = _cpd;
	invcpd = 1.0/cpd;
	cpdhalf = cpd/2;
	np = _np;
	cpd3 = cpd*cpd*cpd;
	pos = (posstruct *)malloc(sizeof(posstruct)*np);
	vel = (velstruct *)malloc(sizeof(velstruct)*np);
	aux = (auxstruct *)malloc(sizeof(auxstruct)*np);
	acc = (accstruct *)malloc(sizeof(accstruct)*np);
	cell = (int *)malloc(sizeof(int)*cpd3);
	n = (int *)malloc(sizeof(int)*cpd3);
	// Now divide up the particles into cells
	// This makes a uniform density, but with a lot more in the last cell
	int each = np/(cpd3+10);
	for (int c=0; c<cpd3; c++) {
	    cell[c] = each*c;
	    n[c] = each;
	}
	n[cpd3-1] = np-each*(cpd3-1);

	int total = 0;
	for (int c=0; c<cpd3; c++) total+=n[c];
	printf("Made %d particles, often %d per cell, %d in the last\n", total, each, np-each*(cpd3-1));

	srand48(123);
	for (int p=0; p<np; p++) aux[p].val = p;
	for (int p=0; p<np; p++) vel[p] = velstruct(0.0,0.0,0.0);
	for (int p=0; p<np; p++) {
	    // Make some cell-centered positions
	    pos[p].x = (drand48()-0.5)*invcpd;
	    pos[p].y = (drand48()-0.5)*invcpd;
	    pos[p].z = (drand48()-0.5)*invcpd;
	}
	return;
    }

    void WrapAndIndex() {
	for (int k=0; k<np; k++) {
	    // We're going to use vel to hold the cell ijk for now
	    vel[k].x = floor(pos[k].x*cpd);
	    vel[k].y = floor(pos[k].y*cpd);
	    vel[k].z = floor(pos[k].z*cpd);
	    pos[k] -= (vel[k]+float3(0.5,0.5,0.5))*invcpd;   // Cell-centered pos
	    // printf("%f %f %f\n", pos[k].x*cpd, pos[k].y*cpd, pos[k].z*cpd);
	    // assert(pos[k].x>=-invcpd/2.0 && pos[k].x<=invcpd/2.0);
	    // assert(pos[k].y>=-invcpd/2.0 && pos[k].y<=invcpd/2.0);
	    // assert(pos[k].z>=-invcpd/2.0 && pos[k].z<=invcpd/2.0);
	    // Now wrap and form the sorting index in aux.
	    while (vel[k].x>=cpd) vel[k].x-=cpd;
	    while (vel[k].y>=cpd) vel[k].y-=cpd;
	    while (vel[k].z>=cpd) vel[k].z-=cpd;
	    while (vel[k].x<0) vel[k].x+=cpd;
	    while (vel[k].y<0) vel[k].y+=cpd;
	    while (vel[k].z<0) vel[k].z+=cpd;
	    aux[k].val = ((int)vel[k].x*cpd+(int)vel[k].y)*cpd+(int)vel[k].z;
	    assert(aux[k].val<cpd3&&aux[k].val>=0);
	    aux[k].val <<= 20;
	    aux[k].val += k;
	    acc[k].x = 0.0;
	    acc[k].y = 0.0;
	    acc[k].z = 0.0;
	}
	// Now we have our particles.  Need to sort them by aux.
	SortStruct *s = new SortStruct[np];
	for (int j=0;j<np;j++) {
	    s[j].aux = aux[j];
	    s[j].acc = acc[j];
	    s[j].pos = pos[j];
	    s[j].vel = vel[j];
	}
	std::sort(s,s+np);
	for (int j=0;j<np;j++) {
	    aux[j] = s[j].aux;
	    acc[j] = s[j].acc;
	    pos[j] = s[j].pos;
	    vel[j] = s[j].vel;
	}
	delete[] s;
	// Now divide up the particles into cells
	cell = (int *)malloc(sizeof(int)*cpd3);
	n = (int *)malloc(sizeof(int)*cpd3);
	int p=0;
	for (int c=0;c<cpd3;c++) {
	    cell[c] = p;
	    while (p<np && (aux[p].val>>20)==c) p++;
	    n[c] = p-cell[c];
	}
	for (int j=0;j<np;j++) aux[j].val &= 0xfffff;  // Remove the cell info
	return;
    }

    void ReadParticles(int _cpd, const char fname[], posstruct offset) {
        // Read a bunch of particle positions from a file.  Grid them up for use.
	// Particles are assumed to be in the unit cell; will be wrapped.
	// This will be an ASCII file; first 3 numbers on each line are used.
        cpd = _cpd;
	invcpd = 1.0/cpd;
	np = 1e6;	// Just pick this as a maximum
	cpd3 = cpd*cpd*cpd;
	pos = (posstruct *)malloc(sizeof(posstruct)*np);
	vel = (velstruct *)malloc(sizeof(velstruct)*np);
	aux = (auxstruct *)malloc(sizeof(auxstruct)*np);
	acc = (accstruct *)malloc(sizeof(accstruct)*np);
	FILE *fp = fopen(fname,"r");
	char line[200];
	int k=0;
	float tmp[3];
	while (fgets(line, 200, fp)!=NULL) {
	    assert(k<np);
	    if (line[0]=='#') continue;
	    int ret = sscanf(line, "%f %f %f", tmp, tmp+1, tmp+2);
	    if (ret!=3) continue;
	    pos[k] = posstruct(tmp[0], tmp[1], tmp[2]);
	    pos[k] += offset;
	    k++;
	}
	np = k;
	fclose(fp);
	WrapAndIndex();
	return;
    }

    Cell GetCell(int i, int j, int k) {
	int id = (i*cpd+j)*cpd+k;
	Cell c;
	c.pos = pos+cell[id];
	c.vel = vel+cell[id];
	c.aux = aux+cell[id];
	c.acc = acc+cell[id];
	c.n = n[id];
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

Particles *PP;

