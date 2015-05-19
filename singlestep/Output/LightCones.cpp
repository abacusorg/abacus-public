/*
 * LightCones.cpp
 *
 *  Created on: Oct 1, 2012
 *      Author: dferrer
 *
 *
 *      Identify particles in the lightcone for a given set of origins and output them.
 *
 *      The particles are currently output in the following scheme:
 *      filename --"LCN/StepI/slabJ.lc"
 *      a
 *      da
 *      NCells
 *      cellheader (center, NP)
 *      *particles (pos, vel, z)*
 *      cellheader (center, NP)
 *      *particles (pos, vel, z)*
 *
 */


#define LCTOLERANCE 0.0

FLOAT3 *LCOrigin;

void getLightConeFN(int i,int slab, char * fn){ //ensure the necessary directories exist and return a filename for this light cone.
		char dir1[1080];
		char dir2[1080];
		CheckDirectoryExists(P.LightConeDirectory);
		sprintf(dir1,"%s/LC%.2d",P.LightConeDirectory,i);
		if(!FileExists(dir1)){ //TODO:This assumes the location is not a file. We may want to improve this
			mkdir(dir1,0777);
		}
		sprintf(dir2,"%s/LC%.2d/Step%.4d",P.LightConeDirectory,i,WriteState.FullStepNumber);
		if (!FileExists(dir2)){
			mkdir(dir2,0777);
		}
		CheckDirectoryExists(dir2);
		sprintf(fn,"%s/LC%.2d/Step%.4d/slab%.5d.lc",P.LightConeDirectory,i,WriteState.FullStepNumber,slab);

}

#define etaktoMpc 4000.
#define c_kms 299792.0
inline int inLightCone(posstruct pos, velstruct vel, int lcn, double tol1,double tol2){ //check if this position is in the current lightcone.
	double r1 = (cosm->today.etaK - cosm->current.etaK)*etaktoMpc/ReadState.BoxSizeMpc;
	double r2 = (cosm->today.etaK - cosm->next.etaK)*etaktoMpc/ReadState.BoxSizeMpc;
	if (r1 ==r2) return false;
	double r = (pos-LCOrigin[lcn]).norm();
	double vronc  = 0;
	if(r >1e-12) vronc = vel.dot(pos-LCOrigin[lcn])/r * ReadState.VelZSpace_to_kms /c_kms ;
	r = (r-r1) *1/(1+vronc) +r1; //apply velocity correction
	int result = ((r < r1+tol1) && (r > r2 - tol2));
	return result;

}

inline void interpolateParticle(posstruct & pos, velstruct &vel,accstruct acc, int lcn){
	double r1 = (cosm->today.etaK - cosm->current.etaK)*etaktoMpc/ReadState.BoxSizeMpc;
	double r2 = (cosm->today.etaK - cosm->next.etaK)*etaktoMpc/ReadState.BoxSizeMpc;

	double r = (pos-LCOrigin[lcn]).norm();
    double vronc  = 0;
    if(r >1e-12) vronc = vel.dot(pos-LCOrigin[lcn])/r * ReadState.VelZSpace_to_kms /c_kms;
	double deltar= (r-r1)/(r2-r1) *1/(1+vronc); //apply velocity correction
	if(deltar >=0 or ReadState.DeltaScaleFactor ==0){ //check that we are doing valid interpolation
		double delta_etaK = (WriteState.FirstHalfEtaKick+WriteState.LastHalfEtaKick)*deltar; //how much we need to kick the particle
		double delta_etaD = WriteState.DeltaEtaDrift*deltar; //how much to drift the particle

		pos = pos + delta_etaD * vel;
		vel = vel + delta_etaK * acc;
	}
	else{ //interpolate using the previous timestep information
		double r0 = (cosm->today.etaK -cosm->current.etaK +cosm->KickFactor(ReadState.ScaleFactor -ReadState.DeltaScaleFactor,ReadState.DeltaScaleFactor))*etaktoMpc/ReadState.BoxSizeMpc;
		double deltar2= (r-r1)/(r0-r1) *1/(1-vronc);

		//assert(deltar2 <=2); //we should not be going back more than one timestep}
		double delta_etaK = (ReadState.FirstHalfEtaKick+ReadState.LastHalfEtaKick)*deltar2; //how much we need to kick the particle
		double delta_etaD = ReadState.DeltaEtaDrift*deltar2; //how much to drift the particle

		pos = pos - delta_etaD * vel;
		vel = vel - delta_etaK * acc;
	}
}

void makeLightCone(int slab,int lcn){ //lcn = Light Cone Number
	//double r1 = cosm->today.H*(cosm->today.etaK - cosm->current.etaK)*4000/P.BoxSize;
        //double r2 = cosm->today.H*(cosm->today.etaK - cosm->next.etaK)*4000/P.BoxSize;
	//printf("r1: %f r2: %f \n",r1,r2);
	AppendArena *AA = new OutputLightcone();

	SlabType lightcone = (SlabType)((int)(LightCone0 + lcn));
	int headersize = 1024*1024;
	long long int slabtotal = 0;

	LBW->AllocateSpecificSize(lightcone, slab,
			Slab->size(slab)*(AA->sizeof_particle())
			+ PP->cpd*(PP->cpd)*(AA->sizeof_cell()) + headersize);
	AA->initialize(lightcone, slab, PP->cpd, ReadState.VelZSpace_to_Canonical);

	if (!P.OmitOutputHeader) {
	        AA->addheader((const char *) P.header());
	        AA->addheader((const char *) ReadState.header());
	        AA->finalize_header();
	}
	auxstruct as;
	uint64 mask = as.lightconemask(lcn);

	double kickfactor = WriteState.FirstHalfEtaKick;
	velstruct vel;
	posstruct pos;
	integer3 ijk(slab,0,0);
	for (ijk.y=0; ijk.y<PP->cpd; ijk.y++){
		for (ijk.z=0;ijk.z<PP->cpd;ijk.z++) {
			Cell c = PP->GetCell(ijk);
			accstruct *acc = PP->AccCell(ijk);
			long long int np;
			if(inLightCone(PP->CellCenter(ijk),0*pos,lcn,3.0/P.cpd, 3.0/P.cpd)){
				AA->addcell(ijk,1);
				np = 0;
				posstruct cc = PP->CellCenter(ijk);
#ifdef GLOBALPOS
				cc= 0*cc;
#endif
				for (int p=0;p<c.count();p++) {
					if(!c.aux[p].lightconedone(mask)){
						vel = (c.vel[p] - acc[p]*kickfactor);
						if(inLightCone((c.pos[p]+cc),vel,lcn,1.0,0)){
							pos = c.pos[p]+cc;
							interpolateParticle(pos,vel,acc[p],lcn);
							AA->addparticle(pos, vel, c.aux[p]);
							c.aux[p].setlightconedone(mask);
							np++;
						}
					}
				}
				AA->endcell();
			}
			slabtotal+=np;
			}

		}
	    // Write out this filename
	    char filename[1024];
	    getLightConeFN(lcn,slab,filename);
	    if(slabtotal) STDLOG(1,"Slab %d had %d particles in lightcone %lld\n",slab,slabtotal,lcn);
	    LBW->ResizeSlab(lightcone, slab, AA->bytes_written());
	    LBW->WriteArena(lightcone, slab, IO_DELETE, IO_NONBLOCKING, filename);
	    delete AA;
}

