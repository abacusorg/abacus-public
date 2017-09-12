#include "appendarena.cpp"

AppendArena *get_AA_by_format(const char* format){
    AppendArena *AA;

    if (strcmp(format,"RVdouble")==0) {
        STDLOG(1,"Using Output Format RVdouble\n");
        AA = new OutputRVdouble();
    } else if (strcmp(format,"Packed")==0) {
        STDLOG(1,"Using Output Format Packed\n");
        AA = new OutputPacked();
    } else if (strcmp(format,"Heitmann")==0) {
        STDLOG(1,"Using Output Format Heitmann\n");
        AA = new OutputHeitmann();
    } else if (strcmp(format,"RVdoubleTag")==0) {
        STDLOG(1,"Using Output Format RVdoubleTag\n");
        AA = new OutputRVdoubleTag();
    }
    else {
        QUIT("Unrecognized case: OutputFormat = %s\n", format);
    }
    
    return AA;
}

uint64 Output_TimeSlice(int slab) {
    AppendArena *AA;
    FLOAT vscale;

    AA = get_AA_by_format(P.OutputFormat);

    // Setup the Arena
    int headersize = 1024*1024;
    LBW->AllocateSpecificSize(TimeSlice, slab, 
    	   Slab->size(slab)*(AA->sizeof_particle())
	+ PP->cpd*(PP->cpd)*(AA->sizeof_cell()) + headersize);
    AA->initialize(TimeSlice, slab, PP->cpd, ReadState.VelZSpace_to_Canonical);

    // Write the header to its own file
    if(slab == 0){
        char filename[1024];
        sprintf(filename, "%s/slice%5.3f/header",  
            P.OutputDirectory, 
            ReadState.Redshift);
        std::ofstream headerfile;
        headerfile.open(filename);
        headerfile << P.header();
        headerfile << ReadState.header();
        headerfile << "\nOutputType = \"TimeSlice\"\n";
        headerfile.close();
    }
        
    // and also add the header to the slab file
    if (!P.OmitOutputHeader) {
        AA->addheader((const char *) P.header());
        AA->addheader((const char *) ReadState.header());
        char head[1024];
        sprintf(head, "\nOutputType = \"TimeSlice\"\n"); 
        AA->addheader((const char *) head);
        sprintf(head, "SlabNumber = %d\n", slab);
        AA->addheader((const char *) head);
        // For sanity, be careful that the previous lines end with a \n!
        AA->finalize_header();
    }

    // Now scan through the cells
    FLOAT kickfactor = WriteState.FirstHalfEtaKick;	   // Amount to unkick.
    velstruct vel;
    integer3 ijk(slab,0,0);
    uint64 n_added = 0;
    for (ijk.y=0; ijk.y<PP->cpd; ijk.y++) 
        for (ijk.z=0;ijk.z<PP->cpd;ijk.z++) {
            Cell c = PP->GetCell(ijk);
            accstruct *acc = PP->AccCell(ijk);

            // We sometimes use the maximum velocity to scale.
            // But we do not yet have the global velocity (slab max will be set in Finish,
            // while the global max has to wait for all slabs to be done).
            // What is available after the kick is the max_component_velocity in each cell.
            vscale = c.ci->max_component_velocity/ReadState.VelZSpace_to_Canonical;	
            // The maximum velocity of this cell, converted to ZSpace unit-box units.

            // Start the cell
            AA->addcell(ijk, vscale);
            // Now pack the particles
            for (int p=0;p<c.count();p++) {
                vel = (c.vel[p] - acc[p]*kickfactor);    // We supply in code units
                AA->addparticle(c.pos[p], vel, c.aux[p]);
            }
            n_added += c.count();
            AA->endcell();
        }

    LBW->ResizeSlab(TimeSlice, slab, AA->bytes_written());

    // Write out this filename
    char filename[1024]; 	
    sprintf(filename, "%s/slice%5.3f/%s.z%5.3f.slab%04d.dat",  
    	P.OutputDirectory, 
	ReadState.Redshift,
	P.SimName,
	ReadState.Redshift,
	slab);
    LBW->WriteArena(TimeSlice, slab, IO_DELETE, IO_NONBLOCKING, filename);
    delete AA;
    
    return n_added;
}
