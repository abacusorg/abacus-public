// Copyright 2012-2025 The Abacus Developers
// SPDX-License-Identifier: GPL-3.0-or-later

/*
 * This class creates an AppendArena object of the desired output type
 * and repeatedly calls "addparticle()" on it.  The actual implementation
 * of the writing functions is in appendarena.cpp.
 *
*/



#include "appendarena.cpp"

AppendArena *get_AA_by_format(const std::string &format){
    AppendArena *AA;

    if (format == "RVdouble") {
        STDLOG(1,"Using Output Format RVdouble\n");
        AA = new OutputRVdouble();

    } else if (format == "Pack9") {
        STDLOG(1,"Using Output Format Pack9\n");
        AA = new OutputPacked<9>();

    } else if (format == "Pack" or format == "Pack14") {
        STDLOG(1,"Using Output Format Pack14\n");
        AA = new OutputPacked<14>();
        
    } else if (format == "Heitmann") {
        STDLOG(1,"Using Output Format Heitmann\n");
        AA = new OutputHeitmann();

    } else if (format == "RVdoublePID") {
        STDLOG(1,"Using Output Format RVdoublePID\n");
        AA = new OutputRVdoublePID();

    } else if (format == "RVZel") {
        STDLOG(1,"Using Output Format RVZel\n");
        AA = new OutputRVZel();
        
    }
    else {
        QUIT("Unrecognized case: OutputFormat = {:s}\n", format);
    }

    return AA;
}

AppendArena *get_PID_AA_by_format(const std::string &format){
    AppendArena *PID_AA;
    if (format == "Pack9") {
        PID_AA = new OutputPID(); 
        STDLOG(2, "Chose PID timeslice append arena to complement pack9 RVs.\n");
    }
    else {
        PID_AA = NULL; 
        STDLOG(2, "No pack 9 timeslice requested; setting PID timeslice append area to NULL.\n");
    }
    return PID_AA;
}

void WriteHeaderFile(const fs::path &fn){
	std::ofstream headerfile;
	headerfile.open(fn);
	headerfile << P.header();
	headerfile << ReadState.header();
	headerfile.close();
}

uint64 Output_TimeSlice(int slab, FLOAT unkickfactor) {
    AppendArena *AA, *PID_AA;

    AA = get_AA_by_format(P.OutputFormat);
    PID_AA = get_PID_AA_by_format(P.OutputFormat);
    
    // Setup the Arena
    int headersize = 1024*1024;
    SB->AllocateSpecificSize(FieldTimeSlice, slab, 
    	   SS->size(slab) * AA->sizeof_particle()
	+ CP->cpd * node_z_size * AA->sizeof_cell() + headersize);

    AA->initialize(FieldTimeSlice, slab, CP->cpd, ReadState.VelZSpace_to_Canonical);

    if (PID_AA != NULL) {
        SB->AllocateSpecificSize(FieldTimeSlicePIDs, slab, 
           SS->size(slab) * PID_AA->sizeof_particle()
    + CP->cpd * node_z_size * PID_AA->sizeof_cell() + headersize);
        
        PID_AA->initialize(FieldTimeSlicePIDs, slab, CP->cpd, ReadState.VelZSpace_to_Canonical);
    }

    STDLOG(4,"Writing header\n");

    // Write the header to its own file
    if(slab == 0){
        WriteHeaderFile(P.OutputDirectory / fmt::format("slice{:5.3f}/header", ReadState.Redshift));
    }

    STDLOG(4,"Adding header to slab file\n");  
    // and also add the header to the slab file
    if (!P.OmitOutputHeader) {
        AA->addheader(P.header());
        AA->addheader(ReadState.header());
        AA->addheader("\nOutputType = \"FieldTimeSlice\"\n");
        AA->addheader(fmt::format("SlabNumber = {:d}\n", slab));
        // For sanity, be careful that the previous lines end with a \n!
        AA->finalize_header();
    }
    STDLOG(2,"Scanning through cells\n");

    // Now scan through the cells
    integer3 ij(slab,0,node_z_start);
    uint64 n_added = 0;
    #pragma omp parallel for schedule(static) reduction(+:n_added)
    for (int y=0; y<CP->cpd; y++) {
        integer3 ijk = ij; ijk.y = y;
        // We are required to provide an offset in bytes for this pencil's portion of the buffer.
        // We can use CellInfo(ijk)->startindex, even though the position slabs have ghosts
        // because we specifically want the cumulative counts without ghosts
    	long long int start = CP->CellInfo(ijk)->startindex;   // Assumes cells are packed in order in the slab
        AA->start_pencil(y, start*AA->sizeof_particle() + AA->sizeof_cell()*node_z_size*y);
        if (PID_AA!=NULL) PID_AA->start_pencil(y, start*PID_AA->sizeof_particle());

        for (ijk.z = node_z_start; ijk.z < node_z_start + node_z_size; ijk.z++) {
            Cell c = CP->GetCell(ijk);
            // We sometimes use the maximum velocity to scale.
            // But we do not yet have the global velocity (slab max will be set in Finish,
            // while the global max has to wait for all slabs to be done).
            // What is available after the kick is the max_component_velocity in each cell.
            FLOAT vscale = c.ci->max_component_velocity/ReadState.VelZSpace_to_Canonical;	
            // The maximum velocity of this cell, converted to ZSpace unit-box units.
            // Start the cell
            AA->addcell(y, ijk, vscale);
            if (PID_AA != NULL) PID_AA->addcell(y, ijk, vscale);
            // Now pack the particles
            accstruct *acc = CP->AccCell(ijk);
            for (int p=0;p<c.count();p++) {
                velstruct vel = (c.vel[p] - static_cast<FLOAT3>(acc[p])*unkickfactor);    // We supply in code units
                // Detail: we write particles with their L0 bits intact.  So if we want to run a non-group-finding step
                // after a group-finding step (e.g. for debugging), we need to know that we can ignore the L0 bit
                if(GFC == NULL || !c.aux[p].is_L0()){
                    c.aux[p].set_compressed_density(acc[p].w);
                    AA->addparticle(y, c.pos[p], vel, c.aux[p]);
                    if (PID_AA != NULL) PID_AA->addparticle(y, c.pos[p], vel, c.aux[p]);
                    n_added++;
                }
            }
            AA->endcell(y);
            if (PID_AA != NULL) PID_AA->endcell(y); 
        }
    }

    STDLOG(2,"Resizing slab\n");
    SB->ResizeSlab(FieldTimeSlice, slab, AA->finalize_arena());
    STDLOG(4,"StoreArenaNonBlocking\n");
    // Write out this time slice
    SB->StoreArenaNonBlocking(FieldTimeSlice, slab);
    delete AA;

    if (PID_AA != NULL) { 
        STDLOG(2,"Resizing slab\n");
        SB->ResizeSlab(FieldTimeSlicePIDs, slab, PID_AA->finalize_arena());
        STDLOG(4,"StoreArenaNonBlocking\n");
        // Write out this time slice
        SB->StoreArenaNonBlocking(FieldTimeSlicePIDs, slab);
        delete PID_AA;
    }
    
    return n_added;
}

