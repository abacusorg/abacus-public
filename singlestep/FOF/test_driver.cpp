// This is a code to provide a test for the group finding.

// There are a couple of other unit tests:
// g++ -DTEST -fopenmp -lgomp -O2 fof_sublist.cpp 
// g++ -DTEST -fopenmp -lgomp -O2 slab_accum.cpp

#include "test_driver_header.cpp"
#include "test_driver_particles.cpp"

#define ASCII_TEST_OUTPUT

#include "groupfinding.cpp"

int main(int argc, char *argv[]) {
    setup_log();
    printf("Running with %d threads\n", omp_get_max_threads());
    int cpd = 9;
    PP = new Particles;
    // PP->MakeParticles(cpd,3000);
    // GFC = new GroupFindingControl(0.04, PP->cpd, PP->invcpd, 2, PP->np);
    float3 offset(0.0,0.0,0.0);
    if (argc==4) offset = float3(atof(argv[1])/cpd,atof(argv[2])/cpd,atof(argv[3])/cpd);
    PP->ReadParticles(cpd, "sample.dat", offset);
    GFC = new GroupFindingControl(0.01, 0.007, 0.005, PP->cpd, PP->invcpd, 2, 20, PP->np);

    printf("CPD = %d, b = %f, b*CPD = %f\n", GFC->cpd, GFC->linking_length, GFC->linking_length*GFC->cpd);

    for (int s=0; s<PP->cpd; s++) GFC->ConstructCellGroups(s);
    printf("Finished Constructing Cell Groups\n"); fflush(NULL);

    for (int s=0; s<PP->cpd; s++) FindGroupLinks(s);
    printf("Finished Finding Group Links\n"); fflush(NULL);

    for (int s=0; s<PP->cpd; s++) FindAndProcessGlobalGroups(s);
    printf("Finished Global Groups\n"); fflush(NULL);

    for (int s=0; s<PP->cpd; s++) GFC->DestroyCellGroups(s);
    printf("Finished Destroying Cell Groups\n"); fflush(NULL);

    GFC->report();

    delete PP;
    delete GFC;
    return 0;
}

/* Here is the reference output:
Made 3000 particles, often 4 per cell, 88 in the last
CPD = 9, b = 0.040000, b*CPD = 0.360000
Finished Constructing Cell Groups
Finished Finding Group Links
Finished Global Groups
Finished Destroying Cell Groups
Found 1929 cell groups (including boundary singlets)
Used 2919 pseudoParticles, 3260 faceParticles, 251 faceGroups
Found 502 links between groups.
Found 553 global groups
*/
