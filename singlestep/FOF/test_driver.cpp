// This is a code to provide a test for the group finding.

// There are a couple of other unit tests:
// make fof_sublist
// make slab_accum

#include "test_driver_header.cpp"
#include "test_driver_particles.cpp"

#define ASCII_TEST_OUTPUT

#include "groupfinding.cpp"

int main(int argc, char *argv[]) {
    setup_log();
    fmt::print("Running with {:d} threads\n", omp_get_max_threads());
    int cpd = 9;
    CP = new CellParticles;
    // CP->MakeParticles(cpd,3000);
    // GFC = new GroupFindingControl(0.04, CP->cpd, CP->invcpd, 2, CP->np);
    float3 offset(0.0,0.0,0.0);
    if (argc==4) offset = float3(atof(argv[1])/cpd,atof(argv[2])/cpd,atof(argv[3])/cpd);
    CP->ReadParticles(cpd, "sample.dat", offset);
    GFC = new GroupFindingControl(0.01, 0.007, 0.005, CP->cpd, CP->invcpd, 2, 20, CP->np);

    fmt::print("CPD = {:d}, b = {:f}, b*CPD = {:f}\n", GFC->cpd, GFC->linking_length, GFC->linking_length*GFC->cpd);

    for (int s=0; s<CP->cpd; s++) GFC->ConstructCellGroups(s);
    fmt::print("Finished Constructing Cell Groups\n"); fflush(NULL);

    for (int s=0; s<CP->cpd; s++) FindGroupLinks(s);
    fmt::print("Finished Finding Group Links\n"); fflush(NULL);

    for (int s=0; s<CP->cpd; s++) FindAndProcessGlobalGroups(s);
    fmt::print("Finished Global Groups\n"); fflush(NULL);

    for (int s=0; s<CP->cpd; s++) GFC->DestroyCellGroups(s);
    fmt::print("Finished Destroying Cell Groups\n"); fflush(NULL);

    GFC->report();

    delete CP;
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
