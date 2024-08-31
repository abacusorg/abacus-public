#define THISCODE std::string("Time-stamp: <DumpSpiral.cc on Saturday, 19 November, 2011 at 21:28:38 EST (bob)>")

/*
compile with:
icc -o dumpspiral -I../ParseHeader -I../ParseHeader/util  -I../util -I../include -I.. DumpSpiral.cc ../util/stringutil.cc -L ../ParseHeader -lparseheader

 */

#include <filesystem>
namespace fs = std::filesystem;

#include "header.h"
#include "util.h"
#include "threevector.hh"
#include "ParseHeader.hh"
#include "timestamp.cc"
#include "stringutil.cc"

typedef std::complex<double> Complex;
typedef ThreeVector<double> double3;
typedef ThreeVector<float> float3;
typedef ThreeVector<int> int3;
typedef long long int LLI;
typedef unsigned long long int ULLI;

using namespace std;

#include "particles.h"
#include "Grid.cc"
#include "Spiral.cc"
#include "spiralparameters.cc"
#include "spiralspec.cc"
#include "dimensionalization.cc"
#include "Cosmology.cc"

void Usage(char *progname) {
    fmt::print("Usage: {:s}: -i infile -r rotfile\n", progname);
    abort();
}

void ParseCommandLine(int argc, char **argv, fs::path &infile, fs::path &rotfile) {
    int c;
    while ((c = getopt (argc, argv, "i:r:")) != -1)
        switch (c)
        {
        case 'i':
            infile = optarg;
            break;
        case 'r':
            rotfile = optarg;
            break;
        case '?':
            Usage(argv[0]);
        default:
            abort ();
        }
     
    int error = 0;
    for (int index = optind; index < argc; index++) {
        fmt::print(std::cerr, "Non-option argument {:s}\n", argv[index]);
        error = 1;
    }
    if(error || infile.empty() || rotfile.empty()) Usage(argv[0]);
}

int main(int argc, char **argv) {

    fs::path infile, rotfile;
    ParseCommandLine(argc, argv, infile, rotfile);
    
    ParseHeader PH;

    codeparam cp;
    PH.register_vars(cp);
    spiralparam sp;
    PH.register_vars(sp);
    HeaderStream instream(infile);
    PH.ReadHeader(instream);

    Spiral SP(cp.N1D, sp.kvec, sp.phase);

    int np;
    double anow;
    fread(&anow, sizeof(double), 1, instream.fp);
    fread(&np, sizeof(LLI), 1, instream.fp);
    assert( np == cp.N1D*cp.N1D*cp.N1D );
    particlestructure *ps = new particlestructure[np];
    fread(ps, sizeof(particlestructure), np, instream.fp);
    instream.Close();

    MyCosmology cosmo;
    cosmo.Omega_m = cp.Omega_M;
    cosmo.Omega_K = cp.Omega_K;
    cosmo.Omega_DE = cp.Omega_DE;
    cosmo.H0 = cp.H0;
    cosmo.w0 = cp.w0;
    cosmo.wa = cp.wa;

    Cosmology C(1e-3,cosmo);

    //ConvertVpecToDimensionless(ps, np, cp, anow);
    ConvertVzToDimensionless(ps, np, cp, anow, C.H(anow));

    SP.Project(ps, np);

    ConvertDimensionlessToVpec(ps, np, cp, anow);
    //ConvertDimensionlessToVz(ps, np, cp, anow, C.H(anow));
    fmt::print("H({}) = {}\n", anow, C.H(anow));

    FILE *out = OpenForWrite(rotfile, true);
    WriteHStream(out, "##\n", "# ");
    WriteHStream(out, " Input parameters from " + instream.name + ":\n", "# ");
    WriteHStream(out, "##\n", "# ");
    WriteHStream(out, instream, "# ");
    WriteHStream(out, "\n", "# " );
    WriteHStream(out, "# by: " + THISCODE + "\n");
    //    FinalizeHeader(out); // don't want the EOH token here...

    for(int p=0; p<np; p++) {
        fmt::print(out, "{: e} {: e}\n", ps[p].position.x, ps[p].velocity.x);
    }

    fclose(out);
}
