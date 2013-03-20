// a second set of data
typedef struct {
    double3 kvec;
    double3 phase;
    double Across;
    char datafile[1024];
} spiralparam;

// and a second specialization
template <>
void ParseHeader::register_vars(spiralparam &P) {
    installvector("kvec", &(P.kvec.x), 3, 1, MUST_DEFINE);
    installvector("phase", &(P.phase.x), 3, 1, MUST_DEFINE);
    installscalar("Across", P.Across, MUST_DEFINE);
    installscalar("datafile", P.datafile, MUST_DEFINE);
}
