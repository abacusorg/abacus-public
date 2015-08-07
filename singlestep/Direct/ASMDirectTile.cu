#ifndef DOUBLEPRECISION
#ifndef FLOAT
#define FLOAT float
#endif
#else
#ifndef FLOAT
#define FLOAT double
#endif
#endif


extern "C" __device__  void FullDirectTile(
        FLOAT *BlockSources_x,
        FLOAT *BlockSources_y,
        FLOAT *BlockSources_z,
        FLOAT *a_x,      FLOAT *a_y,       FLOAT *a_z,
        FLOAT *sink_x,   FLOAT *sink_y,   FLOAT *sink_z,
        FLOAT *eps2)

#ifdef DUMMYPTX
{
    return;
}
#else
;
#endif
