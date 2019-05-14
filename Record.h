#ifndef RECORD_H
#define RECORD_H

enum class HR {U,V,T, Nx, Ny, Nz, Px, Py, Pz};
enum class HI {M, T};

struct ScatterRecord
{
    vec3 o, dir;
    bool is_specular;
    vec3 attenuation;

    int pdfIdx;
};
#endif
