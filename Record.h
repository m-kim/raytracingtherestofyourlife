#ifndef RECORD_H
#define RECORD_H

struct HitRecord
{
  vtkm::Id hitIdx;
    float t;
    float u;
    float v;
    vec3 p;
    vec3 normal;

    int texId, matId;
};


struct ScatterRecord
{
    ray specular_ray;
    bool is_specular;
    vec3 attenuation;

    int pdfIdx;
};
#endif
