//==================================================================================================
// Written in 2016 by Peter Shirley <ptrshrl@gmail.com>
//
// To the extent possible under law, the author(s) have dedicated all copyright and related and
// neighboring rights to this software to the public domain worldwide. This software is distributed
// without any warranty.
//
// You should have received a copy (see file COPYING.txt) of the CC0 Public Domain Dedication along
// with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
//==================================================================================================

#ifndef AABBH
#define AABBH
#include "ray.h"
#include "hitable.h"

inline float ffmin(float a, float b) { return a < b ? a : b; }
inline float ffmax(float a, float b) { return a > b ? a : b; }

class aabb {
    public:
        aabb() {}
        aabb(const vec3& a, const vec3& b) { _min = a; _max = b;}

        vec3 min() const {return _min; }
        vec3 max() const {return _max; }

        bool hit(const ray& r, float tmin, float tmax) const {
            for (int a = 0; a < 3; a++) {
                float t0 = ffmin((_min[a] - r.origin()[a]) / r.direction()[a],
                                (_max[a] - r.origin()[a]) / r.direction()[a]);
                float t1 = ffmax((_min[a] - r.origin()[a]) / r.direction()[a],
                                (_max[a] - r.origin()[a]) / r.direction()[a]);
                tmin = ffmax(t0, tmin);
                tmax = ffmin(t1, tmax);
                if (tmax <= tmin)
                    return false;
            }
            return true;
        }

        float area() const {
               float a = _max[0] - _min[0];
               float b = _max[1] - _min[1];
               float c = _max[2] - _min[2];
               return 2*(a*b + b*c + c*a);
        }

        int longest_axis() const {
               float a = _max[0] - _min[0];
               float b = _max[1] - _min[1];
               float c = _max[2] - _min[2];
               if (a > b && a > c)
                   return 0;
               else if (b > c)
                   return 1;
               else
                   return 2;
        }

        vec3 _min;
        vec3 _max;
};

aabb surrounding_box(aabb box0, aabb box1) {
    vec3 small( fmin(box0.min()[0], box1.min()[0]),
                fmin(box0.min()[1], box1.min()[1]),
                fmin(box0.min()[2], box1.min()[2]));
    vec3 big  ( fmax(box0.max()[0], box1.max()[0]),
                fmax(box0.max()[1], box1.max()[1]),
                fmax(box0.max()[2], box1.max()[2]));
    return aabb(small,big);
}


#endif

