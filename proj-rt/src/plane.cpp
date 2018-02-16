#include "plane.h"
#include "ray.h"
#include <cfloat>


// Intersect with the half space defined by the plane.  The plane's normal
// points outside.  If the ray starts on the "inside" side of the plane, be sure
// to record a hit with t=0 as the first entry in hits.
bool Plane::
Intersection(const Ray& ray, std::vector<Hit>& hits) const
{
    // implicit: dot(x-x0,n)=0
    vec3 u = ray.endpoint;
    vec3 w = ray.direction;
    
    if (dot(normal, w) == 0) return false; // is parallel to plane, ignore inf solutions
    
    // if intersect return true
    float t = dot(-normal, (u - x1))/dot(normal, w);
    //std::cout << "t: " << t << '\n';
    if (t < 0) return false;
    
    hits.emplace_back(this, t, false);
    
    return true;
}

vec3 Plane::
Normal(const vec3& point) const
{
    return normal;
}
