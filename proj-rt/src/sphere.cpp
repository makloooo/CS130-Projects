#include "sphere.h"
#include "ray.h"


// Determine if the ray intersects with the sphere
bool Sphere::Intersection(const Ray& ray, std::vector<Hit>& hits) const
{
    // TODO
	
    vec3 v = ray.endpoint - this->center;
    vec3 w = ray.direction;
    
    //std::cout << "Our sphere's center : " << this->center << '\n';
    //std::cout << "Our vectors : " << v << '|' << w << '\n';
    
    double D = pow(dot(w, v),2) - (v.magnitude_squared()-pow(radius,2));
    //if (D >= 0) std::cout << "Our determinant: " << D << '\n';
    //std::cout << "Our dot(w,v): " << dot(w,v) << '\n';
    if (D < 0) return false;
    else if (D == 0) {
		double t = -dot(w,v);
		hits.emplace_back(this, t, false);
		hits.emplace_back(this, t, true);
	}
	else if (D > 0) {
		double t0 = -dot(w,v)-sqrt(D);
		double t1 = -dot(w,v)+sqrt(D);
		if (t0 >= t1) { double tmp = t0; t0 = t1; t1 = tmp; }
		if (t0 < 0) t0 = 0;
		hits.emplace_back(this, t0, false);
		hits.emplace_back(this, t1, true);
	}
	
	return true;
}

vec3 Sphere::Normal(const vec3& point) const
{
    vec3 normal;
    // TODO: set the normal
    normal = (point - center).normalized();
    return normal;
}
