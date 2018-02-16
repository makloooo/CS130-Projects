#include "reflective_shader.h"
#include "ray.h"
#include "render_world.h"

vec3 Reflective_Shader::
Shade_Surface(const Ray& ray,const vec3& intersection_point,
    const vec3& same_side_normal,int recursion_depth) const
{
	vec3 reflected_color;
	
	// what does this even reference to? itself?
	vec3 shader_color = shader->Shade_Surface(ray, intersection_point, same_side_normal, recursion_depth);

	// skeptical about this
	vec3 reflect_direction = (ray.direction - (2*dot(ray.direction,same_side_normal)*same_side_normal)).normalized();
	Ray reflected_ray = Ray(intersection_point, reflect_direction);

	/*
	std::cout << "Our ray (" << ray.endpoint << ") with direction <" << ray.direction
						<< "> has a reflection (" << reflected_ray.endpoint << ") with direction <"
						<< reflected_ray.direction << ">\n";
	*/

	reflected_color = world.Cast_Ray(reflected_ray, recursion_depth+1);

	return reflectivity*reflected_color+(1-reflectivity)*shader_color;
}
