#include "light.h"
#include "phong_shader.h"
#include "ray.h"
#include "render_world.h"
#include "object.h"

vec3 Phong_Shader::
Shade_Surface(const Ray& ray,const vec3& intersection_point,
    const vec3& same_side_normal,int recursion_depth) const
{
	// This executes per fragment
  //std::cout << "New fragment\n";
  // ray = view ray, intersecton point = intersect with view ray  

  vec3 color;
	
	vec3 ambient = world.ambient_color*world.ambient_intensity*color_ambient;
	vec3 specular = vec3(0,0,0);
	vec3 diffuse = vec3(0,0,0);

	for (unsigned i = 0; i < world.lights.size(); ++i) {
		// l is light direction
		Light* light = world.lights.at(i); // do for each light
		vec3 l = light->position - intersection_point; // vector
		vec3 d = l.normalized();
		if (world.enable_shadows) { // FIXME: we dont do anything with the intersection point here
			// if further along the ray, then make shadow
			Hit hit;
			Ray l_ray = Ray(light->position, -d);
			//std::cout << "position: " << light->position << " | -d: " << -d << '\n';
			world.Closest_Intersection(l_ray, hit);
			if (hit.ray_exiting || l_ray.Point(hit.t) != intersection_point) {
				continue;
			}
			//std::cout << l_ray.Point(hit.t) << " | " << intersection_point << '\n';
		}
		// std::cout << "Moving on!\n";
		// is it behind anything?
		double intensity = std::max(0.0, dot(d, same_side_normal));

		//std::cout << "Direction: " << d << '\n';
		//std::cout << "Normal: " << same_side_normal << '\n';
		//std::cout << "Our intensity : " << intensity << '\n'; // '|' << same_side_normal << '\n';
		//std::cout << "Our Light brightness : " << light->Emitted_Light(ray)/255 << '\n';

		vec3 light_color = light->Emitted_Light(ray)/l.magnitude_squared();

		diffuse += color_diffuse*intensity*light_color;
		//std::cout << "Our diffuse: " << diff << '\n';

		vec3 R = (2*dot(d,same_side_normal)*same_side_normal-d).normalized();
		vec3 C = (world.camera.position - intersection_point).normalized();
		double blinding = std::max(0.0, dot(R,C));

		specular += color_specular*pow(blinding,specular_power)*light_color;
	}
	color = ambient + diffuse + specular;
  return color;
}
