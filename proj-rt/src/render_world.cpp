#include "render_world.h"
#include "flat_shader.h"
#include "object.h"
#include "light.h"
#include "ray.h"

Render_World::Render_World()
    :background_shader(0),ambient_intensity(0),enable_shadows(true),
    recursion_depth_limit(3)
{}

Render_World::~Render_World()
{
    delete background_shader;
    for(size_t i=0;i<objects.size();i++) delete objects[i];
    for(size_t i=0;i<lights.size();i++) delete lights[i];
}

// Find the closest object of intersection and return the object that was
// intersected.  Record the Hit structure in hit.  If no intersection occurred,
// return NULL.  Note that in the case of a Boolean, the object returned will be
// the Boolean, but the object stored in hit will be the underlying primitive.
// Any intersection with t<=small_t should be ignored.
Object* Render_World::Closest_Intersection(const Ray& ray, Hit& hit)
{
    float min_t = INFINITY;
    
    Object* closest_object = 0;
    // For each object in objects
    for (unsigned int i = 0; i < objects.size(); ++i) {
		// create empty list of hits
		std::vector<Hit> hits;
		// use object.intersect to fill hits
		objects.at(i)->Intersection(ray, hits);
		// for each h in list of hits
		for (unsigned int j = 0; j < hits.size(); ++j) {
			// if h is closest so far (with smallest t, that is larger than small_t)
			//std::cout << "We hit at : " << hits.at(j).t << '\n';
			if (hits.at(j).t < min_t && hits.at(j).t > small_t) {
				// set object as closest object, set hit to h and update min_t
				closest_object = objects.at(i);
				hit = hits.at(j);
				min_t = hits.at(j).t;
			}
		}
		hits.clear();
	}
	return closest_object; // return the closest obj
}

// set up the initial view ray and call
void Render_World::Render_Pixel(const ivec2& pixel_index)
{
    Ray ray; // set up the initial view ray here
    ray.endpoint = camera.position;
    ray.direction = (camera.World_Position(pixel_index)-camera.position).normalized();
    //std::cout << "Our pixel : " << pixel_index << '\n';
    vec3 color=Cast_Ray(ray,1);
    camera.Set_Pixel(pixel_index,Pixel_Color(color));
    
}

void Render_World::Render()
{
    for(int j=0;j<camera.number_pixels[1];j++)
        for(int i=0;i<camera.number_pixels[0];i++)
            Render_Pixel(ivec2(i,j));
     
    
    //std::cout << "Total number of hits: " << debug_hits << '\n';
    //std::cout << "Total number of non-hits: " << debug_nohits << '\n';
}

// cast ray and return the color of the closest intersected surface point,
// or the background color if there is no object intersection
vec3 Render_World::Cast_Ray(const Ray& ray,int recursion_depth)
{
	vec3 color;
	vec3 intersect;
	vec3 normal;
	
	if (recursion_depth > recursion_depth_limit) 
		return background_shader->Shade_Surface(ray, intersect, normal, 1);

	Hit hit;
	// if hit an object, set color to that
	Object* obj = Closest_Intersection(ray, hit);
	if (obj) {
		++debug_hits;
		intersect = ray.Point(hit.t);
		color = obj->material_shader->Shade_Surface(ray, intersect, obj->Normal(intersect), recursion_depth);
	}
	else {
		++debug_nohits;
		color = background_shader->Shade_Surface(ray, intersect, normal, 1); // if it returns null, no intersection
	}
	
	return color; // later distort color with other reflected surfaces
}
