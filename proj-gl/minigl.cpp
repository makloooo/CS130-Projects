/**
 * minigl.cpp
 * -------------------------------
 * Implement miniGL here.
 *
 * You may include minigl.h and any of the standard C++ libraries.
 * No other includes are permitted.  Other preprocessing directives
 * are also not permitted.  These requirements are strictly
 * enforced.  Be sure to run a test grading to make sure your file
 * passes the sanity tests.
 *
 * The behavior of the routines your are implenting is documented here:
 * https://www.opengl.org/sdk/docs/man2/
 * Note that you will only be implementing a subset of this.  In particular,
 * you only need to implement enough to pass the tests in the suite.
 */

#include "minigl.h"
#include "vec.h"
#include "mat.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>
#include <cstdio>
#include <limits>
#include <cfloat>

using namespace std;

/**
 * User defined data variables/structures
 */

bool mgl_Enabled = false;

struct vertex {
	vec3 rgb;
	vec4 pos;
	MGLfloat z;
	MGLfloat w;
	vertex() {}
	vertex(vec3 rgb, vec4 pos, MGLfloat z, MGLfloat w) {
		this->rgb = rgb;
		this->pos = pos;
		this->z = z;
		this->w = w;
	}
};

struct triangle {
	vertex a, b, c;
	float area;
	triangle() {}
	triangle(vertex a, vertex b, vertex c) {
		this->a = a;
		this->b = b;
		this->c = c;
	}
};

struct dataOut {
	MGLfloat depth;
	MGLfloat alpha;
	MGLfloat beta;
	MGLfloat gamma;
	dataOut(MGLfloat d, MGLfloat a, MGLfloat b, MGLfloat g) {
		this->depth = d;
		this->alpha = a;
		this->beta = b;
		this->gamma = g;
	}
};

/* Declarations */
vector<triangle> triangles;
vector<vertex> vertices;
vec3 curr_color;
MGLpoly_mode curr_geometry;
MGLmatrix_mode curr_matrix_mode = MGL_MODELVIEW;

enum Bounds {LEFT, RIGHT, BOTTOM, TOP, NEAR, FAR};
MGLfloat bound[6] = {0, 0, 0, 0, 1, -1};
MGLint window[6] = {-1, 1, -1, 1, 1, -1};

MGLfloat near_plane = 0.0;
MGLfloat far_plane = 0.0;

mat4 modelview_matrix;
vector<mat4> modelview_stack;
mat4 proj_matrix;
vector<mat4> proj_stack;

vector<mat4>* curr_stack;
mat4* curr_matrix = &modelview_matrix;

// Epsilon function adapted from floating-point-gui.de/errors/comparison
bool eq(float a, float b, float epsilon=FLT_EPSILON) {
	a = abs(a); b = abs(b); float diff = a-b;
	if (a == b) return true;
	else if (a == 0.0 || b == 0.0 || diff < FLT_MIN) return diff < (epsilon * FLT_MIN);
	else return diff / min(a+b, FLT_MAX) < epsilon;
}

bool inTriangle(float alpha, float beta, float gamma) {
	if (eq(alpha, 0)) alpha = 0;
	if (eq(beta, 0)) beta = 0;
	if (eq(gamma, 0)) gamma = 0;

	if (alpha < 0) return false;
	if (beta < 0) return false;
	if (gamma < 0) return false;
	return true;
}

bool inViewPlane(vec3 p) {
	for (int i = 0; i < 3; ++i) if (eq(p[i], int(p[i]))) p[i] = int(p[i]);
	if (p[2] > bound[NEAR] || p[2] < bound[FAR]) return false;
	if (p[1] < bound[BOTTOM] || p[1] > bound[TOP]) return false;
	if (p[0] < bound[LEFT] || p[0] > bound[RIGHT]) return false;
	return true;
}

bool inWindow(vec4 p) {
	if (p[2] <= window[NEAR] && p[2] >= window[FAR] &&
			p[1] >= window[BOTTOM] && p[1] <= window[TOP] &&
			p[0] >= window[LEFT] && p[0] <= window[RIGHT]) return true;
	return false;
}

bool inWindow(triangle t) {
	// Simple case: At least one vertex enclosed in window
	if (inWindow(t.a.pos) || inWindow(t.b.pos) || inWindow(t.c.pos)) return true;

	// Other case: All triangle vertices outside of window, but 
	// lines between vertices intersect window
	//
	// If both coordinates of a certain dimension are beyond the same 
	// bound, they will never intersect the window
	// Implementing an elegant way to iterate through enums is
	// too much. Just brute force.
	
	vec4 x0 = t.a.pos; vec4 x1 = t.b.pos; vec4 x2 = t.c.pos;
	if (!(x0[0] < window[LEFT] && x1[0] < window[LEFT] && x2[0] < window[LEFT]) ||
			!(x0[0] > window[RIGHT] && x1[0] > window[RIGHT] && x2[0] > window[RIGHT]) ||
			!(x0[1] > window[TOP] && x1[1] > window[TOP] && x2[0] > window[TOP]) ||
			!(x0[1] < window[BOTTOM] && x1[1] < window[BOTTOM] && x2[2] < window[BOTTOM]) ||
			!(x0[2] > window[NEAR] && x1[2] > window[NEAR] && x2[2] < window[NEAR]) ||
			!(x0[2] < window[FAR] && x1[2] < window[FAR] && x2[2] < window[FAR])) return true;

	return false;
}

vec4 screen(vec4 v, MGLsize width, MGLsize height) {
	vec4 result = vec4((v[0]+1)*width/2, (v[1]+1)*height/2, v[2], v[3]);
	return result;
}
	
/** Note this isn't the actual area of it, but only used for
 *  the purposes of finding the ratio between two areas
 */
float areaOf(vec4 a, vec4 b, vec4 c) {
	float result = a[0]*(b[1]-c[1]);
	result += a[1]*(c[0]-b[0]);
	result += (b[0]*c[1]-b[1]*c[0]);
	//cout << "Calculated Area: " << result << endl;
	return result;
}

/** Get the z-coordinate of three vertices using correct
 *  perspective interpolation
 */
dataOut perspective_interpolation(triangle t, float alpha, float beta, float gamma) {
	MGLfloat d = alpha*t.a.z/t.a.w + beta*t.b.z/t.b.w + gamma*t.c.z/t.c.w;
	// Color interpolation (This is only a copy, real attributes not changed.
	float k = alpha/t.a.w + beta/t.b.w + gamma/t.c.w;
	alpha = alpha/t.a.w/k; beta = beta/t.b.w/k; gamma = gamma/t.c.w/k;

	dataOut info(d, alpha, beta, gamma);
	return info;
}

/**
 * Standard macro to report errors
 */
inline void MGL_ERROR(const char* description) {
    printf("%s\n", description);
    exit(1);
}


/**
 * Read pixel data starting with the pixel at coordinates
 * (0, 0), up to (width,  height), into the array
 * pointed to by data.  The boundaries are lower-inclusive,
 * that is, a call with width = height = 1 would just read
 * the pixel at (0, 0).
 *
 * Rasterization and z-buffering should be performed when
 * this function is called, so that the data array is filled
 * with the actual pixel values that should be displayed on
 * the two-dimensional screen.
 */
void mglReadPixels(MGLsize width,
                   MGLsize height,
                   MGLpixel* data)
{
	cout << "Begin reading pixels\n";
	
	bound[LEFT] = 0;
	bound[RIGHT] = width-1;
	bound[BOTTOM] = 0;
	bound[TOP] = height-1;
	
	cout << "Bounds: " << bound[LEFT] << '|' << bound[RIGHT] << '|' 
		 << bound[BOTTOM] << '|' << bound[TOP] << '|' << bound[NEAR] 
		 << '|' << bound[FAR] << endl;
	
	float depth[width][height];
	for (MGLsize i = 0; i < width; ++i) {
		for (MGLsize j = 0; j < height; ++j) {
			depth[i][j] = +INFINITY;
		}
	}

	// width is basically (0, width-1)
	// height is (0, height-1)
	// iterate through vertices, draw pixel per triangle
	for (unsigned int i = 0; i < triangles.size(); ++i) {
		triangle t = triangles.at(i);
		
		vec4 a = t.a.pos;
		vec4 b = t.b.pos;
		vec4 c = t.c.pos;

		cout << "a: " << a << endl;
		cout << "b: " << b << endl;
		cout << "c: " << c << endl;

		// When we clip we use z values of actual and not projected object

		/* This is an already working solution, saved just in case.*/
		// screen projection
		a = screen(a, width, height);
		b = screen(b, width, height);
		c = screen(c, width, height);
		cout << "Screen projection done\n";
		
		cout << "a: " << a << endl;
		cout << "b: " << b << endl;
		cout << "c: " << c << endl;

		cout << "our window width is: " << width << endl;
		cout << "our window length is: " << height << endl;	
	
		// finding bounding box
		int x_min = min(a[0], min(b[0], c[0]));
		int x_max = max(a[0], max(b[0], c[0]));
		int y_min = min(a[1], min(b[1], c[1]));
		int y_max = max(a[1], max(b[1], c[1]));
		cout << "Found bounding box @ " << x_min << '-' << x_max << ' ' << y_min << '-' << y_max << '\n';

		// Snap these back to NDC vertices
	  a = t.a.pos;
		b = t.b.pos;
		c = t.c.pos;

		t.area = areaOf(a, b, c);

		float dx = 2.0/width;
		float dy = 2.0/height;
		cout << "\"Area\" of triangle: " << t.area << '\n';
		for (int x = x_min; x <= x_max; ++x) {
			if (x < 0 || x >= int(width)) continue; // If pixel isn't going to be in the window, don't bother
			for (int y = y_min; y <= y_max; ++y) {
				if (y < 0 || y >= int(height)) continue;
				vec4 p = vec4((x+0.5)*dx-1, (y+0.5)*dy-1, depth[x][y], 1);
				float alpha = areaOf(p, b, c)/t.area;
				float beta = areaOf(a, p, c)/t.area;
				float gamma = areaOf(a, b, p)/t.area;
	
				if (inTriangle(alpha, beta, gamma)) {
					dataOut d = perspective_interpolation(t, alpha, beta, gamma);
					if (inViewPlane(vec3(x, y, d.depth)) && d.depth < depth[x][y]) {
						//cout << "Drawing Pixel\n";
						//cout << "Replacing depth " << depth[x][y] << " with " << d.depth << endl;
						depth[x][y] = d.depth;
						int offs = x+y*width;
						vec3 rgb = t.a.rgb*d.alpha*255 + t.b.rgb*d.beta*255 + t.c.rgb*d.gamma*255;
						*(data+offs) = Make_Pixel(rgb[0], rgb[1], rgb[2]);
					}
				}
			}
		}
	}

	cout << "Pixels Read\n";
}

/**
 * Start specifying the vertices for a group of primitives,
 * whose type is specified by the given mode.
 */
void mglBegin(MGLpoly_mode mode)
{
	//if (mode != GL_TRIANGLES || mode != MGL_QUADS) MGL_ERROR("GL_INVALID_ENUM");
	if (mgl_Enabled) MGL_ERROR("GL_INVALID_OPERATION");
	mgl_Enabled = true;
	curr_geometry = mode;
	cout << "Begin Mini_GL with mode ";
	(mode == MGL_QUADS) ? cout << "GL_QUADS" : cout << "GL_TRIANGLES";
	cout << '\n';
}

/**
 * Stop specifying the vertices for a group of primitives.
 */
void mglEnd()
{
	if (!mgl_Enabled) MGL_ERROR("GL_INVALID_OPERATION");
	MGLint N = (curr_geometry == MGL_TRIANGLES) ? 3 : 4;
	for (unsigned i = 0; i < vertices.size(); ++i) {
		if ((i+1)%N == 0) {
			triangle tmp = triangle(vertices[i], vertices[i-1], vertices[i-2]);
			if (inWindow(tmp)) triangles.emplace_back(tmp);

			if (N == 4) {
				tmp = triangle(vertices[i], vertices[i-2], vertices[i-3]);
				if (inWindow(tmp)) triangles.emplace_back(tmp);
			}
		}
	}
	cout << "MGLint N = " << N << endl;
	cout << "Total number of triangles: " << triangles.size() << endl;
	mgl_Enabled = false;
	vertices.clear();
	cout << "Ended Mini_GL\n";
}

/**
 * Specify a two-dimensional vertex; the x- and y-coordinates
 * are explicitly specified, while the z-coordinate is assumed
 * to be zero.  Must appear between calls to mglBegin() and
 * mglEnd().
 */
void mglVertex2(MGLfloat x,
                MGLfloat y)
{
	if (!mgl_Enabled) MGL_ERROR("GL_INVALID_OPERATION");
	vec4 pos = vec4(x, y, 0, 1);
	cout << "Modelview Matrix: " << modelview_matrix << endl;
	cout << "Projection Matrix: " << proj_matrix << endl;
	pos = modelview_matrix*pos;
	pos = proj_matrix*pos;
	MGLfloat depth = pos[2];
	MGLfloat w = pos[3];
	pos /= pos[3];
	vertices.emplace_back(curr_color, pos, depth, w);
	cout << "2D-Vertex Pushed: (" << x << ',' << y << ")\n";
	cout << "Transformed Coordinates: (" << pos[0] << ',' << pos[1] << ")\n";
	cout << "Number of Vertices: " << vertices.size() << '\n';
}

/**
 * Specify a three-dimensional vertex.  Must appear between
 * calls to mglBegin() and mglEnd().
 */
void mglVertex3(MGLfloat x,
                MGLfloat y,
                MGLfloat z)
{
	if (!mgl_Enabled) MGL_ERROR("GL_INVALID_OPERATION");
	vec4 pos = vec4(x, y, z, 1);
	// divide in read pixel
	cout << "Modelview Matrix: " << modelview_matrix << endl;
	cout << "Projection Matrix: " << proj_matrix << endl;
	pos = modelview_matrix*pos;
	pos = proj_matrix*pos;
	MGLfloat depth = pos[2];
	MGLfloat w = pos[3];
	pos /= pos[3];
	vertices.emplace_back(curr_color, pos, depth, w);
	cout << "3D-Vertex Pushed: (" << x << ',' << y << ',' << z << ")\n";
	cout << "Transformed Coordinates: (" << pos << ")\n";
	cout << "Number of Vertices: " << vertices.size() << '\n';
}

/**
 * Set the current matrix mode (modelview or projection).
 */
void mglMatrixMode(MGLmatrix_mode mode)
{
	if (mgl_Enabled) MGL_ERROR("MGL_INVALID_OPERATION");
	
	if (mode == MGL_PROJECTION) {
		curr_matrix = &proj_matrix;
		curr_stack = &proj_stack;
	}
	else if (mode == MGL_MODELVIEW) {
		curr_matrix = &modelview_matrix;
		curr_stack = &modelview_stack;
	}
	else MGL_ERROR("MGL_INVALID_ENUM");
	curr_matrix_mode = mode;
	cout << "Current Mode set to ";
	(mode == MGL_PROJECTION) ? cout << "MGL_PROJECTION\n" : cout << "MGL_MODELVIEW\n";
}

/**
 * Push a copy of the current matrix onto the stack for the
 * current matrix mode.
 */
void mglPushMatrix()
{
	curr_stack->emplace_back(*curr_matrix);
	cout << "Current Matrix Pushed\n";
	cout << "Top of Stack is now " << *curr_matrix << '\n';
}

/**
 * Pop the top matrix from the stack for the current matrix
 * mode.
 */
void mglPopMatrix()
{
	cout << "Top of Stack Popped\n";
	*curr_matrix = curr_stack->back();
	curr_stack->pop_back();
	cout << "Current Matrix now " << *curr_matrix << endl;
}

/**
 * Replace the current matrix with the identity.
 */
void mglLoadIdentity()
{
	/*
	if (curr_matrix_mode == MGL_PROJECTION) {}
	else if (curr_matrix_mode == MGL_MODELVIEW) {}
	else MGL_ERROR("GL_INVALID_OPERATION"); */
	*curr_matrix = {1, 0, 0, 0,
				   0, 1, 0, 0,
				   0, 0, 1, 0,
				   0, 0, 0, 1};
	cout << "Identity Loaded\n";
}

/**
 * Replace the current matrix with an arbitrary 4x4 matrix,
 * specified in column-major order.  That is, the matrix
 * is stored as:
 *
 *   ( a0  a4  a8  a12 )
 *   ( a1  a5  a9  a13 )
 *   ( a2  a6  a10 a14 )
 *   ( a3  a7  a11 a15 )
 *
 * where ai is the i'th entry of the array.
 */
void mglLoadMatrix(const MGLfloat *matrix)
{
	// Since matrices don't have copy constructors.
	for (int i = 0; i < curr_matrix->cols(); ++i) {
		for (int j = 0; j < curr_matrix->rows(); ++j) {
			(*curr_matrix)(i, j) = *(matrix+(i+j*curr_matrix->cols()));
		}
	}
}

/**
 * Multiply the current matrix by an arbitrary 4x4 matrix,
 * specified in column-major order.  That is, the matrix
 * is stored as:
 *
 *   ( a0  a4  a8  a12 )
 *   ( a1  a5  a9  a13 )
 *   ( a2  a6  a10 a14 )
 *   ( a3  a7  a11 a15 )
 *
 * where ai is the i'th entry of the array.
 */
void mglMultMatrix(const MGLfloat* matrix)
{
	// Description says only for matrix*matrix operations
	mat4 m; m.make_zero();
	// Since matrices don't have copy constructors.
	for (int i = 0; i < m.cols(); ++i) {
		for (int j = 0; j < m.rows(); ++j) {
			m(i, j) = *(matrix+(i+j*m.cols()));
		}
	}
	
	*curr_matrix = *curr_matrix * m;
	cout << "Result was : " << *curr_matrix << endl;
}

/**
 * Multiply the current matrix by the translation matrix
 * for the translation vector given by (x, y, z).
 */
void mglTranslate(MGLfloat x,
                  MGLfloat y,
                  MGLfloat z)
{
	mat4 translate = {1,0,0,0,
					  0,1,0,0,
					  0,0,1,0,
					  x,y,z,1};
	cout << "Translating Matrix by: " << translate << endl;
	mglMultMatrix(translate.values);
}

/**
 * Multiply the current matrix by the rotation matrix
 * for a rotation of (angle) degrees about the vector
 * from the origin to the point (x, y, z).
 */
void mglRotate(MGLfloat angle,
               MGLfloat x,
               MGLfloat y,
               MGLfloat z)
{
	
	float c = cos(angle*PI/180);
	float s = sin(angle*PI/180);
	vec3 v = vec3(x, y, z).normalized();	
	x = v[0]; y = v[1]; z = v[2];

	cout << x << ' ' << y << ' ' << z << endl;
	cout << angle << ' ' << c << ' ' << s << endl;

	mat4 rotate = {x*x*(1-c)+c, y*x*(1-c)+z*s, x*z*(1-c)-y*s, 0,
                 x*y*(1-c)-z*s, y*y*(1-c)+c, y*z*(1-c)+x*s, 0,
                 x*z*(1-c)+y*s, y*z*(1-c)-x*s, z*z*(1-c)+c, 0,
                 0  , 0  , 0  , 1};
	cout << "Rotating Matrix by: " << rotate << endl;
	mglMultMatrix(rotate.values);
}

/**
 * Multiply the current matrix by the scale matrix
 * for the given scale factors.
 */
void mglScale(MGLfloat x,
              MGLfloat y,
              MGLfloat z)
{
	mat4 scale = {x,0,0,0,
				  0,y,0,0,
				  0,0,z,0,
				  0,0,0,1};
	cout << "Scaling Matrix by: " << scale << endl;
	mglMultMatrix(scale.values);
}

/**
 * Multiply the current matrix by the perspective matrix
 * with the given clipping plane coordinates.
 */
void mglFrustum(MGLfloat left,
                MGLfloat right,
                MGLfloat bottom,
                MGLfloat top,
                MGLfloat near,
                MGLfloat far)
{
	float width = right-left;
	float height = top-bottom;
	float depth = far-near;
	mat4 frustrum = {2*near/width,0,0,0,
					 0,2*near/height,0,0,
					 (right+left)/width,(top+bottom)/height,-(far+near)/depth,-1,
					 0,0,-2*far*near/depth,0};
	cout << "Frustrum: " << frustrum << endl;
	mglMultMatrix(frustrum.values);
	cout << "GL_Frustrum completed\n";
}

/**
 * Multiply the current matrix by the orthographic matrix
 * with the given clipping plane coordinates.
 */
void mglOrtho(MGLfloat left,
              MGLfloat right,
              MGLfloat bottom,
              MGLfloat top,
              MGLfloat near,
              MGLfloat far)
{
	if (left == right || bottom == top || near == far) 
		MGL_ERROR("GL_INVALID_VALUE");
	if (mgl_Enabled) MGL_ERROR("GL_INVALID_OPERATION");
	MGLfloat width = right-left;
	MGLfloat height = top-bottom;
	MGLfloat depth = far-near;
	const mat4 ortho = {2/width, 0, 0, 0,
				        0, 2/height, 0, 0,
				        0, 0, -2/depth, 0,
				        -(right+left)/width, -(top+bottom)/height, -(far+near)/depth, 1};
	cout << "Ortho: " << ortho << '\n';
	mglMultMatrix(ortho.values);
	cout << "GL_Ortho done\n";
}

/**
 * Set the current color for drawn shapes.
 */
void mglColor(MGLfloat red,
              MGLfloat green,
              MGLfloat blue)
{
	curr_color = vec3(red, green, blue);
}
