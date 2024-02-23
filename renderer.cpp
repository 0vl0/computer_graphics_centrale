#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
 
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include <stdio.h>
 
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <limits>
#include <algorithm>

#include <string>

const double PI = 3.1415926535897932; 

 
class Vector{
public:
    explicit Vector(double x=0., double y=0., double z=0.){
        coords[0] = x;
        coords[1] = y;
        coords[2] = z;
    };

    Vector& operator+=(const Vector& b){
        coords[0] += b[0];
        coords[1] += b[1];
        coords[2] += b[2];
        return *this;
    }

    Vector operator*=(const double x){
        coords[0] *= x;
        coords[1] *= x;
        coords[2] *= x;
        return *this;
    }

    Vector operator/=(const double x){
        coords[0] /= x;
        coords[1] /= x;
        coords[2] /= x;
        return *this;
    }

    double norm2(){
        /*
        Return ||v||_2^2, since ||v||_2 is rarely used in the code compared to the power of two of the norm.
        */
        return coords[0]*coords[0] + coords[1]*coords[1] + coords[2]*coords[2];
    }

    void normalize(){
        double n2 = sqrt(norm2());
        coords[0] /= n2;
        coords[1] /= n2;
        coords[2] /= n2;
    }

    int argmin() {
        int minIndex = 0;
        int minValue = abs(coords[0]);
        if (abs(coords[1]) < minValue){minIndex = 1;} 
        if (abs(coords[2]) < minValue){minIndex = 2;} 
        return minIndex;
    }

    const double& operator[](int i) const {return coords[i];}
    double& operator[](int i) {return coords[i];}

private:
    double coords[3];
};

Vector operator+(const Vector& a, const Vector& b){
    return Vector(a[0]+b[0], a[1]+b[1], a[2]+b[2]);
}

Vector operator-(const Vector& a, const Vector& b){
    return Vector(a[0]-b[0], a[1]-b[1], a[2]-b[2]);
}

// vectorial product
Vector cross(const Vector& a, const Vector& b){
    return Vector(a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]);
}

Vector operator*(const double x, const Vector& a){
    return Vector(x*a[0], x*a[1], x*a[2]);
}

Vector operator*(const Vector& a, const Vector& b){
    return Vector(a[0]*b[0], a[1]*b[1], a[2]*b[2]);
}


double dot(const Vector& a, const Vector& b){
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

double gammaCorrection(const double gamma, const double x){
    return std::min(255., pow(x, 1/gamma));
}

double getRandom(float maxi){
     return maxi*rand()/RAND_MAX;
}

class Ray{
public:
    explicit Ray(Vector origin, Vector direction, int d=5, bool fromDiffuse=false){
        O = origin;
        u = direction;
        depth = d;
        lastSphereDiffuse = fromDiffuse;
    }

    void normalize(){
        u.normalize();
    }

    Vector getDirection(){return u;}
    Vector getOrigin(){return O;}

    const Vector getDirection() const {return u;}
    const Vector getOrigin() const {return O;}

    const int getDepth() const {return depth;}

    bool fromDiffuseObject() const {
        return lastSphereDiffuse;
    }

private:
    Vector O; 
    Vector u;
    int depth;
    bool lastSphereDiffuse;
};

class Geometry {
public:
    Geometry(Vector color, bool reflection=false, bool refraction=false, float refraction_index=1, bool light_source=false){
        albedo = color;
        reflective = reflection;
        refractive = refraction;
        n = refraction_index;
        light = light_source;
    }

    virtual bool intersect(const Ray& ray, double& t, Vector& N, Vector& P) = 0;

    bool isReflective(){
        return reflective;
    }

    bool isRefractive(){
        return refractive;
    }

    float getN(){
        return n;
    }
       
    Vector getAlbedo(){
        return albedo;
    }

    bool isLightSource(){
        return light;
    }

private:
    Vector albedo;
    bool reflective;
    bool refractive;
    float n;
    bool light;
};


class Sphere : public Geometry  {
public:
    explicit Sphere(Vector center, double radius, Vector color, bool reflection=false, bool refraction=false, float refraction_index=1, bool source=false)
    : Geometry(color, reflection, refraction, refraction_index, source){
        C = center; 
        R = radius;
    }

    bool intersect(const Ray& ray, double& t, Vector& N, Vector& P){
        /*
        Return true if the ray intersects the sphere.
        The parametric equation of the Ray is O(t) = O+tu.
        The intersection double t is referenced, so that it can be accessed outside of the function.
        The function only returns true or false, depending if intersection occurs.
        */
        const Vector O = ray.getOrigin();
        const Vector u = ray.getDirection();

        // (b/2)^2
        double b = dot(u, O-C);
        double c = (O-C).norm2() - pow(R,2);

        // a = 0, because u is unitary
        double delta = pow(b,2) - c;

        if (delta < 0){return false;}

        double sqrt_delta = sqrt(delta);
        t = -b+sqrt_delta;
        if (t < 0){
            return false;
        }
        double t1 = -b-sqrt_delta;
        if (t1 > 0){
            t = t1;
        }

        // intersection point
        P = O+t*u;
        N = P-getCenter();
        return true;
    }


    Vector getCenter(){
        return C;
    }

    double getRadius(){
        return R;
    }


private:
    Vector C; 
    double R;
};

class BoundingBox{
public:
Vector m;
Vector M;

    BoundingBox(Vector mini, Vector maxi){
        m = mini;
        M = maxi;
    };

    BoundingBox(){};

    bool intersect(const Ray& R){
        // whether or not there is an intersection
        Vector inv_u(1./R.getDirection()[0], 1./R.getDirection()[1], 1./R.getDirection()[2]);
        double tmx = (m[0]-R.getOrigin()[0])*inv_u[0];
        double tMx = (M[0]-R.getOrigin()[0])*inv_u[0];
        double t1x = std::min(tmx, tMx);
        double t2x = std::max(tmx, tMx);

        double tmy = (m[1]-R.getOrigin()[1])*inv_u[1];
        double tMy = (M[1]-R.getOrigin()[1])*inv_u[1];
        double t1y = std::min(tmy, tMy);
        double t2y = std::max(tmy, tMy);

        double tmz = (m[2]-R.getOrigin()[2])*inv_u[2];
        double tMz = (M[2]-R.getOrigin()[2])*inv_u[2];
        double t1z = std::min(tmz, tMz);
        double t2z = std::max(tmz, tMz);

        double t_entry = std::max(std::max(t1x, t1y), t1z);
        double t_exit = std::min(std::min(t2x, t2y), t2z);

        if (t_exit < 0){return false;}
        if (t_exit > t_entry){return true;}

        return false;
    }

};

// BVH node
class Node{
public:
    Node* left;
    Node* right;
    BoundingBox bbox;
    int start;
    int end;

    Node(Node* l, Node* r, BoundingBox box, int index_start, int index_end){
        left = l;
        right = r;
        bbox = box;
        start = index_start;
        end = index_end;
    }

    Node(){}
};

class TriangleIndices {
public:
    TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {
    };
    int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
    int uvi, uvj, uvk;  // indices within the uv coordinates array
    int ni, nj, nk;  // indices within the normals array
    int group;       // face group
};

class BVH{
public:
Node* root;
std::vector<Vector> vector_vertices; // vertices
std::vector<int> vector_indices; // sorted indices
std::vector<TriangleIndices> vector_indices_triangles; // triangles

    BVH(std::vector<Vector> v_vertices){
        vector_vertices = v_vertices;
        Node* root = new Node();
        for (int i = 0; i<vector_vertices.size(); i++){
            vector_indices.push_back(i); 
        }
    }

    BVH(){}

    int get_index(int i){
        return vector_indices[i];
    }

    void init(std::vector<Vector> v_vertices, std::vector<TriangleIndices> indices){
        vector_vertices = v_vertices;
        vector_indices_triangles = indices;
        Node* root = new Node();
        for (int i = 0; i<vector_indices_triangles.size(); i++){
            vector_indices.push_back(i); 
        }
    }

    Vector get_barycenter(int start, int end){
        Vector barycenter(0., 0., 0.);
        for (int i=start; i<end; i++){
            barycenter += compute_barycenter(i);
            // vector_vertices[vector_indices[i]];
        }
        barycenter /= (end-start+1);
        return barycenter;
    }

    void construct(int start, int end){
        root = construct_bvh(start, end);
    }

    Vector compute_barycenter(int index_triangle){
        // compute barycenter of index_triangle-th triangle
        TriangleIndices indices_triangle_i = vector_indices_triangles[index_triangle];
        Vector b(0., 0., 0.);
        int index_1 = indices_triangle_i.vtxi;
        int index_2 = indices_triangle_i.vtxj;
        int index_3 = indices_triangle_i.vtxk;

        b += vector_vertices[index_1];
        b += vector_vertices[index_2];
        b += vector_vertices[index_3];

        b /= 3.0;

        return b;
    }

    Node* construct_bvh(int start, int end){
        // construct the BVH tree.
        // Vector barycenter = get_barycenter(start, end);
        BoundingBox bbox = compute_bbox(start, end);
        // stop if there are less than 5 triangles in the bounding box
        Node* node = new Node(NULL, NULL, bbox, start, end);

        Vector size_bbox = bbox.M-bbox.m;
        int dimension_to_split = 0;
        if (size_bbox[1] > size_bbox[0] && size_bbox[1] > size_bbox[2]){
            dimension_to_split = 1;
        } 
        if (size_bbox[2] > size_bbox[0] && size_bbox[2] > size_bbox[1]){
            dimension_to_split = 2;
        } 
        double threshold = bbox.m[dimension_to_split] + size_bbox[dimension_to_split]/2;

        int index_pivot = start;
        for (int i = start; i<end; i++){
            if (compute_barycenter(vector_indices[i])[dimension_to_split] < threshold){
                std::swap(vector_indices[i], vector_indices[index_pivot]);
                index_pivot += 1;
            }
        }
        
        if (end-start+1 <= 3 || index_pivot <= start || index_pivot >= end){
            return node;
        }
        node->left = construct_bvh(start, index_pivot-1); // index_pivot excluded
        node->right = construct_bvh(index_pivot, end); // end+1 excluded

        return node;
    }

    BoundingBox compute_bbox(int index_triangle_start, int index_triangle_end){
        BoundingBox bbox;
        bbox.m = Vector(1E9, 1E9, 1E9);
        bbox.M = Vector(-1E9, -1E9, -1E9);

        // for(int i = index_triangle_start; i<index_triangle_end; i++){
        //    for(int j = 0; j<3; j++){
        //     bbox.m[j] = std::min(bbox.m[j], vector_vertices[vector_indices[i]][j]);
        //     bbox.M[j] = std::max(bbox.M[j], vector_vertices[vector_indices[i]][j]);
        //    } 
        // }
        for(int i = index_triangle_start; i<index_triangle_end+1; i++){
            TriangleIndices indices_triangle_i = vector_indices_triangles[vector_indices[i]];
            int index_1 = indices_triangle_i.vtxi;
            int index_2 = indices_triangle_i.vtxj;
            int index_3 = indices_triangle_i.vtxk;
                
            for(int j = 0; j<3; j++){
                bbox.m[j] = std::min(bbox.m[j], vector_vertices[index_1][j]);
                bbox.m[j] = std::min(bbox.m[j], vector_vertices[index_2][j]);
                bbox.m[j] = std::min(bbox.m[j], vector_vertices[index_3][j]);

                bbox.M[j] = std::max(bbox.M[j], vector_vertices[index_1][j]);
                bbox.M[j] = std::max(bbox.M[j], vector_vertices[index_2][j]);
                bbox.M[j] = std::max(bbox.M[j], vector_vertices[index_3][j]);
            } 
        }

        return bbox;
    }

};
 
 
class TriangleMesh : public Geometry{
public:
  ~TriangleMesh() {}
  TriangleMesh(Vector color, bool reflection=false, bool refraction=false, float refraction_index=1, bool light=false) : Geometry(color, reflection, refraction, refraction_index, light){
    bbox = BoundingBox();
    bvh = BVH();
  };
    
    void readOBJ(const char* obj) {
 
        char matfile[255];
        char grp[255];
 
        FILE* f;
        f = fopen(obj, "r");
        int curGroup = -1;
        while (!feof(f)) {
            char line[255];
            if (!fgets(line, 255, f)) break;
 
            std::string linetrim(line);
            linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
            strcpy(line, linetrim.c_str());
 
            if (line[0] == 'u' && line[1] == 's') {
                sscanf(line, "usemtl %[^\n]\n", grp);
                curGroup++;
            }
 
            if (line[0] == 'v' && line[1] == ' ') {
                Vector vec;
 
                Vector col;
                if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6) {
                    col[0] = std::min(1., std::max(0., col[0]));
                    col[1] = std::min(1., std::max(0., col[1]));
                    col[2] = std::min(1., std::max(0., col[2]));
 
                    vertices.push_back(vec);
                    vertexcolors.push_back(col);
 
                } else {
                    sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                    vertices.push_back(vec);
                }
            }
            if (line[0] == 'v' && line[1] == 'n') {
                Vector vec;
                sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                normals.push_back(vec);
            }
            if (line[0] == 'v' && line[1] == 't') {
                Vector vec;
                sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
                uvs.push_back(vec);
            }
            if (line[0] == 'f') {
                TriangleIndices t;
                int i0, i1, i2, i3;
                int j0, j1, j2, j3;
                int k0, k1, k2, k3;
                int nn;
                t.group = curGroup;
 
                char* consumedline = line + 1;
                int offset;
 
                nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
                if (nn == 9) {
                    if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                    if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                    if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                    if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
                    if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
                    if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
                    if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
                    if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
                    if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
                    indices.push_back(t);
                } else {
                    nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
                    if (nn == 6) {
                        if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                        if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                        if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                        if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
                        if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
                        if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
                        indices.push_back(t);
                    } else {
                        nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
                        if (nn == 3) {
                            if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                            if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                            if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                            indices.push_back(t);
                        } else {
                            nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
                            if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                            if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                            if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                            if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
                            if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
                            if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
                            indices.push_back(t);
                        }
                    }
                }
 
                consumedline = consumedline + offset;
 
                while (true) {
                    if (consumedline[0] == '\n') break;
                    if (consumedline[0] == '\0') break;
                    nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
                    TriangleIndices t2;
                    t2.group = curGroup;
                    if (nn == 3) {
                        if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                        if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                        if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                        if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
                        if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
                        if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
                        if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
                        if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
                        if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;
                        indices.push_back(t2);
                        consumedline = consumedline + offset;
                        i2 = i3;
                        j2 = j3;
                        k2 = k3;
                    } else {
                        nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
                        if (nn == 2) {
                            if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                            if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                            if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                            if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
                            if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
                            if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
                            consumedline = consumedline + offset;
                            i2 = i3;
                            j2 = j3;
                            indices.push_back(t2);
                        } else {
                            nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
                            if (nn == 2) {
                                if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                                if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                                if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                                if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
                                if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
                                if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;                             
                                consumedline = consumedline + offset;
                                i2 = i3;
                                k2 = k3;
                                indices.push_back(t2);
                            } else {
                                nn = sscanf(consumedline, "%u%n", &i3, &offset);
                                if (nn == 1) {
                                    if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                                    if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                                    if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                                    consumedline = consumedline + offset;
                                    i2 = i3;
                                    indices.push_back(t2);
                                } else {
                                    consumedline = consumedline + 1;
                                }
                            }
                        }
                    }
                }
 
            }
 
        }
        fclose(f);
        // init_bounding_box();
    }
 
    // bool intersect(const Ray& ray, double& t, Vector& N, Vector& P){
    //     t = 1E19;
    //     if (bbox.intersect(ray)){
    //         for(int i =0; i<indices.size(); i++){
    //             TriangleIndices indices_triangle_i = indices[i];
    //             int index_1 = indices_triangle_i.vtxi;
    //             int index_2 = indices_triangle_i.vtxj;
    //             int index_3 = indices_triangle_i.vtxk;

    //             Vector A = vertices[index_1];
    //             Vector B = vertices[index_2];
    //             Vector C = vertices[index_3];
    //             Vector e1 = B-A;
    //             Vector e2 = C-A;
    //             N = cross(e1, e2);

    //             Vector u = ray.getDirection();
    //             Vector O = ray.getOrigin();

    //             double beta = dot(e2, cross(A-O, u))/dot(u, N);
    //             double gamma = -dot(e1, cross(A-O, u))/dot(u, N);
    //             double alpha = 1-beta-gamma;
    //             t = dot(A-O, N)/dot(u,N);
    //             P = A + beta*e1 + gamma*e2;
    //             if (t >= 0 && beta >= 0 && gamma >= 0 && alpha >= 0 && beta <= 1 && gamma <= 1 && alpha <= 1){return true;}
    //         }
    //     }
    //     return false;
    // }

    void init_bounding_box(){
        bbox.m = Vector(1E9, 1E9, 1E9);
        bbox.M = Vector(-1E9, -1E9, -1E9);

        for(int i = 0; i<vertices.size(); i++){
           for(int j = 0; j<3; j++){
                bbox.m[j] = std::min(bbox.m[j], vertices[i][j]);
                bbox.M[j] = std::max(bbox.M[j], vertices[i][j]);
           } 
        }
    }

    void translation(Vector translation_vector, double scale_factor){
        for(int i = 0; i<vertices.size(); i++){
            vertices[i] += translation_vector;
            vertices[i] *= scale_factor;
        }
        // init_bounding_box();
        bvh.init(vertices, indices);
        // int a = indices.size();
        bvh.construct(0, indices.size()-1);
    }

    bool intersect(const Ray& ray, double& t, Vector& N, Vector& P){
        std::vector<Node*> list_nodes;
        list_nodes.push_back(bvh.root);
        bool intersection = false;
        Vector N_min(0., 0., 0.);
        Vector P_min(0., 0., 0.);
        double t_min = 1E19;
        t = 1E19;
        while (!list_nodes.empty()){
            Node* node = list_nodes.back();
            list_nodes.pop_back();
            if (node->right != nullptr){
                if (node->right->bbox.intersect(ray)){
                    list_nodes.push_back(node->right);
                }
            }
            if (node->left != nullptr){
                if (node->left->bbox.intersect(ray)){
                    list_nodes.push_back(node->left);
                }
            }
            // check intersection with ray
            // leaf, excluding node.left in leaves
            if (node->right == nullptr && node->left == nullptr){
                // bbox = node.bbox;
                for(int i = node->start; i<node->end+1; i++){
                    TriangleIndices indices_triangle_i = indices[bvh.get_index(i)];
                    int index_1 = indices_triangle_i.vtxi;
                    int index_2 = indices_triangle_i.vtxj;
                    int index_3 = indices_triangle_i.vtxk;

                    int index_normal_1 = indices_triangle_i.ni;
                    int index_normal_2 = indices_triangle_i.nj;
                    int index_normal_3 = indices_triangle_i.nk;

                    Vector N_A = normals[index_normal_1];
                    Vector N_B = normals[index_normal_2];
                    Vector N_C = normals[index_normal_3];


                    Vector A = vertices[index_1];
                    Vector B = vertices[index_2];
                    Vector C = vertices[index_3];
                    Vector e1 = B-A;
                    Vector e2 = C-A;
                    N = cross(e1, e2);

                    Vector u = ray.getDirection();
                    Vector O = ray.getOrigin();

                    double beta = dot(e2, cross(A-O, u))/dot(u, N);
                    double gamma = -dot(e1, cross(A-O, u))/dot(u, N);
                    double alpha = 1-beta-gamma;
                    t = dot(A-O, N)/dot(u,N);
                    P = A + beta*e1 + gamma*e2;

                    Vector shading_normal = (alpha*N_A + beta*N_B + gamma*N_C);
                    shading_normal.normalize();
                    N = shading_normal;

                    if (t >= 0 && beta >= 0 && gamma >= 0 && alpha >= 0 && beta <= 1 && gamma <= 1 && alpha <= 1){
                        intersection = true;
                        if (t < t_min){
                            t_min = t;
                            N_min = N;
                            P_min = P;
                        }
                    }
                }
            }
        }
        N = N_min;
        t = t_min;
        P = P_min;
        return intersection;
    }

    std::vector<TriangleIndices> indices;
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    std::vector<Vector> uvs;
    std::vector<Vector> vertexcolors;

    BoundingBox bbox;

    BVH bvh;
};




class Scene{
public:
    Scene(Vector source, double I, float refraction_index = 1.0, int path_per_pixel = 32, int index_light=0){
        lightSource = source;
        lightIntensity = I;
        list_objects = {};
        n = refraction_index;
        rayPerPixel = path_per_pixel;
        index_light_sphere = index_light;
    }

    int intersect(const Ray& ray, double& t, Vector& N, Vector &P){
        double t_min = std::numeric_limits<int>::max(); 
        int i_min = -1;
        Vector N_min(0., 0., 0.);
        Vector P_min(0., 0., 0.);
        for (int i=0; i<list_objects.size(); i++){
            if (list_objects[i]->intersect(ray, t, N, P)){
                if (t < t_min){
                    t_min = t;
                    i_min = i;
                    N_min = N;
                    P_min = P;
                }
            }
        }
        t = t_min;
        N = N_min;
        P = P_min;
        return i_min;
    }

    Vector get_random_cos(Vector N){
        // random_cos
        double r1 = getRandom(1);
        double r2 = getRandom(1);

        float x = cos(2*PI*r1)*sqrt(1-r2);
        float y = sin(2*PI*r1)*sqrt(1-r2);
        float z = sqrt(r2);

        // random cos??
        // <w,N>/PI

        int minIndex = N.argmin();
        Vector T1, T2;
        if (minIndex == 0){
            T1 = Vector(0, -N[2], N[1]);
        }
        else if (minIndex == 1){
            T1 = Vector(N[2], 0, -N[0]);
        }
        else if (minIndex == 2){
            T1 = Vector(-N[1], N[0], 0);
        }

        // P': random cos
        // omega_i(P'): de P a P', choisi aleatoirement
        // N': normale au point P'
        // normale au point P': vecteur u, vecteur generee = N'

        T1.normalize();
        T2 = cross(N,T1);
        Vector random_direction = z*N + x*T1 + y*T2;
        return random_direction;
    }

    Vector getColor(const Ray& R){
        Vector color(0., 0., 0.);
        if (R.getDepth() < 0){return color;}

        double t = 0.;
        Vector N(0., 0., 0.);
        Vector P(0., 0., 0.);
        int i_min = intersect(R, t, N, P);
        if (i_min >= 0){
            // intersection point
            // Vector P = R.getOrigin() + t*R.getDirection();

            // intersection object closest to camera
            Geometry* sphereIntersected = getObject(i_min);
            //Vector C = sphereIntersected->getCenter();

            // normal vector at the intersection point
            //Vector N = P-C;
            N.normalize();

            float offset = 0.001;

            if (sphereIntersected->isReflective()){
                /* --- Reflective surfaces ---
                   Reflected direction is just incident direction with the normal direction reversed.
                   (think of a ball hitting a wall) 
                */
                Ray reflectedRay(P+offset*N, R.getDirection() - 2*dot(R.getDirection(),N)*N, R.getDepth()-1);
                reflectedRay.normalize();
                color = getColor(reflectedRay);
            }else if(sphereIntersected->isRefractive()){
                /* --- Refractive surfaces ---
                Fresnel law is used to compute the refraction R and transmission T coefficient.
                Instead of computing both the reflected and transmitted light for each ray, a random number
                is chosen and only one ray is sent.
                We get the same resulting light by averaging the rays.
                */

                float n1, n2;
                bool exiting =  dot(R.getDirection(), N) > 0;
                offset *= -1;
                if (exiting){
                    N = (-1)*N;
                    n1 = sphereIntersected->getN();
                    n2 = n;
                } 
                else{
                    n2 = sphereIntersected->getN();
                    n1 = n;
                }
                
                double k0 = pow((n1-n2),2)/pow((n1+n2),2);
                double reflection_coef = k0 + (1-k0)*pow((1-abs(dot(N,R.getDirection()))),5); 
                double rnd_number = rand() / double(RAND_MAX);
                //rnd_number = 0;
                if (rnd_number < reflection_coef){
                    // reflection
                    Ray reflectedRay(P-offset*N, R.getDirection() - 2*dot(R.getDirection(),N)*N, R.getDepth()-1);
                    reflectedRay.normalize();
                    color = getColor(reflectedRay);     
                }
                else{
                    // refraction 
                    float R_normal_scalar = 1-pow((n1/n2),2)*(1-pow(dot(R.getDirection(), N),2));
                    // Check if total reflection occurs
                    // If that's the case, black color is returned
                    if (R_normal_scalar < 0){
                        color = Vector(0.,0.,0.);
                    }
                    else{
                        // Check if Ray is entering or exiting the object, and change offset / normal direction accordingly
                        Vector R_normal = -sqrt(R_normal_scalar)*N;
                        Vector R_tangential = (n1/n2)*(R.getDirection()-dot(R.getDirection(),N)*N); 
                        Ray refractedRay(P+offset*N, R_normal+R_tangential, R.getDepth()-1);
                        // New direction is R_normal+R_tangential and is already normalized
                        // It can be seens by computing the equation for the norm of refractedRay.
                        refractedRay.normalize();
                        color = getColor(refractedRay);
                    }
                }
                
            }else{
                /* --- Diffuse surfaces ---
                 Add a visiblity term, depending on position of the light and the intersection point.
                 It creates shadows if light is not visible from intersection point.
                */
                if (sphereIntersected->isLightSource()){
                    // double a = 1/pow(2*PI*sphereIntersected.getRadius(),2);
                    if (R.fromDiffuseObject()){
                        color = Vector(0., 0., 0.);
                    }
                    else{
                        color = (getLightIntensity()/pow(2*PI*getLightSphere()->getRadius(),2))*Vector(1., 1., 1.);
                    }
                    //color = Vector(1., 1., 1.);
                    // int a = 1;
                }
                else{
                    // Direct contribution
                    // Vector lightVector = getLightSource()-P;
                    // double d2 = lightVector.norm2();
                    // lightVector.normalize();

                    // Ray PLight(P+offset*N, lightVector);
                    // // Add normal light if light source is visible from P
                    // int i_min_shadow = intersect(PLight, t);
                    // if (i_min_shadow == -1 || t*t > d2){
                    //     color = (getLightIntensity()/(4*pow(PI,2)*d2))*std::max(0.,dot(N,lightVector))*sphereIntersected.getAlbedo();
                    // }

                    // D in the pdf document
                    Vector lightVector = getLightSphere()->getCenter()-P;
                    double d2 = lightVector.norm2();
                    lightVector.normalize();

                    // V in pdf
                    Vector random_direction = get_random_cos((-1)*lightVector);

                    // sample point on sphere, x' in pdf
                    Vector P_prime = getLightSphere()->getCenter() + getLightSphere()->getRadius()*random_direction;

                    Vector omega_i = P_prime-P; 
                    d2 = omega_i.norm2();
                    omega_i.normalize();
                    
                    // normal at the surface of light source sphere
                    Vector N_prime = P_prime-getLightSphere()->getCenter();
                    N_prime.normalize();

                    // proba to sample point on the light sphere
                    //double densite_proba = dot(N_prime, (-1)*lightVector)/(PI*pow(getLightSphere().getRadius(),2));
                    double densite_proba = dot(omega_i, lightVector)/(PI*pow(getLightSphere()->getRadius(),2));
                    

                    double L_wi = getLightIntensity()/(4.*pow(PI*getLightSphere()->getRadius(),2));

                    double dot1 = dot(N,omega_i);
                    double dot2 = dot(N_prime,(-1)*omega_i);

                    color =  L_wi*((std::max(0., dot(N,omega_i))*std::max(0., dot(N_prime,(-1)*omega_i)))/(d2*densite_proba))*sphereIntersected->getAlbedo();
                    // if (color[0] > 0 && color[1] > 0 && color[2] > 0 ){
                    //     int a = 1;
                    // }

                    Vector omegaIndirect = get_random_cos(N);
                    Ray indirectRay = Ray(P+offset*N, omegaIndirect, R.getDepth()-1, true);
                    Vector colorIndirect = sphereIntersected->getAlbedo()*getColor(indirectRay);
                    color = color + colorIndirect;
                }
            }
        }
        return color;
    }

    Geometry* getObject(int index){
        return list_objects[index];
    }

    void add(Geometry* S){
        list_objects.push_back(S);
    }

    Vector getLightSource(){
        return lightSource;
    }

    double getLightIntensity(){
        return lightIntensity;
    }

    int getRayPerPixel(){
        return rayPerPixel;
    }

    int get_index_light_sphere(){
        return index_light_sphere;
    }

    // cast to Sphere
    Sphere* getLightSphere(){
        return dynamic_cast<Sphere*>(list_objects[index_light_sphere]);
    }


private:
    std::vector<Geometry*> list_objects; 
    Vector lightSource;
    double lightIntensity;
    float n;
    int rayPerPixel;
    int index_light_sphere;
};

int main() {
    int W = 128;
    int H = 128;
    const double gamma = 2.2;

    TriangleMesh* m = new TriangleMesh(Vector(1., 1., 1.), false, false);
    m->readOBJ("cadnav.com_model/Models_F0202A090/cat.obj");
    m->translation(Vector(0., -30., 0.), 0.6);

    // horizontal field of view
    const double fov = 60*(PI/180);
    const Vector camera(0,0,55);
    
    Scene scene(Vector(-10, 20, 40), 1e10);
    Sphere* lightSource = new Sphere(Vector(-10, 25, 40), 5, Vector(1.,1.,1.), true, false, 1.5, true);
    Sphere* centralSphere = new Sphere(Vector(0., 0., 0.), 10, Vector(1.,1.,1.), false, true, 1.5);
    Sphere* centralSphereLeft = new Sphere(Vector(22., 0., 0.), 10, Vector(1.,1.,1.), false, false, 1.5);
    Sphere* centralSphereRight = new Sphere(Vector(-22., 0., 0.), 10, Vector(1.,1.,1.), true, false, 1.5);
    Sphere* leftSphere = new Sphere(Vector(0., 0., 1000.), 940, Vector(1., 0., 1.));
    Sphere* rearSphere = new Sphere(Vector(0., -1000., 0.), 990, Vector(0., 0., 1.));
    Sphere* rightSphere = new Sphere(Vector(0.,0.,-1000.), 940, Vector(0., 1., 0.));
    Sphere* topSphere = new Sphere(Vector(0., 1000., 0.), 940, Vector(1., 0., 0.));
//
    scene.add(lightSource);
    // scene.add(centralSphere);
    // scene.add(centralSphereLeft);
    // scene.add(centralSphereRight);
    
    scene.add(leftSphere);
    scene.add(rearSphere);
    scene.add(rightSphere);
    scene.add(topSphere);

    scene.add(m);
    
    std::vector<unsigned char> image(W*H*3, 0);
    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            // ray vector

            double red = 0; 
            double green = 0;
            double blue = 0;
            #pragma omp parallel for num_threads(6)
            for (int k=0; k<scene.getRayPerPixel(); k++){
                //anti aliasing
                double r1 = getRandom(1);
                double r2 = getRandom(1);
                double g1 = log(-2*log(r1))*cos(2*PI*r2)*0.25;
                double g2 = log(-2*log(r1))*sin(2*PI*r2)*0.25;
                Ray R(camera, Vector(camera[0]+j-W/2+0.5+g1, camera[1]+H/2-i-0.5+g2, camera[2]-W/(2*tan(fov/2))));
                R.normalize();

                Vector color = scene.getColor(R);
                red += color[0];
                green += color[1];
                blue += color[2];
            }
            red = gammaCorrection(gamma, (red/scene.getRayPerPixel()));
            green = gammaCorrection(gamma, (green/scene.getRayPerPixel()));
            blue = gammaCorrection(gamma, (blue/scene.getRayPerPixel()));

            image[(i*W + j) * 3 + 0] = red;
            image[(i*W + j) * 3 + 1] = green;
            image[(i*W + j) * 3 + 2] = blue;
        }
    }
    stbi_write_png("test_shading_normal.png", W, H, 3, &image[0], 0);
 
    return 0;
}