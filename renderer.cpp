#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
 
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
 
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
 
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

    const double& operator[](int i) const {return coords[i];}
    double& operator[](int i) {return coords[i];}

private:
    double coords[3];
};

Vector operator+(const Vector& a, const Vector& b){
    return Vector(a[0]+b[0], a[1]+b[1], a[2]+b[2]);
}

double dot(const Vector& a, const Vector& b){
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}


class Ray{
public:
    explicit Ray(Vector origin, Vector direction){
        O = origin;
        u = direction;
    }

private:
    Vector O; 
    Vector u;
};


class Sphere{
public:
    explicit Sphere(Vector center, double radius, Vector color){
        C = center; 
        R = radius;
        albedo = color;
    }


    bool intersect(const Ray& ray){
    /*
    Returns true if the ray intersects the sphere.
    Important parameters are returned with references.
    */
        return true;
    }

private:
    Vector C; 
    double R;
    Vector albedo;
};

class Scene{
public:
    Scene(){
        list_objects = {};
    }

    bool intersect(const Ray& ray){
        return true;
    }

    void add(Sphere S){
        list_objects.push_back(S);
    }

private:
    std::vector<Sphere> list_objects; 
};

 
int main() {
    int W = 512;
    int H = 512;

    Vector camera(0,0,55);
    Vector light(-10, 20, 40);
    
    Scene scene;
    Sphere centralSphere(Vector(0.,0.,0.), 10.0, Vector(255.,255.,255.));
    scene.add(centralSphere);
    
    std::vector<unsigned char> image(W*H * 3, 0);
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
 
            image[(i*W + j) * 3 + 0] = 255;
            image[(i*W + j) * 3 + 1] = 0;
            image[(i*W + j) * 3 + 2] = 0;
        }
    }
    stbi_write_png("image.png", W, H, 3, &image[0], 0);
 
    return 0;
}