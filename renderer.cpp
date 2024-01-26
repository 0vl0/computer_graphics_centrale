#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
 
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
 
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

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

double dot(const Vector& a, const Vector& b){
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}


class Ray{
public:
    explicit Ray(Vector origin, Vector direction){
        O = origin;
        u = direction;
    }

    void normalize(){
        u.normalize();
    }

    Vector getDirection(){return u;}
    Vector getOrigin(){return O;}

    const Vector getDirection() const {return u;}
    const Vector getOrigin() const {return O;}

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


    bool intersect(const Ray& ray, double& t){
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


    // horizontal field of view
    const double fov = 60*(PI/180);
    const Vector camera(0,0,55);
    Vector light(-10, 20, 40);
    //light.normalize();
    
    Scene scene;
    Sphere centralSphere(Vector(0.,0.,0.), 10, Vector(255.,255.,255.));
    scene.add(centralSphere);
    
    std::vector<unsigned char> image(W*H*3, 0);
    double t = 0.;
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            // ray vector
            Ray R(camera, Vector(camera[0]+j-W/2+0.5, camera[1]+H/2-i-0.5, camera[2]-W/(2*tan(fov/2))));
            R.normalize();

            if (centralSphere.intersect(R,t)){
                image[(i*W + j) * 3 + 0] = 255;
                image[(i*W + j) * 3 + 1] = 255;
                image[(i*W + j) * 3 + 2] = 255;
            }
 
        }
    }
    stbi_write_png("image_sphere.png", W, H, 3, &image[0], 0);
 
    return 0;
}