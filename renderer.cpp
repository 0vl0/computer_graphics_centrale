#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
 
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
 
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <limits>

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

Vector operator*(const double x, const Vector& a){
    return Vector(x*a[0], x*a[1], x*a[2]);
}

double dot(const Vector& a, const Vector& b){
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

double gammaCorrection(const double gamma, const double x){
    return std::min(255., pow(x, 1/gamma));
}


class Ray{
public:
    explicit Ray(Vector origin, Vector direction, int d=5){
        O = origin;
        u = direction;
        depth = d;
    }

    void normalize(){
        u.normalize();
    }

    Vector getDirection(){return u;}
    Vector getOrigin(){return O;}

    const Vector getDirection() const {return u;}
    const Vector getOrigin() const {return O;}

    const int getDepth() const {return depth;}

private:
    Vector O; 
    Vector u;
    int depth;
};


class Sphere{
public:
    explicit Sphere(Vector center, double radius, Vector color, bool reflection=false, bool refraction=false, float refraction_index=1.){
        C = center; 
        R = radius;
        albedo = color;
        reflective = reflection;
        refractive = refraction;
        n = refraction_index;
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

    Vector getAlbedo(){
        return albedo;
    }

    Vector getCenter(){
        return C;
    }

    bool isReflective(){
        return reflective;
    }

    bool isRefractive(){
        return refractive;
    }

    float getN(){
        return n;
    }


private:
    Vector C; 
    double R;
    Vector albedo;
    bool reflective;
    bool refractive;
    float n;
};

class Scene{
public:
    Scene(Vector source, double I, float refraction_index = 1.0){
        lightSource = source;
        lightIntensity = I;
        list_objects = {};
        n = refraction_index;
    }

    int intersect(const Ray& ray, double& t){
        double t_min = std::numeric_limits<int>::max(); 
        int i_min = -1;
        for (int i=0; i<list_objects.size(); i++){
            if (list_objects[i].intersect(ray, t)){
                if (t < t_min){
                    t_min = t;
                    i_min = i;
                }
            }
        }
        t = t_min;
        return i_min;
    }

    Vector getColor(const Ray& R){
        Vector color(0., 0., 0.);
        if (R.getDepth() < 0){return color;}

        double t = 0.;
        int i_min = intersect(R, t);
        if (i_min >= 0){
            // intersection point
            Vector P = R.getOrigin() + t*R.getDirection();

            // intersection object closest to camera
            Sphere sphereIntersected = getObject(i_min);
            Vector C = sphereIntersected.getCenter();

            // normal vector at the intersection point
            Vector N = P-C;
            N.normalize();

            float offset = 0.001;

            if (sphereIntersected.isReflective()){
                /* --- Reflective surfaces ---
                   Reflected direction is just incident direction with the normal direction reversed.
                   (think of a ball hitting a wall) 
                */
                Ray reflectedRay(P+offset*N, R.getDirection() - 2*dot(R.getDirection(),N)*N, R.getDepth()-1);
                reflectedRay.normalize();
                color = getColor(reflectedRay);
            }else if(sphereIntersected.isRefractive()){
                /* --- Refractive surfaces ---
                */
                float n1, n2;
                bool exiting =  dot(R.getDirection(), N) > 0;
                offset *= -1;
                if (exiting){
                    N = (-1)*N;
                    n1 = sphereIntersected.getN();
                    n2 = n;
                } 
                else{
                    n2 = sphereIntersected.getN();
                    n1 = n;
                }
                
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
            }else{
                /* --- Diffuse surfaces ---
                 Add a visiblity term, depending on position of the light and the intersection point.
                 It creates shadows if light is not visible from intersection point.
                */
                Vector lightVector = getLightSource()-P;
                double d2 = lightVector.norm2();
                lightVector.normalize();

                Ray PLight(P+offset*N, lightVector);
                // Add normal light if light source is visible from P
                int i_min_shadow = intersect(PLight, t);
                if (i_min_shadow == -1 || t*t > d2){
                    color = (getLightIntensity()/(4*pow(PI,2)*d2))*std::max(0.,dot(N,lightVector))*sphereIntersected.getAlbedo();
                }
            }
        }
        return color;
    }

    Sphere getObject(int index){
        return list_objects[index];
    }

    void add(Sphere S){
        list_objects.push_back(S);
    }

    Vector getLightSource(){
        return lightSource;
    }

    double getLightIntensity(){
        return lightIntensity;
    }

private:
    std::vector<Sphere> list_objects; 
    Vector lightSource;
    double lightIntensity;
    float n;
};


int main() {
    int W = 512;
    int H = 512;
    const double gamma = 2.2;

    // horizontal field of view
    const double fov = 60*(PI/180);
    const Vector camera(0,0,55);
    //light.normalize();
    
    Scene scene(Vector(-10, 20, 40), 1e8);
    Sphere centralSphere(Vector(0., 0., 0.), 10, Vector(255.,255.,255.), false, true, 1.5);
    Sphere centralSphereLeft(Vector(22., 0., 0.), 10, Vector(255.,255.,255.), false, false, 1.5);
    Sphere centralSphereRight(Vector(-22., 0., 0.), 10, Vector(255.,255.,255.), true, false, 1.5);
    Sphere leftSphere(Vector(0., 0., 1000.), 940, Vector(255., 0., 255.));
    Sphere rearSphere(Vector(0., -1000., 0.), 990, Vector(0., 0., 255.));
    Sphere rightSphere(Vector(0.,0.,-1000.), 940, Vector(0., 255., 0.));
    Sphere topSphere(Vector(0., 1000., 0.), 940, Vector(255., 0., 0.));

    scene.add(centralSphere);
    scene.add(centralSphereLeft);
    scene.add(centralSphereRight);
    scene.add(leftSphere);
    scene.add(rearSphere);
    scene.add(rightSphere);
    scene.add(topSphere);
    
    std::vector<unsigned char> image(W*H*3, 0);
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            // ray vector
            Ray R(camera, Vector(camera[0]+j-W/2+0.5, camera[1]+H/2-i-0.5, camera[2]-W/(2*tan(fov/2))));
            R.normalize();

            Vector color = scene.getColor(R);
            image[(i*W + j) * 3 + 0] = gammaCorrection(gamma, color[0]);
            image[(i*W + j) * 3 + 1] = gammaCorrection(gamma, color[1]);
            image[(i*W + j) * 3 + 2] = gammaCorrection(gamma, color[2]);    
        }
    }
    stbi_write_png("image_v4_refraction.png", W, H, 3, &image[0], 0);
 
    return 0;
}