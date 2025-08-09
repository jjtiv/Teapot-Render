#include <iostream>
#include <fstream>
#include <strstream>
#include <math.h>
#include <cmath>
#include <C:\Toolbox\Eigen3\Eigen\Dense>
//#include <Eigen/Dense>
#define MAX std::numeric_limits<double>::max()

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


using namespace std;

bool EC = false;
double shadowBias = 1e-6;

double det(const Eigen::Vector3d &a, const Eigen::Vector3d &b, const Eigen::Vector3d &c) {
    return a[0]* (b[1] * c[2] - c[1] * b[2]) +
      b[0] * (c[1] * a[2] - a[1] * c[2]) +
      c[0] * (a[1] * b[2] - b[1] * a[2]);
}

struct Fill{
    Eigen::Vector3d rgb;
    double kd, ks, e, kt, ir;
};

class Ray{
    public:
        Eigen::Vector3d e;
        Eigen::Vector3d d;
        int depth;
        Ray(const Eigen::Vector3d &_e, const Eigen::Vector3d &_d, int _depth = 0) : e(_e), d(_d), depth(_depth) {};
};

class HitRecord {
    public:
      double t, alpha, beta, gamma;
      Eigen::Vector3d p, n, v;
      Fill f;
      int raydepth;
};

class Surface{
    public:
        virtual bool intersect(const Ray &r, double t0, double t1, HitRecord &hr) const = 0;
};

class Triangle : public Surface{
    public:
        Eigen::Vector3d a,b,c;
        Triangle(Eigen::Vector3d ai, Eigen::Vector3d bi, Eigen::Vector3d ci){
            a = ai;
            b = bi;
            c = ci;
        }
        bool intersect(const Ray &r, double t0, double t1, HitRecord &hr) const {
            Eigen::Vector3d ba = a-b;
            Eigen::Vector3d ca = a-c;
            Eigen::Vector3d ea = a-r.e;
            double detA = det(ba, ca, r.d);
            double beta = det(ea, ca, r.d)/detA;
            if (beta < 0 || beta > 1) return false;
            double gamma = det(ba, ea, r.d)/detA;
            if (gamma < 0.0 || gamma > 1.0-beta) return false;
            double t = det(ba, ca, ea)/detA;
            if (t < t0 || t > t1) return false;
            hr.t = t;
            hr.p = r.e + t * r.d;
            hr.n = ba.cross(ca);
            hr.n.normalize();
            hr.alpha = 1.0 - beta - gamma;
            hr.beta = beta;
            hr.gamma = gamma;
            return true;
        };
};

class TrianglePatch : public Triangle{
    public:
        Eigen::Vector3d n1,n2,n3;
        TrianglePatch(Eigen::Vector3d ai, Eigen::Vector3d bi, Eigen::Vector3d ci, Eigen::Vector3d n1i, Eigen::Vector3d n2i, Eigen::Vector3d n3i) : Triangle(ai, bi, ci), n1(n1i), n2(n2i), n3(n3i) {}
        bool intersect(const Ray &r, double t0, double t1, HitRecord &hr) const override {
            if(!EC){
                Eigen::Vector3d ba = a-b;
                Eigen::Vector3d ca = a-c;
                Eigen::Vector3d ea = a-r.e;
                double detA = det(ba, ca, r.d);
                double beta = det(ea, ca, r.d)/detA;
                if (beta < 0 || beta > 1) return false;
                double gamma = det(ba, ea, r.d)/detA;
                if (gamma < 0.0 || gamma > 1.0-beta) return false;
                double t = det(ba, ca, ea)/detA;
                if (t < t0 || t > t1) return false;
                hr.t = t;
                hr.p = r.e + t * r.d;
                hr.n = ba.cross(ca);
                hr.n.normalize();
                hr.alpha = 1.0 - beta - gamma;
                hr.beta = beta;
                hr.gamma = gamma;
                return true;
            }
        
            Eigen::Vector3d ba = a - b;
            Eigen::Vector3d ca = a - c;
            Eigen::Vector3d ea = a - r.e;
        
            double detA = det(ba, ca, r.d);
            if (fabs(detA) < 1e-8) return false;  // Avoid division by zero (parallel ray case)
        
            double beta = det(ea, ca, r.d) / detA;
            if (beta < 0 || beta > 1) return false;
        
            double gamma = det(ba, ea, r.d) / detA;
            if (gamma < 0 || gamma > 1 - beta) return false;
        
            double t = det(ba, ca, ea) / detA;
            if (t < t0 || t > t1) return false;

            // std::cout << "Intersection Found!\n";
            // std::cout << "Ray Origin: " << r.e.transpose() << "\n";
            // std::cout << "Ray Direction: " << r.d.transpose() << "\n";
            // std::cout << "Triangle Vertices: " << a.transpose() << " | " << b.transpose() << " | " << c.transpose() << "\n";
            // std::cout << "Barycentric Coords: alpha=" << (1.0 - beta - gamma) 
            //         << ", beta=" << beta << ", gamma=" << gamma << "\n";
        
            // Store intersection information
            hr.t = t;
            hr.p = r.e + t * r.d;
        
            // Interpolate normal using barycentric coordinates
            hr.n = (1.0 - beta - gamma) * n1 + beta * n2 + gamma * n3;
            hr.n.normalize();
        
            hr.alpha = 1.0 - beta - gamma;
            hr.beta = beta;
            hr.gamma = gamma;
        
            return true;
        }
        
};

// class Poly : public Surface{

// }

// class Polypatch : public Poly{

// }




struct Light{
    Eigen::Vector3d p, c;
};



struct Tracer{
    Eigen::Vector3d bcolor;
    Eigen::Vector3d from, at, up;
    double angle, hither;
    Eigen::Vector2d resolution;
    vector<Light> lights;
    vector<pair<Surface *, Fill> > surfaces;
    bool ec;

    Eigen::Vector3d* im;
};



void parseFile(Tracer& tracer, char* inputFile, char* outputFile);
void createImg(Tracer& tracer);
void writeImg(Tracer& tracer, char* outputFile);
Eigen::Vector3d castRay(Tracer& tracer, const Ray &r, double t0, double t1, int depth);

void parseFile(Tracer& tracer, char* inputFile, char* outputFile){
    ifstream file(inputFile);
    if(!file.is_open()){
        cout << "Unable to open file." << endl;
        return;
    }


    string line;
    Fill currFill;
    string ch;
    while(getline(file, line)){
        switch(line[0]){
            case 'b':{
                stringstream ss(line);
                string junk;
                ss >> junk >> tracer.bcolor[0] >> tracer.bcolor[1] >> tracer.bcolor[2];
            
                break;
            }
            case 'v':{
                getline(file, line);
                string junk;
                stringstream fromss(line);
                fromss >> junk >> tracer.from[0] >> tracer.from[1] >> tracer.from[2];
      
                getline(file, line);
                stringstream atss(line);
                atss >> junk >> tracer.at[0] >> tracer.at[1] >> tracer.at[2];
        
                getline(file, line);
                stringstream upss(line);
                upss >> junk >> tracer.up[0] >> tracer.up[1] >> tracer.up[2];
        
                getline(file, line);
                stringstream angless(line);
                angless >> junk >> tracer.angle;
        
                getline(file, line);
                stringstream hitherss(line);
                hitherss >> junk >> tracer.hither;
        
                getline(file, line);
                stringstream resolutionss(line);
                resolutionss >> junk >> tracer.resolution[0] >> tracer.resolution[1];


                break;
            }
            case 'f':{
                stringstream ss(line);
                string junk;
                ss >> junk >> currFill.rgb[0] >> currFill.rgb[1] >> currFill.rgb[2] >> currFill.kd >> currFill.ks >> currFill.e >> currFill.kt >> currFill.ir;
                break;
            }

            case 'p': {
                bool patch = false;
                stringstream ssn(line);
                unsigned int nverts;
                
                if (line[1] == 'p') {
                    patch = true; 
                    //ssn >> ch; 
                }
                
                ssn >> ch >> nverts;
            
                vector<Eigen::Vector3d> vertices;
                vector<Eigen::Vector3d> normals;
                
                for (unsigned int i = 0; i < nverts; i++) {
                    getline(file, line);
                    stringstream ss(line);
                    Eigen::Vector3d v, n;
                    
                    if (patch) {
                        ss >> v[0] >> v[1] >> v[2] >> n[0] >> n[1] >> n[2];
                    } else {
                        ss >> v[0] >> v[1] >> v[2]; 
                    }
                    
                    vertices.push_back(v);
                    if (patch) {
                        n.normalize();  
                        normals.push_back(n);
                    }
                }
            
                // Debugging: Print the line and patch status
                // if(vertices.size() == 3){
                //     cout << "Line read: " << line << endl;
                //     cout << "Patch flag: " << patch << " Number of vertices: " << nverts << endl;
                // }
            
                // Handling the pp 3 case
                if (vertices.size() == 3) {
                    if (patch) {
                        // Debugging: Ensure this section is triggered for PP 3
                        // cout << "Adding PP 3 patch: "
                        //      << vertices[0].transpose() << ", "
                        //      << vertices[1].transpose() << ", "
                        //      << vertices[2].transpose() << endl;
                        
                        tracer.surfaces.push_back(pair<Surface *, Fill>(
                            new TrianglePatch(vertices[0], vertices[1], vertices[2], 
                                              normals[0], normals[1], normals[2]), currFill));
                    }
                }
                
                // Handling the p 4 case
                else if (vertices.size() == 4) {
                    if (!patch) {
                        Eigen::Vector3d n0 = (vertices[1] - vertices[0]).cross(vertices[2] - vertices[0]);
                        Eigen::Vector3d n1 = (vertices[2] - vertices[1]).cross(vertices[3] - vertices[1]);
                        Eigen::Vector3d n2 = (vertices[3] - vertices[2]).cross(vertices[0] - vertices[2]);
                        Eigen::Vector3d n3 = (vertices[0] - vertices[3]).cross(vertices[1] - vertices[3]);
            
                        if (n0.dot(n1) > 0 && n0.dot(n2) > 0 && n0.dot(n3) > 0) {
                            tracer.surfaces.push_back(pair<Surface *, Fill>(
                                new Triangle(vertices[0], vertices[1], vertices[2]), currFill));
                            tracer.surfaces.push_back(pair<Surface *, Fill>(
                                new Triangle(vertices[0], vertices[2], vertices[3]), currFill));
                        }
                    }
                }
                
                break;
            }
            
            

            case 'l':{
                stringstream ss(line);
                Light l;
                string junk;
                ss >> junk >> l.p[0] >> l.p[1] >> l.p[2];
                if (!ss.eof()) {
                    ss>>l.c[0]>>l.c[1]>>l.c[2];
                    //coloredlights = true;
                }
                tracer.lights.push_back(l);
                break;
            }

            default:{
                break;
            }
        }
    }


    file.close();
    int inputSize = tracer.resolution[0] * tracer.resolution[1];
    tracer.im = new Eigen::Vector3d[inputSize];
    

    createImg(tracer);
    writeImg(tracer, outputFile);
}


Eigen::Vector3d shade(Tracer& tracer, const HitRecord &hr, int depth){
    Eigen::Vector3d surfaceColor = hr.f.rgb;

    double k_d = hr.f.kd;  
    double k_s = hr.f.ks;  
    double p = hr.f.e;     
    double k_t = hr.f.kt;  //unneeded for base
    double i_r = hr.f.ir;  //unneeded for base

    Eigen::Vector3d normal = hr.n; 
    normal.normalize();
    
    Eigen::Vector3d viewDir =  hr.v; 
    viewDir.normalize();
    
    Eigen::Vector3d ambientLight(0.0, 0.0, 0.0);

    Eigen::Vector3d shadedColor = ambientLight;

    

    for (const auto &light : tracer.lights) {
        Eigen::Vector3d lightDir = (light.p - hr.p).normalized();
        Eigen::Vector3d H = (viewDir + lightDir).normalized();

        Ray shadowRay(hr.p + normal * shadowBias, lightDir);
        HitRecord shadowHR;

        bool isShadowed = false;

        for (unsigned int i = 0; i < tracer.surfaces.size(); i++) {
            const auto &surfacePair = tracer.surfaces[i];
            const Surface *surface = surfacePair.first;
            
            if (surface->intersect(shadowRay, shadowBias, numeric_limits<double>::infinity(), shadowHR)) {
                isShadowed = true;
                break;
            }
        }


        if(!isShadowed){
            double diffuse = max(0.0, normal.dot(lightDir));

            double specular = pow(max(0.0, normal.dot(H)), p);

            double lightIntensity = 1.0 / sqrt(tracer.lights.size()); // maybe remove -1

            shadedColor[0] += (k_d * surfaceColor[0] * diffuse + k_s * specular) * lightIntensity;
            shadedColor[1] += (k_d * surfaceColor[1] * diffuse + k_s * specular) * lightIntensity;
            shadedColor[2] += (k_d * surfaceColor[2] * diffuse + k_s * specular) * lightIntensity;
            shadedColor[0] = std::min(1.0, std::max(0.0, shadedColor[0]));
            shadedColor[1] = std::min(1.0, std::max(0.0, shadedColor[1]));
            shadedColor[2] = std::min(1.0, std::max(0.0, shadedColor[2]));

        }
    }

        if (depth < 5 && k_s > 0.0) {
            Eigen::Vector3d reflectionDir =  2.0 * normal.dot(viewDir) * normal - viewDir;
            Ray reflectionRay(hr.p + normal * 0.001, reflectionDir); 
            HitRecord reflectionHR;

            Eigen::Vector3d reflectedColor = castRay(tracer, reflectionRay, 0.001, numeric_limits<double>::infinity(), depth + 1);

            shadedColor += k_s * reflectedColor;
        }



    return shadedColor;
}
  
Eigen::Vector3d castRay(Tracer& tracer, const Ray &r, double t0, double t1, int depth){
    HitRecord hr;
    Eigen::Vector3d color(tracer.bcolor);
    
    bool hit = false;
    for (unsigned int k=0; k<tracer.surfaces.size(); k++) {
      const pair<Surface *, Fill> &s  = tracer.surfaces[k];
      if (s.first->intersect(r, t0, t1, hr)) {
        t1 = hr.t;
        hr.f = s.second;
        hr.raydepth = r.depth;
        hr.v = r.e - hr.p;
        hr.v.normalize();
        hit = true;
      }
    }
  
    if (hit) {
        color = shade(tracer, hr, depth);

        // if (hr.f.ks > 0.0 && depth < 5) {
        //     Eigen::Vector3d normal = hr.n;
        //     Eigen::Vector3d viewDir = -hr.v;
        //     Eigen::Vector3d reflectionDir = viewDir - 2.0 * normal.dot(viewDir) * normal;
        //     Ray reflectionRay(hr.p + normal * 1e-6, reflectionDir); 
            
        //     Eigen::Vector3d reflectedColor = castRay(tracer, reflectionRay, 1e-6, numeric_limits<double>::infinity(), depth + 1);

        //     //cout << hr.f.ks << endl;
        //     color += hr.f.ks * reflectedColor; //hr.f.ks
        // }
    }

    return color;
  }

void createImg(Tracer& tracer){
    Eigen::Vector3d W = (tracer.from - tracer.at).normalized();
    Eigen::Vector3d U = tracer.up.cross(W).normalized();
    Eigen::Vector3d V = W.cross(U).normalized();

    double d = (tracer.from - tracer.at).norm();
    double h = tan((M_PI/180.0) * (tracer.angle/2.0)) * d;
    double increment = (2*h) / tracer.resolution[0];
    double l = -h + 0.5*increment;
    double t = h*(((double)tracer.resolution[1])/tracer.resolution[0]) - 0.5*increment;

    int samples = 1;

    Eigen::Vector3d *pixel = tracer.im;
    int counter = 0;

    for (unsigned int j=0; j<tracer.resolution[1]; j++) {
        for (unsigned int i=0; i<tracer.resolution[0]; i++, pixel++) {
          double x = l + i*increment;
          double y = t - j*increment;
          Eigen::Vector3d dir = x*U + y*V - d*W;
          (*pixel) = Eigen::Vector3d(0.0,0.0,0.0);
          for (int k=0; k<samples; k++) {
            for (int l=0; l<samples; l++) {
            Eigen::Vector3d origin = tracer.from;
            Eigen::Vector3d imagept = tracer.from + dir;
            Ray r(origin, imagept-origin);
            r.d.normalize();
            (*pixel) += castRay(tracer, r, tracer.hither, MAX, 0) / (samples*samples);
            //test
            counter++;
            }
          }
        }
        cout << "Percent Complete: " << double(counter / (tracer.resolution[1] * tracer.resolution[0])) << endl;
    }

}

void writeImg(Tracer& tracer, char* outputFile){
    std::ofstream outfile(outputFile, std::ios::binary); // Open in binary mode
    if (outfile.is_open()) {
        outfile << "P6\n" << tracer.resolution[0] << " " << tracer.resolution[1] << "\n255\n";
        for (int i = 0; i < (tracer.resolution[0] * tracer.resolution[1]); i++) {
            unsigned char r = static_cast<unsigned char>(min(tracer.im[i][0] * 255, 255.0));
            unsigned char g = static_cast<unsigned char>(min(tracer.im[i][1] * 255, 255.0));
            unsigned char b = static_cast<unsigned char>(min(tracer.im[i][2] * 255, 255.0));
            outfile.write(reinterpret_cast<char*>(&r), 1);
            outfile.write(reinterpret_cast<char*>(&g), 1);
            outfile.write(reinterpret_cast<char*>(&b), 1);
        }
    }
    outfile.close();
}


void baseProj(char* inputFile, char* outputFile){
    Tracer tracer;
    //tracer.ec = ec;
    parseFile(tracer, inputFile, outputFile);
}

int main(int argc, char* argv[]){
    

    if(argc == 4 && !strcmp(argv[1], "-p")){
        EC = true;
        
        baseProj(argv[2], argv[3]);
        
    }else{
        baseProj(argv[1], argv[2]);
    }
    
    cout << "Program finished." << endl;
    return 0;
}