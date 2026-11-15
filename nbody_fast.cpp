#include <cassert>
#include <cstdlib>
#include <ios>
#include <iostream>
#include <cmath>
#include <random>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <filesystem>

namespace Random{
    const int rand_seed = 1;
    std::default_random_engine engine(rand_seed); 
    std::uniform_real_distribution<double> uni_0_1(0., 1.);
}

class Particle{
public:
    double m, x, y, z, vx, vy, vz, ax, ay, az, r;

    Particle():
        m{}, x{}, y{}, z{}, vx{}, vy{}, vz{}, ax{}, ay{}, az{}, r{}
    {}
};

class ParticleTree{
private:
    const int N;
    const int num_particle;
    std::vector<Particle> data;

    int l_child_id(int parent_id)const{
        //return (parent+1)*2-1;
        return 2*parent_id+1;
    }
    int r_child_id(int parent_id)const{
        //return (parent+1)*2+1-1;
        return 2*parent_id+2;
    }
    int parent_id(int child_id)const{
        return (child_id-1)/2;
    }
public:
    ParticleTree(const std::vector<Particle>& particles):
        N(
            //粒子数以上の最小の2のべき乗
            [&](){
                int size=1;
                while(size<particles.size())size*=2;
                return size;
            }()
        ),
        data(2*N),
        num_particle(particles.size())
    {
        //[0]集団の質量和と重心、平均速度を計算
        for(int i=0;i<num_particle;i++){
            data[N+i] = particles[i];
        }
        for(int i=N-1;i>=0;i--){
            Particle& parent = data[i];
            const Particle& left = data[l_child_id(i)];
            const Particle& right = data[r_child_id(i)];
            const double ml=left.m;
            const double mr=right.m;
            const double total = ml+mr;
            if(total<1e-40)continue;//子供が両方粒子を持っていない場合、こうしないとバグる
            const double l=ml/total;
            const double r=mr/total;
            parent.m = total;
            parent.x = left.x * l + right.x * r;
            parent.y = left.y * l + right.y * r;
            parent.z = left.z * l + right.z * r;
            parent.vx = left.vx * l + right.vx * r;
            parent.vy = left.vy * l + right.vy * r;
            parent.vz = left.vz * l + right.vz * r;
        }
        //[1]親と比較した相対速度、相対位置を計算
        for(int i=2*N-1;i>0;i--){
            Particle& child = data[i];
            const Particle& parent = data[parent_id(i)];

            child.x -= parent.x;
            child.y -= parent.y;
            child.z -= parent.z;
            child.vx -= parent.vx;
            child.vy -= parent.vy;
            child.vz -= parent.vz;
        }
        for(int i=0;i<N;i++){
            if(i%20==0){
                std::cout<<"\n";
            }
            std::cout<<data[i].vx<<" ";
            
        }
    }

    //絶対位置・絶対速度を取り出す
    void extract_status(std::vector<Particle>& result)const{
        assert(num_particle==result.size());
        std::vector<Particle> copy_data(data);
        for(int i=1;i<2*N;i++){
            Particle& child = copy_data[i];
            const Particle& parent = copy_data[parent_id(i)];
            child.x += parent.x;
            child.y += parent.y;
            child.z += parent.z;
            child.vx += parent.vx;
            child.vy += parent.vy;
            child.vz += parent.vz;
        }
        for(int i=0;i<num_particle;i++){
            result[i]=copy_data[N+i];
        }
    }

    void calc_acceleration(const double eps){
        for(int i=0;i<N;i++){
            Particle& left = data[l_child_id(i)];
            Particle& right = data[r_child_id(i)];
            const double ml = left.m;
            const double mr = right.m;
            if(ml+mr < 1e-40)continue;

            const double dx = right.x - left.x;
            const double dy = right.y - left.y;
            const double dz = right.z - left.z;
            const double r = sqrt(dx*dx + dy*dy + dz*dz + eps*eps);
            
            const double fx = dx/(r*r*r);
            const double fy = dy/(r*r*r);
            const double fz = dz/(r*r*r);

            left.ax = mr*fx;
            left.ay = mr*fy;
            left.az = mr*fz;
            right.ax = -ml*fx;
            right.ay = -ml*fy;
            right.az = -ml*fz;
        }
    }

    void push_position(const double dt){
        for(int i=0;i<2*N;i++){
            Particle& target = data[i];
            target.x += target.vx * dt;
            target.y += target.vy * dt;
            target.z += target.vz * dt;
        }
    }

    void push_velocity(const double dt){
        for(int i=0;i<2*N;i++){
            Particle& target = data[i];
            target.vx += target.ax * dt;
            target.vy += target.ay * dt;
            target.vz += target.az * dt;
        }
    }

};

// Generate a Gaussian random number using the Box Muller method
double bmrand() {
    const double x = Random::uni_0_1(Random::engine);
    const double y = Random::uni_0_1(Random::engine);
    // two numbers are generated, but use only one of them
    return std::sqrt(-2.0*std::log(x))*std::cos(2.0*M_PI*y);
}

// Calculate the total energy
double total_ekinetic(const std::vector<Particle>& data) {
    double ekin = 0.0;
    const int n = data.size();
    for (int i = 0; i < n; ++i) {
        ekin += 0.5*data[i].m*(data[i].vx*data[i].vx
            +data[i].vy*data[i].vy+data[i].vz*data[i].vz);
    }
    return ekin;
}

// Calculate the total potential energy
double total_potential(const std::vector<Particle>& data, const double eps) {
    double pot = 0.0;
    const int n = data.size();
    // Accumulate the potential for all the pairs of the particles
    // Be careful of the loop limits not to double count
    for (int i = 0; i < n-1; ++i) {
        for (int j = i+1; j < n; ++j) {
            double dx = data[j].x - data[i].x;
            double dy = data[j].y - data[i].y;
            double dz = data[j].z - data[i].z;
            double r = sqrt(dx*dx + dy*dy + dz*dz + eps*eps);
            pot -= data[i].m*data[j].m/r;
        }
    }
    return pot;
}

// Calculate the initial particle distribution
// Particles are randomly distributed within a uniform sphere
void init_position(std::vector<Particle>& data, const double rad) {
    assert(rad>=0.);
    const int n = data.size();
    std::uniform_real_distribution<double> dist(-rad, rad);
    for (int i = 0; i < n; ++i) {
        double x, y, z;
        do { // the rejection method
            x = dist(Random::engine);
            y = dist(Random::engine);
            z = dist(Random::engine);
        } while(x*x+y*y+z*z > rad);
        data[i].x = x;
        data[i].y = y;
        data[i].z = z;

        /*
        data[i].x = 0.;
        data[i].y = 0.;
        data[i].z = 0.;*/
    }
    return;
}
// Calculate the initial velocity distribution
void init_velocity(std::vector<Particle>& data, const double rv, const double pot) {
    const double sigma = sqrt(2.0*rv*fabs(pot)/3.0);
    // vx, vy, vz follows the Maxwellian distribution
    // with the average = 0 and dispersion = sigma
    const int n = data.size();
    for (int i = 0; i < n; ++i) {
        data[i].vx = sigma*bmrand();
        data[i].vy = sigma*bmrand();
        data[i].vz = sigma*bmrand();

        /*data[i].vx = 0.;
        data[i].vy = 0.;
        data[i].vz = 0.;*/
    }
    return;
}

// Generate the initial particle data
void init_particles(std::vector<Particle>& data, const double eps,
                    double rad, double rv) {
    const int n = data.size();
    for (int i = 0; i < n; ++i)
        data[i].m = 1.0/(double)n; // set the mass - equal mass
    init_position(data, rad);
    double pot = total_potential(data, eps);
    init_velocity(data, rv, pot);
    return;
}

// The leapfrog integrator
// Use the Kick-Drift-Kick formula
void leapfrog(ParticleTree& tree, const double dt, const double eps) {
    // Assume that ax, ay, az contain the accleration at t=tn calculated in the last step
    // Note: we need to set them by calc_acceleration for the very first step (see main())

    // push the particle velocities
    tree.push_velocity(dt/2.);

    // push the particle positions
    tree.push_position(dt);

    // calculate the acceleration at t=tn+1
    tree.calc_acceleration(eps);

    // push the particle velocities again
    tree.push_velocity(dt/2.);
    return;
}

// data output 
void output(const std::vector<Particle>& data, int step, std::string& folder_name) {
    std::string filename = [&](){
        char fn[32];
        sprintf(fn, "output%04d.dat", step);
        return folder_name + "/" + std::string(fn);
    }();
    std::ofstream ofs(filename,std::ios::binary);

    const int n = data.size();
    for (int i = 0; i < n; ++i) {
        const Particle& target = data[i];
        ofs
        <<target.m
        <<target.x
        <<target.y
        <<target.z
        <<target.vx
        <<target.vy
        <<target.vz;
    }

    ofs.close();
    return;
}

// calculate the particle radii from the center of the domain for analysis
void calc_radii(std::vector<Particle>& data) {
    const int n = data.size();
    for (int i = 0; i < n; ++i)
        data[i].r = sqrt(data[i].x*data[i].x
                + data[i].y*data[i].y + data[i].z*data[i].z);
    return;
}

void radial_profile(std::vector<Particle>& data, const int step, std::string& folder_name){
    
    std::string filename = [&](){
        char fn[32];
        sprintf(fn, "radprof%04d.dat", step);  
        return folder_name + "/" + std::string(fn);
    }();

    std::ofstream ofs(filename,std::ios::binary);

    std::sort(data.begin(), data.end(), 
    [](const Particle& l, const Particle& r){
            return l.r < r.r;
        }
    );
    int d = 64;
    const int n = data.size();
    for (int i = 0; i < n-d; ++i) {
        // Pick a shell containing i-th to (i+d)th particles and
        // calculate the radius and density of this shell 
        double radius = 0.5 * (data[i].r+data[i+d].r);  // the radius of the shell
        double volume = 4.0*M_PI/3.0*(data[i+d].r*data[i+d].r*data[i+d].r
                    - data[i].r*data[i].r*data[i].r); // the volume of the shell
        double mass = d * data[i].m; // the mass of the shell
        double rho = mass / volume;  // the density of the shell
        ofs<<radius<<rho;
    }
}


// Calculate the total energy and virial ratio for analysis
void analyze(const std::vector<Particle>& data, const double eps, const double t, std::ofstream& ofs) {
    double ekin = total_ekinetic(data);
    double pot  = total_potential(data, eps);
    double etot = ekin + pot;
    double rv = ekin / fabs(pot);
    ofs<<
        [&](){
            char fn[128];
            sprintf(fn, "time = %f, Etot = %.16f, Rvir = %g\n", t, etot, rv);
            return std::string(fn);
        }();
    return;
}

void run(
    const int N, // the number of the particles
    const int NSTEP, // the number of the steps
    const int OSTEP,   // the output inteval
    const double eps,  // the softerning parameter
    const double rad,  // the radius of the initial sphere
    const double rv,  // the initial virial ratio
    const double tend,
    std::string folder_name
){
    Random::engine.seed(Random::rand_seed);
    std::filesystem::create_directory(folder_name);
    // Declare and allocate the data structure array
    std::vector<Particle> data(N);

    init_particles(data, eps, rad, rv);
    double t = 0.0, dt = tend / NSTEP;
    std::ofstream ofs(folder_name + "/analysis.txt");
    analyze(data, eps, t, ofs);

    double sum = 0.;
    for(int i=0;i<data.size();i++){
        sum+=data[i].vx;
    }
    std::cout<<"---"<<sum/(double)data.size()<<"--------\n\n\n";

    ParticleTree tree(data);

    // calculate the acceleration at t=0 before starting the main loop
    tree.calc_acceleration(eps);

    for (int i = 0; i < NSTEP; ++i) { // main loop
        if ((i % OSTEP) == 0) { // output the data every OSTEP step
            tree.extract_status(data);
            output(data, i/OSTEP, folder_name);
            radial_profile(data, i/OSTEP, folder_name);
        }

        leapfrog(tree, dt, eps); // integrate dt using leapfrog
        t += dt;

        tree.extract_status(data);
        analyze(data, eps, t, ofs); // do analysis
    }

    // output the final data
    output(data, NSTEP/OSTEP, folder_name); 
    radial_profile(data, NSTEP/OSTEP, folder_name);
}


int main(){
    const int N       = 4096; // the number of the particles
    const int NSTEP   = 1024; // the number of the steps
    const int OSTEP   = 64;   // the output inteval
    const double eps  = 0.1;  // the softerning parameter
    const double rad  = 1.0;  // the radius of the initial sphere
    const double rv   = 0.1;  // the initial virial ratio
    const double tend = 5.0;  // the time of the end of the simulation

    for(int i=1;i<17;i*=2){
        run(N, NSTEP*i, OSTEP*i, eps, rad, rv, tend, "fast"+std::to_string(i));
    }
    Random::engine.seed(1);

    return 0;
}

