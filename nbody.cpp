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

// Calculate the acceleration
void calc_acceleration(std::vector<Particle>& data, const double eps) {
    // Clear the acceleration with zero
    const int n = data.size();
    for (int i = 0; i < n; ++i) {
            data[i].ax = 0.0;
            data[i].ay = 0.0;
            data[i].az = 0.0;
    }
    double m  = data[0].m; // assume equal mass
    for (int i = 0; i < n; ++i) {  // calculate acceleration for all i-particles
        for (int j = 0; j < n; ++j) { // accumulate over all j-particles
        // *** implement by yourself ***
            if(i==j)continue;
            double dx = data[j].x - data[i].x;
            double dy = data[j].y - data[i].y;
            double dz = data[j].z - data[i].z;
            double r = sqrt(dx*dx + dy*dy + dz*dz + eps*eps);
            data[i].ax += dx*data[j].m / (r*r*r);
            data[i].ay += dy*data[j].m / (r*r*r);
            data[i].az += dz*data[j].m / (r*r*r);
        }
    }
    return;
}


// The leapfrog integrator
// Use the Kick-Drift-Kick formula
void leapfrog(std::vector<Particle>& data, const double dt, const double eps) {
    // Assume that ax, ay, az contain the accleration at t=tn calculated in the last step
    // Note: we need to set them by calc_acceleration for the very first step (see main())
    const int n = data.size();

    // push the particle velocities
    for (int i = 0; i < n; ++i) {
        // *** implement by yourself ***
        data[i].vx += data[i].ax * dt/2.;
        data[i].vy += data[i].ay * dt/2.;
        data[i].vz += data[i].az * dt/2.;
    }

    // push the particle positions
    for (int i = 0; i < n; ++i) {
        // *** implement by yourself ***
        data[i].x += data[i].vx * dt;
        data[i].y += data[i].vy * dt;
        data[i].z += data[i].vz * dt;
    }

    // calculate the acceleration at t=tn+1
    calc_acceleration(data, eps);

    // push the particle velocities again
    for (int i = 0; i < n; ++i) {
        // *** implement by yourself ***
        data[i].vx += data[i].ax * dt/2.;
        data[i].vy += data[i].ay * dt/2.;
        data[i].vz += data[i].az * dt/2.;
    }
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

    // calculate the acceleration at t=0 before starting the main loop
    calc_acceleration(data, eps);

    for (int i = 0; i < NSTEP; ++i) { // main loop
        if ((i % OSTEP) == 0) { // output the data every OSTEP step
            output(data, i/OSTEP, folder_name);
            radial_profile(data, i/OSTEP, folder_name);
        }

        leapfrog(data, dt, eps); // integrate dt using leapfrog
        t += dt;

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
        run(N, NSTEP*i, OSTEP*i, eps, rad, rv, tend, std::to_string(i));
    }
    Random::engine.seed(1);

    return 0;
}

