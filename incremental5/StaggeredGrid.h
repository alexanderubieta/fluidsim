
#include <iostream>
#include <math.h>

#include <vector>

#include <Eigen/Dense>

#include "Particle.h"

#include "PressureSolver.h"

class StaggeredGrid {
  public:

    // dimensions
    Eigen::Vector3d dims;

    // pressure
    Array3D<double> sp;

    // vel
    Array3D<double> su;
    Array3D<double> sv;
    Array3D<double> sw;

    // accum
    Array3D<double> fu;
    Array3D<double> fv;
    Array3D<double> fw;

    const Eigen::Vector3d lc;

    const double dx;

    // Half-grid-cell-width shifts for splatting particle values onto the grid
    const Eigen::Vector3d shift_yz_;  // (0, dx_/2, dx_/2)
    const Eigen::Vector3d shift_xz_;  // (dx_/2, 0, dx_/2)
    const Eigen::Vector3d shift_xy_;  // (dx_/2, dx_/2, 0)

    Array3D<MaterialType> labels;

    Array3D<unsigned short> neighbors;

    StaggeredGrid(int gx, int gy, int gz, Eigen::Vector3d glc, double gdx) : 
    dims(gx,gy,gz), sp(gx,gy,gz), 
    su(gx+1,gy,gz), sv(gx,gy+1,gz), 
    sw(gx,gy,gz+1), fu(gx+1,gy,gz), 
    fv(gx,gy+1,gz), fw(gx,gy,gz+1),
    lc(glc), dx(gdx), labels(gx, gy, gz),
    shift_yz_(0.0,dx/2.0,dx/2.0), shift_xz_(dx/2.0,0.0,dx/2.0),
    shift_xy_(dx/2.0, dx/2.0, 0.0), neighbors(gx,gy,gz) {
        // reset();
    }

    const Array3D<double>& p() const{
        return sp;
    }

    const Array3D<double>& u() const{
        return su;
    }

    const Array3D<double>& v() const{
        return sv;
    }

    const Array3D<double>& w() const{
        return sw;
    }

    const Array3D<MaterialType>& cell_labels() const{
        return labels;
    }

    const int nx(){
        return dims[0];
    }

    const int ny(){
        return dims[1];
    }

    const int nz(){
        return dims[2];
    }

    void ParticlesToGrid(const std::vector<Particle>& particles){
        reset();
        for (std::vector<Particle>::const_iterator p = particles.begin();
                p != particles.end(); p++) {
            Eigen::Vector3d pos = p->pos - lc;
            Eigen::Vector3d vel = p->vel;
            Eigen::Vector3d ijk = floor(pos,dx);
            labels(ijk[0],ijk[1],ijk[2]) = MaterialType::FLUID;

            Eigen::Vector3d pshift = pos - shift_yz_;
            Eigen::Vector3d pshift1 = pos - shift_xz_;
            Eigen::Vector3d pshift2 = pos - shift_xy_;
            
            Eigen::Vector3d cell = floor(pshift,dx);
            Eigen::Vector3d bary = weights(pshift/dx,cell);
            splat(vel[0],bary,cell,&su,&fu);
            cell = floor(pshift1,dx);
            bary = weights(pshift1/dx,cell);
            splat(vel[1],bary,cell,&sv,&fv);
            cell = floor(pshift2,dx);
            bary = weights(pshift2/dx,cell);
            splat(vel[2],bary,cell,&sw,&fw);
        }
        normalize();
        boundarySplat();
    }

    const Eigen::Vector3d Advect(Eigen::Vector3d pos, double t){
        Eigen::Vector3d pshift = pos - shift_yz_;
        Eigen::Vector3d pshift1 = pos - shift_xz_;
        Eigen::Vector3d pshift2 = pos - shift_xy_;

        Eigen::Vector3d cell = floor(pshift - lc,dx);
        Eigen::Vector3d bary = weights(pshift/dx,cell);

        double ru = (1-bary[0])*(1-bary[1])*(1-bary[2])*u()(cell[0],cell[1],cell[2]) + 
        (1-bary[0])*(1-bary[1])*(bary[2])*u()(cell[0],cell[1],cell[2]+1) + 
        (1-bary[0])*(bary[1])*(1-bary[2])*u()(cell[0],cell[1]+1,cell[2]) + 
        (1-bary[0])*(bary[1])*(bary[2])*u()(cell[0],cell[1]+1,cell[2]+1) + 
        (bary[0])*(1-bary[1])*(1-bary[2])*u()(cell[0]+1,cell[1],cell[2]) +
        (bary[0])*(1-bary[1])*(bary[2])*u()(cell[0]+1,cell[1],cell[2]+1) + //-
        (bary[0])*(bary[1])*(1-bary[2])*u()(cell[0]+1,cell[1]+1,cell[2]) + 
        (bary[0])*(bary[1])*(bary[2])*u()(cell[0]+1,cell[1]+1,cell[2]+1);

        cell = floor(pshift1 - lc,dx);
        bary = weights(pshift1/dx,cell);
        double rv = (1-bary[0])*(1-bary[1])*(1-bary[2])*v()(cell[0],cell[1],cell[2]) + 
        (1-bary[0])*(1-bary[1])*(bary[2])*v()(cell[0],cell[1],cell[2]+1) + 
        (1-bary[0])*(bary[1])*(1-bary[2])*v()(cell[0],cell[1]+1,cell[2]) + 
        (1-bary[0])*(bary[1])*(bary[2])*v()(cell[0],cell[1]+1,cell[2]+1) + 
        (bary[0])*(1-bary[1])*(1-bary[2])*v()(cell[0]+1,cell[1],cell[2]) +
        (bary[0])*(1-bary[1])*(bary[2])*v()(cell[0]+1,cell[1],cell[2]+1) + //-
        (bary[0])*(bary[1])*(1-bary[2])*v()(cell[0]+1,cell[1]+1,cell[2]) + 
        (bary[0])*(bary[1])*(bary[2])*v()(cell[0]+1,cell[1]+1,cell[2]+1);

        cell = floor(pshift2 - lc,dx);
        bary = weights(pshift2/dx,cell);
        double rw = (1-bary[0])*(1-bary[1])*(1-bary[2])*w()(cell[0],cell[1],cell[2]) + 
        (1-bary[0])*(1-bary[1])*(bary[2])*w()(cell[0],cell[1],cell[2]+1) + 
        (1-bary[0])*(bary[1])*(1-bary[2])*w()(cell[0],cell[1]+1,cell[2]) + 
        (1-bary[0])*(bary[1])*(bary[2])*w()(cell[0],cell[1]+1,cell[2]+1) + 
        (bary[0])*(1-bary[1])*(1-bary[2])*w()(cell[0]+1,cell[1],cell[2]) +
        (bary[0])*(1-bary[1])*(bary[2])*w()(cell[0]+1,cell[1],cell[2]+1) + //-
        (bary[0])*(bary[1])*(1-bary[2])*w()(cell[0]+1,cell[1]+1,cell[2]) + 
        (bary[0])*(bary[1])*(bary[2])*w()(cell[0]+1,cell[1]+1,cell[2]+1);

        Eigen::Vector3d res(ru,rv,rw);

        Eigen::Vector3d n = pos + t*res;
        if(n[0] == 0){
            n[0] += 1;
        }
        if(n[0] == nx() - 1){
            n[0] -= 1;
        }
        if(n[1] == 0){
            n[1] += 1;
        }
        if(n[1] == ny() - 1){
            n[1] -= 1;
        }
        if(n[2] == 0){
            n[2] += 1;
        }
        if(n[2] == nz() - 1){
            n[2] -= 1;
        }
        return n;
    }

    void ApplyGravity(double t){
        for(int i = 0; i < nx(); i++){
            for(int j = 0; j <= ny(); j++){
                for(int k = 0; k < nz(); k++){
                    sv(i,j,k) -= t*9.80665;
                }
            }
        }
        boundarySplat();
    }

    void ProjectPressure(){
        // use pressure sim
        sp = 0;
        updateNeighbors();
        PressureSolver pSolve;
        pSolve.ProjectPressure(labels, neighbors, su, sv, sw, &sp);
        ApplyPressure();
    }
  private:
    // set all to zero and set material type to default values
    void reset(){
        su = 0.0;
        sv = 0.0;
        sw = 0.0;
        fu = 0.0;
        fv = 0.0;
        fw = 0.0;
        for(int i = 0; i<nx(); i++){
            for(int j = 0; j<ny(); j++){
                for(int k = 0; k<nz(); k++){
                    if(i == 0 || i == nx() - 1 || j == 0 || j == ny() - 1 || k == 0 || k == nz() - 1){
                        labels(i,j,k) = MaterialType::SOLID;
                    }
                    else{
                        labels(i,j,k) = MaterialType::EMPTY;
                    }
                }
            }
        }
    }

    inline Eigen::Vector3d floor(Eigen::Vector3d p_lc, double dx){
        assert(dx>0.0);
        Eigen::Vector3d vec = p_lc / dx;
        assert(vec[0]>=0.0);
        assert(vec[1]>=0.0);
        assert(vec[2]>=0.0);
        // std::cout <<  "pos " << vec << std::endl;  
        return vec;
    }

    inline Eigen::Vector3d weights(Eigen::Vector3d p_lc_dx, Eigen::Vector3d cell){
       
        // std::cout <<  "weig p-lc" << (p-lc) << std::endl;  
        // std::cout <<  "weig dx" << dx << std::endl;  
        // std::cout <<  "weigh" << p_lc_dx - cell << std::endl;  
        int i = cell[0];
        int j = cell[1];
        int k = cell[2];
        Eigen::Vector3d rhs(i,j,k);
        return p_lc_dx - rhs;
    }

    inline void splat(double v, Eigen::Vector3d& bary,Eigen::Vector3d& cell, Array3D<double>* grid, Array3D<double>* f){
        // std::cout <<  "ori " << (1-bary[0])*(1-bary[1])*(1-bary[2])*v << "i"<<cell[0] << "j"<<cell[1] << "k"<<cell[2] << std::endl;  
        (*grid)(cell[0],cell[1],cell[2]) += (1-bary[0])*(1-bary[1])*(1-bary[2])*v;
        (*f)(cell[0],cell[1],cell[2]) += (1-bary[0])*(1-bary[1])*(1-bary[2]);
        (*grid)(cell[0] + 1,cell[1],cell[2]) += (bary[0])*(1-bary[1])*(1-bary[2])*v;
        (*f)(cell[0] + 1,cell[1],cell[2]) += (bary[0])*(1-bary[1])*(1-bary[2]);
        (*grid)(cell[0],cell[1] + 1,cell[2]) += (1-bary[0])*(bary[1])*(1-bary[2])*v;
        (*f)(cell[0],cell[1] + 1,cell[2]) += (1-bary[0])*(bary[1])*(1-bary[2]);
        (*grid)(cell[0] + 1,cell[1] + 1,cell[2]) += (bary[0])*(bary[1])*(1-bary[2])*v;
        (*f)(cell[0] + 1,cell[1] + 1,cell[2]) += (bary[0])*(bary[1])*(1-bary[2]);
        (*grid)(cell[0],cell[1],cell[2] + 1) += (1-bary[0])*(1-bary[1])*(bary[2])*v;
        (*f)(cell[0],cell[1],cell[2] + 1) += (1-bary[0])*(1-bary[1])*(bary[2]);
        (*grid)(cell[0] + 1,cell[1],cell[2] + 1) += (bary[0])*(1-bary[1])*(bary[2])*v;
        (*f)(cell[0] + 1,cell[1],cell[2] + 1) += (bary[0])*(1-bary[1])*(bary[2]);
        (*grid)(cell[0],cell[1] + 1,cell[2] + 1) += (1-bary[0])*(bary[1])*(bary[2])*v;
        (*f)(cell[0],cell[1] + 1,cell[2] + 1) += (1-bary[0])*(bary[1])*(bary[2]);
        (*grid)(cell[0] + 1,cell[1] + 1,cell[2] + 1) += (bary[0])*(bary[1])*(bary[2])*v;
        (*f)(cell[0] + 1,cell[1] + 1,cell[2] + 1) += (bary[0])*(bary[1])*(bary[2]);
    }

    inline void normalize(){
        for(int i = 0; i<=nx(); i++){
            for(int j = 0; j<=ny(); j++){
                for(int k = 0; k<=nz(); k++){
                    if(j != ny() && k != nz()){
                        if(i != 0 && i != 1 && i != nx() - 1 && i != nx()){
                            if(fu(i,j,k) > 0.00001){
                                su(i,j,k) /= fu(i,j,k);
                            }
                        }
                    }
                    if(i != nx() && k != nz()){
                        if(j != 0 && j != 1 && j != ny() - 1 && j != ny()){
                            if(fv(i,j,k) > 0.00001){
                                sv(i,j,k) /= fv(i,j,k);
                            }
                        }
                    }
                    if(i != nx() && j != ny()){
                        if(k != 0 && k != 1 && k != nz() - 1 && k != nz()){
                            if(fw(i,j,k) > 0.00001){
                                sw(i,j,k) /= fw(i,j,k);
                            }
                        }
                    }

                    // if these ifs fail turn to 0
                }
            }
        }
    }

    inline void boundarySplat(){
        for(int j = 0; j <= ny(); j++){
            for(int k = 0; k <= nz(); k++){
                if(k != nz()){
                    // v
                    sv(0,j,k) = sv(1,j,k);
                    sv(nx()-1,j,k) = sv(nx()-2,j,k);
                }
                if(j != ny()){
                    // w
                    sw(0,j,k) = sw(1,j,k);
                    sw(nx()-1,j,k) = sw(nx()-2,j,k);
                }
                if(j != ny() && k != nz()){
                    su(0,j,k) = 0.0;
                    su(1,j,k) = 0.0;
                    su(nx(),j,k) = 0.0;
                    su(nx()-1,j,k) = 0.0;
                }
            }
        }
        for(int i = 0; i <= nx(); i++){
            for(int k = 0; k <= nz(); k++){
                if(k != nz()){
                    // u
                    su(i,0,k) = su(i,1,k);
                    su(i,ny()-1,k) = su(i,ny()-2,k);
                }
                if(i != nx()){
                    // w
                    sw(i,0,k) = sw(i,1,k);
                    sw(i,ny()-1,k) = sw(i,ny()-2,k);
                }
                if(i != nx() && k != nz()){
                    sv(i,0,k) = 0.0;
                    sv(i,1,k) = 0.0;
                    sv(i,ny(),k) = 0.0;
                    sv(i,ny()-1,k) = 0.0;
                }
            }
        }
        for(int i = 0; i <= nx(); i++){
            for(int j = 0; j <= ny(); j++){
                if(j != ny()){
                    // u
                    su(i,j,0) = su(i,j,1);
                    su(i,j,nz()-1) = su(i,j,nz()-2);
                }
                if(i != nx()){
                    // v
                    sv(i,j,0) = sv(i,j,1);
                    sv(i,j,nz()-1) = sv(i,j,nz()-2);
                }
                if(i != nx() && j != ny()){
                    sw(i,j,0) = 0.0;
                    sw(i,j,1) = 0.0;
                    sw(i,j,nz()) = 0.0;
                    sw(i,j,nz()-1) = 0.0;
                }
            }
        }
    }

    inline void updateNeighbors(){
        neighbors = 0u;
        for(int i = 1; i<nx()-1; i++){
            for(int j = 1; j<ny()-1; j++){
                for(int k = 1; k<nz()-1; k++){
                    unsigned short store = 0u;
                    int count = 0;
                    if(labels(i,j,k+1) != SOLID) count++;
                    if(labels(i,j,k+1) == FLUID){
                        store += 1<<8;
                    }
                    if(labels(i,j+1,k) != SOLID) count++;
                    if(labels(i,j+1,k) == FLUID){
                        store += 1<<7;
                    }
                    if(labels(i+1,j,k) != SOLID) count++;
                    if(labels(i+1,j,k) == FLUID){
                        store += 1<<6;
                    }
                    if(labels(i,j,k-1) != SOLID) count++;
                    if(labels(i,j,k-1) == FLUID){
                        store += 1<<5;
                    }
                    if(labels(i,j-1,k) != SOLID) count++;
                    if(labels(i,j-1,k) == FLUID){
                        store += 1<<4;
                    }
                    if(labels(i-1,j,k) != SOLID) count++;
                    if(labels(i-1,j,k) == FLUID){
                        store += 1<<3;
                    }
                    neighbors(i,j,k) = store | count;
                }
            }
        }
    }

    inline void ApplyPressure(){
        for(int i = 1; i<nx()-1; i++){
            for(int j = 1; j<ny()-1; j++){
                for(int k = 1; k<nz()-1; k++){
                    unsigned short store = 0;
                    if(labels(i,j,k) != SOLID){
                        if(labels(i-1,j,k) != SOLID){
                            su(i,j,k) -= sp(i,j,k) - sp(i-1,j,k);
                        }
                        if(labels(i,j-1,k) != SOLID){
                            sv(i,j,k) -= sp(i,j,k) - sp(i,j-1,k);
                        }
                        if(labels(i,j,k-1) != SOLID){
                            sw(i,j,k) -= sp(i,j,k) - sp(i,j,k-1);
                        }
                    }
                }
            }
        }
    }
};