
#include <iostream>
#include <math.h>

#include <vector>

#include "Array3D.h"

#include <Eigen/Dense>

#include "Particle.h"

class StaggeredGrid {
  public:
    enum MaterialType{
        SOLID, FLUID, EMPTY
    };

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

    Eigen::Vector3d lc;

    double dx;

    // Half-grid-cell-width shifts for splatting particle values onto the grid
    const Eigen::Vector3d shift_yz_;  // (0, dx_/2, dx_/2)
    const Eigen::Vector3d shift_xz_;  // (dx_/2, 0, dx_/2)
    const Eigen::Vector3d shift_xy_;  // (dx_/2, dx_/2, 0)

    Array3D<MaterialType> labels;

    StaggeredGrid(int gx, int gy, int gz, Eigen::Vector3d glc, double gdx) : 
    dims(gx,gy,gz), sp(gx,gy,gz), 
    su(gx+1,gy,gz), sv(gx,gy+1,gz), 
    sw(gx,gy,gz+1), fu(gx+1,gy,gz), 
    fv(gx,gy+1,gz), fw(gx,gy,gz+1),
    lc(glc), dx(gdx), labels(gx, gy, gz),
    shift_yz_(0.0,dx/2.0,dx/2.0), shift_xz_(dx/2.0,0.0,dx/2.0),
    shift_xy_(dx/2.0, dx/2.0, 0.0) {
        reset();
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

    // set all to zero and set material type to default values
    void reset(){
        // for(int i = 0; i<=nx(); i++){
        //     for(int j = 0; j<=ny(); j++){
        //         for(int k = 0; k<=nz(); k++){
        //             if(i == nx()){
        //                 if(j != ny() && k != nz()){
        //                     u(i,j,k) = 0.0;
        //                     (fu)*(i,j,k) = 0.0;
        //                 }
        //             }
        //             else if(j == ny()){
        //                 if(k != nz()){
        //                     v(i,j,k) = 0.0;
        //                     (fv)*(i,j,k) = 0.0;
        //                 }
        //             }
        //             else if(k == nz()){
        //                 w(i,j,k) = 0.0;
        //                 (fw)*(i,j,k) = 0.0;
        //             }
        //             else{
        //                 u(i,j,k) = 0.0;
        //                 (fu)*(i,j,k) = 0.0;
        //                 v(i,j,k) = 0.0;
        //                 (fv)*(i,j,k) = 0.0;
        //                 w(i,j,k) = 0.0;
        //                 (fw)*(i,j,k) = 0.0;
        //                 if(i == 0 || i == nx() - 1 || j == 0 || j == ny() - 1 || k == 0 || k == nz() - 1){
        //                     (labels)*(i,j,k) = SOLID;
        //                 }
        //                 else{
        //                     (labels)*(i,j,k) = EMPTY;
        //                 }
        //             }
        //         }
        //     }
        // }
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
            Eigen::Vector3d pos = p->pos;
            Eigen::Vector3d vel = p->vel;
            Eigen::Vector3d cell = floor(pos,lc,dx);
            labels(cell[0],cell[1],cell[2]) = MaterialType::FLUID;

            Eigen::Vector3d pshift(pos[0], pos[1] - dx / 2, pos[2] - dx / 2);
            Eigen::Vector3d pshift1(pos[0]- dx / 2, pos[1], pos[2] - dx / 2);
            Eigen::Vector3d pshift2(pos[0]- dx / 2, pos[1] - dx / 2, pos[2]);
            
            cell = floor(pshift,lc,dx);
            Eigen::Vector3d bary = weights(pshift,lc,cell,dx);
            splat(vel[0],bary,cell,&su,&fu);
            cell = floor(pshift1,lc,dx);
            bary = weights(pshift1,lc,cell,dx);
            splat(vel[1],bary,cell,&sv,&fv);
            cell = floor(pshift2,lc,dx);
            bary = weights(pshift2,lc,cell,dx);
            splat(vel[2],bary,cell,&sw,&fw);
        }
        normalize();
        boundarySplat();
    }
  private:
    inline Eigen::Vector3d floor(Eigen::Vector3d p, Eigen::Vector3d lc, double dx){
        assert(dx>0.0);
        Eigen::Vector3d vec = (p-lc) / dx;
        assert(vec[0]>=0.0);
        assert(vec[1]>=0.0);
        assert(vec[2]>=0.0);
        return vec;
    }

    inline Eigen::Vector3d weights(Eigen::Vector3d p, Eigen::Vector3d lc,Eigen::Vector3d cell, double dx){
        return (p-lc) / dx - cell;
    }

    inline void splat(double v, Eigen::Vector3d& bary,Eigen::Vector3d& cell, Array3D<double>* grid, Array3D<double>* f){
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
                    su(0,j,k) = 0;
                    su(1,j,k) = 0;
                    su(nx()-1,j,k) = 0;
                    su(nx()-2,j,k) = 0;
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
                    sv(i,0,k) = 0;
                    sv(i,1,k) = 0;
                    sv(i,ny()-1,k) = 0;
                    sv(i,ny()-2,k) = 0;
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
                    sw(i,j,0) = 0;
                    sw(i,j,1) = 0;
                    sw(i,j,nz()-1) = 0;
                    sw(i,j,nz()-2) = 0;
                }
            }
        }
    }
};