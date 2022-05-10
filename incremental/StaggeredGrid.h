#include <GL/freeglut.h>

#include <iostream>
#include <math.h>

#include <vector>

#include "Array3D.h"

class StaggeredGrid {
  public:
    // dims
    GLint sx = 0.0;
    GLint sy = 0.0;
    GLint sz = 0.0;

    // p, vel
    Array3D<double>* sp;
    Array3D<double>* su;
    Array3D<double>* sv;
    Array3D<double>* sw;

    StaggeredGrid(int gx, int gy, int gz) : sx(gx), sy(gy), sz(gz), 
    sp(new Array3D<double>(gx,gy,gz)), su(new Array3D<double>(gx+1,gy,gz)), 
    sv(new Array3D<double>(gx,gy+1,gz)), sw(new Array3D<double>(gx,gy,gz+1)) {}

    // Returns the element at index (|i|, |j|, |k|) of this array.
    // inline const T& operator()(std::size_t i, std::size_t j, std::size_t k) const{
    //   return data[i*y*z + j*z + k];
    // }

    // Returns a modifiable reference to the element at index (|i|, |j|,
    // |k|) of this array.
    // inline T& operator()(std::size_t i, std::size_t j, std::size_t k){
    //   return data[i*y*z + j*z + k];
    // }

    Array3D<double> p(){
        return *sp;
    }

    Array3D<double> u(){
        return *su;
    }

    Array3D<double> v(){
        return *sv;
    }

    Array3D<double> w(){
        return *sw;
    }

    const GLint nx(){
        return sx;
    }

    const GLint ny(){
        return sy;
    }

    const GLint nz(){
        return sz;
    }
  private:
};