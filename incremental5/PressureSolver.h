#include <GL/freeglut.h>

#include <iostream>
#include <math.h>

#include <vector>

#include <cassert>

#include "Array3D.h"
#include "MaterialType.h"

class PressureSolver {
  public:

    PressureSolver()  {
      
    }
    void ProjectPressure(Array3D<MaterialType>& labels, Array3D<unsigned short>& neighbors,
    Array3D<double>& u, Array3D<double>& v,Array3D<double>& w, Array3D<double>* p){
      Array3D<double> r(p->nx(),p->ny(),p->nz());
      Array3D<double> q(p->nx(),p->ny(),p->nz());
      Array3D<double> d(p->nx(),p->ny(),p->nz());
      for(int i = 1; i<p->nx()-1; i++){
        for(int j = 1; j<p->ny()-1; j++){
          for(int k = 1; k<p->nz()-1; k++){
            if (labels(i, j, k) != FLUID) {
              r(i, j, k) = 0.0;
              continue;
            }

            double du_dx = u(i + 1, j, k) - u(i, j, k);
            double dv_dy = v(i, j + 1, k) - v(i, j, k);
            double dw_dz = w(i, j, k + 1) - w(i, j, k);
            double velocity_divergence_of_cell_ijk = du_dx + dv_dy + dw_dz;
            r(i, j, k) = -velocity_divergence_of_cell_ijk;
          }
        }
      }
      d.SetEqualTo(r);
      double sigma = Dot(r, r);
      double tolerance = 1.0e-6 * sigma;
      int thresh = 1000;
      for(int iter = 0; iter < thresh && sigma > tolerance; iter++){
        for(int i = 1; i<p->nx()-1; i++){
          for(int j = 1; j<p->ny()-1; j++){
            for(int k = 1; k<p->nz()-1; k++){
              if (!neighbors(i, j, k)) {
                q(i, j, k) = 0.0;
                continue;
              }

              // Multiply A * d for the row of A corresponding to cell (i, j, k).
              // Store the result in q(i, j, k).

              // Bits & numbers correspond to the matching spot for that calculation
              q(i, j, k) =
                  ((neighbors(i, j, k) & 7) * d(i, j, k)) -
                  ((neighbors(i, j, k) & 8) ? d(i - 1, j, k) : 0) -
                  ((neighbors(i, j, k) & 16 ? d(i, j - 1, k) : 0) -
                  ((neighbors(i, j, k) & 32) ? d(i, j, k - 1) : 0) -
                  ((neighbors(i, j, k) & 64) ? d(i + 1, j, k) : 0) -
                  ((neighbors(i, j, k) & 128) ? d(i, j + 1, k) : 0) -
                  ((neighbors(i, j, k) & 256) ? d(i, j, k + 1) : 0));
            }
          }
        }
        double alpha = sigma / Dot(d,q);
        p->PlusEquals(alpha, d);   // *p += alpha * d_
        r.PlusEquals(-alpha, q);  // r_ -= alpha * q_
        double sigma_old = sigma;
        sigma = Dot(r,r);
        double beta = sigma / sigma_old;
        d.EqualsPlusTimes(r, beta, d);  // d_ = r_ + beta * d_
      }
    }
    /*
    */
    
  private:
   
};