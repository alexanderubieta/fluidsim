#include <GL/freeglut.h>

#include <iostream>
#include <math.h>

#include <vector>

template <class T> class Array3D {
  public:
    // dims
    GLint x = 0.0;
    GLint y = 0.0;
    GLint z = 0.0;
  
    T* data;

    Array3D(int gx, int gy, int gz) : x(gx), y(gy), z(gz), data(new T[gx*gy*gz]) {}

    // Returns the element at index (|i|, |j|, |k|) of this array.
    inline const T& operator()(std::size_t i, std::size_t j, std::size_t k) const{
      return data[i*y*z + j*z + k];
    }

    // Returns a modifiable reference to the element at index (|i|, |j|,
    // |k|) of this array.
    inline T& operator()(std::size_t i, std::size_t j, std::size_t k){
      return data[i*y*z + j*z + k];
    }

    inline const T& operator=(const T& value){
      for(int i = 0; i < x*y*z; i++){
        data[i] = value;
      }
      return value;
    }

    GLint nx() const{
        return x;
    }

    GLint ny() const{
        return y;
    }

    GLint nz() const{
        return z;
    }
  private:
    // Don't allow copy constructor to be called.
    Array3D(const Array3D& other);

    // Don't allow copy-assignment operator to be called.
    Array3D& operator=(const Array3D& other);
};