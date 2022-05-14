#include <GL/freeglut.h>

#include <iostream>
#include <math.h>

#include <vector>

#include <cassert>

template <class T> class Array3D {
  public:
    // dims
    GLint x;
    GLint y;
    GLint z;
  
    T* data;

    inline Array3D(int gx, int gy, int gz) : x(gx), y(gy), z(gz), data(new T[gx*gy*gz]) {
      // operator=(0.0);
    }

    // Returns the element at index (|i|, |j|, |k|) of this array.
    inline const T& operator()(std::size_t i, std::size_t j, std::size_t k) const{
      return data[i*y*z + j*z + k];
    }

    // Returns a modifiable reference to the element at index (|i|, |j|,
    // |k|) of this array.
    inline T& operator()(std::size_t i, std::size_t j, std::size_t k){
      return data[i*y*z + j*z + k];
    }

    // assigns all vals in data to val
    inline const T& operator=(const T& value){
      T* data_pointer = data_;
      for (std::size_t i = 0; i < x*y*z; i++) {
        (*data_pointer) = value;
        data_pointer++;
      }
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

      // Adds |scalar| * |arr| to |*this|.
    // |arr| and |*this| must have identical dimensions.
    inline void PlusEquals(double scalar, const Array3D<T>& arr){
      assert(nx()==arr.nx());
      assert(ny()==arr.ny());
      assert(nz()==arr.nz());
      for(int i = 0; i< nx(); i++){
        for(int j = 0; j< ny(); j++){
          for(int k = 0; k< nz(); k++){
            data[i*y*z + j*z + k] += scalar*arr(i,j,k);
          }
        }
      }
    }

    // Sets |*this| = |arr1| + |scalar| * |arr2|.
    //
    // |arr1|, |arr2|, and |*this| must have identical dimensions.
    //
    // It's okay if |arr1| and/or |arr2| and/or |*this| are the same array since
    // each element will be modified one at a time, not affecting any other
    // element.
    inline void EqualsPlusTimes(const Array3D<T>& arr1, double scalar, 
    const Array3D<T>& arr2){
      assert(nx()==arr1.nx());
      assert(ny()==arr1.ny());
      assert(nz()==arr1.nz());
      assert(nx()==arr2.nx());
      assert(ny()==arr2.ny());
      assert(nz()==arr2.nz());
      for(size_t i = 0; i< nx(); i++){
        for(size_t j = 0; j< ny(); j++){
          for(size_t k = 0; k< nz(); k++){
            data[i*y*z + j*z + k] = arr1(i,j,k) + scalar*arr2(i,j,k);
          }
        }
      }
    }

    // Sets |*this| = |other|.
    // We use this function instead of operator= to prevent any automatic calls to
    // operator= that may occur e.g. when returning an Array3D from a function and
    // assigning the returned Array3D to an already instantiated Array3D object.
    inline void SetEqualTo(const Array3D& other){
      assert(nx()==other.nx());
      assert(ny()==other.ny());
      assert(nz()==other.nz());
      for(int i = 0; i< nx(); i++){
        for(int j = 0; j< ny(); j++){
          for(int k = 0; k< nz(); k++){
            data[i*y*z + j*z + k] += other(i,j,k);
          }
        }
      }
    }
  private:
    // Don't allow copy constructor to be called.
    Array3D(const Array3D& other);

    // Don't allow copy-assignment operator to be called.
    Array3D& operator=(const Array3D& other);
};

// Returns the element-wise "dot product" of |a1| and |a2|.
// |a1| and |a2| must have identical dimensions.
inline double Dot(const Array3D<double>& a1, const Array3D<double>& a2) {
  assert(a1.nx()==a2.nx());
  assert(a1.ny()==a2.ny());
  assert(a1.nz()==a2.nz());
  double ret = 0.0;
  for(int i = 0; i< a1.nx(); i++){
    for(int j = 0; j< a1.ny(); j++){
      for(int k = 0; k< a1.nz(); k++){
        ret += a1(i,j,k)*a2(i,j,k);
      }
    }
  }
  return ret;
}