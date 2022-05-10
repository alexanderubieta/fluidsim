#include <Eigen/Dense>
#include <cassert>
#include <iostream>

typedef Eigen::Matrix<std::size_t, 3, 1> Cell;

inline Cell floor(Eigen::Vector3d p, Eigen::Vector3d lc, double dx){
  assert(dx>0.0);
  Eigen::Vector3d vec = (p-lc) / dx;
  assert(vec[0]>=0.0);
  assert(vec[1]>=0.0);
  assert(vec[2]>=0.0);
  return vec.cast<std::size_t>();
}

inline Eigen::Vector3d weights(Eigen::Vector3d p, Eigen::Vector3d lc,Cell& a, double dx){
  return (p-lc) / dx - a.cast<double>();
}

void print(Cell indices, const Eigen::Vector3d& p, const Eigen::Vector3d& lc, double dx) {
  std::cout << "Indices: " << std::endl;
  std::cout << indices << std::endl;
  Eigen::Vector3d w = weights(p, lc, indices, dx);
  std::cout << "Weights = " << std::endl;
  std::cout << w << std::endl << std::endl;
}

int main(int argc, char** argv) {
  // Make a grid with spacing of 2 with lower corner (1, 2, 3).
  double dx = 2.0;
  Eigen::Vector3d lc(3, 4, 5);

  // (i, j, k) should be (floor((6-1)/2), floor((4-2)/2), floor((3-3)/2)),
  // which is (floor(2.5), floor(1), floor(0)) = (2, 1, 0).
  Eigen::Vector3d p1(4, 10, 14);
  Cell p1_indices = floor(p1, lc, dx);
  print(p1_indices, p1, lc, dx);

  // (i, j, k) should be (3, 8, 9).
  Eigen::Vector3d p2(10, 13, 37);
  Cell p2_indices = floor(p2, lc, dx);
  print(p2_indices, p2, lc, dx);

  // Should lead to an assertion failure for index 0 (x-coordinate) being
  // less than the x-coordinate of |lc|. The program should crash before
  // it even gets to the Print command below.
  Eigen::Vector3d p3(3, 4, 5);
  Cell p3_indices = floor(p3, lc, dx);
  print(p3_indices, p3, lc, dx);

  return 0;
}