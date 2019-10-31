#include <iostream>

#include "optimizer.hh"

int main(int argc, char **argv) {
  BSplineCurve curve;
  curve.p = 3;
  curve.n = 5;
  curve.knots = { 0, 0, 0, 0, 15.0/40.0, 20.0/40.0, 1, 1, 1, 1 };
  curve.cp = {
    { 0, 10,  0},
    {10, 10,  0},
    {15,  0,  0},
    {20,  0,  0},
    {30, 10,  0},
    {40, 10,  0}
  };
  VectorVector dirs(6, {1, 2, 0});
  optimize(curve, dirs, 100, 2);
  for (const auto &cp : curve.cp)
    std::cout << cp.x << ' ' << cp.y << ' ' << cp.z << std::endl;
}
