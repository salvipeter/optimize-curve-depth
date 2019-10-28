#include <iostream>

#include "optimizer.hh"

int main(int argc, char **argv) {
  BSplineCurve curve;
  curve.p = 3;
  curve.n = 4;
  curve.knots = { 0, 0, 0, 0, 0.5, 1, 1, 1, 1 };
  curve.cp = { {0.3239,   -1.54945,  0},
               {0.712463, -0.347985, 0},
               {1.84883,   0.421245, 0},
               {3.71833,  -0.428571, 0},
               {4.40015,  -1.41758,  1}
  };
  VectorVector dirs(5, {0, 0, 1});
  optimize(curve, dirs, 100);
  for (const auto &cp : curve.cp)
    std::cout << cp.x << ' ' << cp.y << ' ' << cp.z << std::endl;
}
