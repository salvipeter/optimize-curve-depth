#include <fstream>
#include <iostream>

#include "optimizer.hh"
#include "powell.hh"

double curvature(const BSplineCurve &c, double u) {
  VectorVector der;
  c.derivatives(u, 2, der);
  return (der[1] ^ der[2]).norm() / std::pow(der[1].norm(), 3);
}

double curvature_error(const BSplineCurve &c, size_t resolution) {
  double error = 0;
  for (size_t i = 0; i < resolution; ++i) {
    double u = (double)i / (resolution - 1);
    error += std::pow(curvature(c, u), 2);
  }
  return error;
}

int main(int argc, char **argv) {
  BSplineCurve curve;
  Vector dir;

  size_t fix = 2;
  size_t example = 2;
  bool exact_curvature = true;

  switch(example) {
  case 1: // Example 1 - line
    curve.p = 3;
    curve.n = 4;
    curve.knots = { 0, 0, 0, 0, 0.5, 1, 1, 1, 1 };
    curve.cp = {
                {-4, 2, 0},
                {-2, 2, 0},
                {0, 0, 0},
                {2, 2, 0},
                {4, 2, 0}
    };
    dir = {1, 2, 0};
    break;
  case 2: // Example 2 - approximate circle
    curve.p = 3;
    curve.n = 6;
    curve.knots = { 0, 0, 0, 0, 0.25, 0.5, 0.75, 1, 1, 1, 1 };
    curve.cp = {
                {-4, 0,  0},
                {-4, 2,  0},
                {-2, 0,  0},
                { 0, 0,  0},
                { 2, 0,  0},
                { 4, 2,  0},
                { 4, 0,  0}
    };
    dir = {0, 1, 0};
    break;
  case 3: // Example 3 - same, but with slanted lines
    curve.p = 3;
    curve.n = 6;
    curve.knots = { 0, 0, 0, 0, 0.25, 0.5, 0.75, 1, 1, 1, 1 };
    curve.cp = {
                {-4, 0,  0},
                {-4, 2,  0},
                {-2, 0,  0},
                { 0, 0,  0},
                { 2, 0,  0},
                { 4, 2,  0},
                { 4, 0,  0}
    };
    dir = {1, 2, 0};
    break;
  }

  VectorVector dirs(7, dir);

  if (!exact_curvature)
    optimize(curve, dirs, 100, fix);
  else {
    auto minimizer =
      [&](const std::vector<double> &w) {
        auto c = curve;
        for (size_t i = 0; i < w.size(); ++i)
          c.cp[i+fix] += dirs[i+fix] * w[i];
        return curvature_error(c, 100);
      };
    std::vector<double> x(curve.n + 1 - 2 * fix, 0);

    Powell::optimize(minimizer, x, 1000, 1.0e-8, 100, 1.0e-4);

    for (size_t i = 0; i < x.size(); ++i)
      curve.cp[i+fix] += dirs[i+fix] * x[i];
  }

  // Print the results
  for (const auto &cp : curve.cp)
    std::cout << cp.x << ' ' << cp.y << ' ' << cp.z << std::endl;

  // Also save the results in a format loadable in this applet:
  // http://math.bme.hu/~kovacsi/cagd-applets/Curves.html
  std::ofstream f("/tmp/curve.txt");
  f << "TYPE B-spline" << std::endl;
  f << "P " << curve.p << std::endl;
  for (double knot : curve.knots)
    f << "KNOT " << knot << std::endl;
  for (auto p : curve.cp)
    f << "CP " << p.x * 100 << ' ' << p.y * 100 << ' ' << p.z * 100 << std::endl;
}
