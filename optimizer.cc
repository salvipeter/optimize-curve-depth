#include "optimizer.hh"

#include <Eigen/Dense>

void optimize(BSplineCurve &curve, const VectorVector &directions, size_t n_samples, size_t fix) {
  size_t n = curve.n + 1, p = curve.p, m = n_samples;

  // Boundary constraint setup
  std::vector<size_t> boundary;
  for (size_t i = 0; i < fix; ++i) {
    boundary.push_back(i);
    boundary.push_back(n - i - 1);
  }
  size_t nb = boundary.size();
  auto boundary_weights = DoubleVector(nb, 0);

  double knot_start = curve.knots.front();
  double knot_length = curve.knots.back() - curve.knots.front();

  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n + nb, n + nb);
  Eigen::VectorXd b = Eigen::VectorXd::Zero(n + nb);

  for (size_t i = 0; i < m; ++i) { // sampled data point at u_i
    // Compute the basis function derivatives
    double u = (double)i / (m - 1);
    u = knot_start + u * knot_length;
    size_t span = curve.findSpan(u);
    DoubleMatrix der;
    curve.basisFunctionDerivatives(span, u, 2, der);

    for (size_t j = span - p; j <= span; ++j) // derivative by w_j
      for (size_t k = span - p; k <= span; ++k) { // coefficient of w_k
        A(j, k) += directions[k] * directions[j] * 2 * der[2][k+p-span] * der[2][j+p-span];
        b(j) -= curve.cp[k] * directions[j] * 2 * der[2][k+p-span] * der[2][j+p-span];
      }
  }

  // Boundary constraints
  for (size_t i = 0; i < nb; ++i) {
    A(boundary[i], n + i) = 1;
    A(n + i, boundary[i]) = 1;
    b(n + i) = boundary_weights[i];
  }

  // Solve the system
  Eigen::VectorXd x = A.colPivHouseholderQr().solve(b);
  for (size_t i = 0; i < nb; ++i) // to avoid minimal errors
    x(boundary[i]) = boundary_weights[i];
  for (size_t i = 0; i < n; ++i)
    curve.cp[i] += directions[i] * x(i);
}
