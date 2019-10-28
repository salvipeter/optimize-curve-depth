#include "optimizer.hh"

#include <Eigen/Dense>

void optimize(BSplineCurve &curve, const VectorVector &directions, size_t n_samples) {
  size_t n = curve.n + 1, p = curve.p, m = n_samples;

  // Boundary constraint setup
  std::vector<size_t> boundary = { 0, n - 1 };
  size_t nb = boundary.size();
  auto boundary_weights = DoubleVector(nb, 0);

  double knot_start = curve.knots.front();
  double knot_length = curve.knots.back() - curve.knots.front();

  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(m * n + nb, n + nb);
  Eigen::VectorXd b = Eigen::VectorXd::Zero(m * n + nb);

  for (size_t i = 0; i < m; ++i) { // sampled data point at u_i
    // Compute the basis function derivatives
    double u = (double)i / (m - 1);
    u = knot_start + u * knot_length;
    size_t span = curve.findSpan(u);
    DoubleMatrix der;
    curve.basisFunctionDerivatives(span, u, 2, der);

    for (size_t j = 0; j < n; ++j) { // derivative by w_j
      // Lagrange multipliers
      auto iter = std::find(boundary.begin(), boundary.end(), j);
      if (iter != boundary.end()) {
        size_t k = iter - boundary.begin();
        A(i * n + j, n + k) = 1;
      }

      // Outside of the span?
      if (span - p > j || j > span)
        continue;

      auto vj = directions[j];
      for (size_t k = 0; k < n + nb; ++k) // coefficient of w_k
        if (span - p <= k && k <= span) {
          A(i * n + j, k) = der[2][k+p-span] * der[2][j+p-span];
          b(i * n + j) -= curve.cp[k] * vj * der[2][k+p-span] * der[2][j+p-span];
        }
    }
  }

  // Boundary constraints
  for (size_t i = 0; i < nb; ++i) {
    A(m * n + i, boundary[i]) = 1;
    b(m * n + i) = boundary_weights[i];
  }

  // Solve the system
  Eigen::VectorXd x = A.colPivHouseholderQr().solve(b);
  for (size_t i = 0; i < nb; ++i) // to avoid minimal errors
    x(boundary[i]) = boundary_weights[i];
  for (size_t i = 0; i < n; ++i)
    curve.cp[i] += directions[i] * x(i);
}
