#include "optimizer.hh"

#include <Eigen/Dense>

void optimize(BSplineCurve &curve, const VectorVector &directions, size_t n_samples, size_t fix) {
  size_t n = curve.n + 1, p = curve.p, m = n_samples;
  size_t nv = n - fix * 2;

  double knot_start = curve.knots.front();
  double knot_length = curve.knots.back() - curve.knots.front();

  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(m * nv, nv);
  Eigen::VectorXd b = Eigen::VectorXd::Zero(m * nv);

  for (size_t i = 0; i < m; ++i) { // sampled data point at u_i
    // Compute the basis function derivatives
    double u = (double)i / (m - 1);
    u = knot_start + u * knot_length;
    size_t span = curve.findSpan(u);
    DoubleMatrix der;
    curve.basisFunctionDerivatives(span, u, 2, der);

    for (size_t j = fix; j < n - fix; ++j) { // derivative by w_j
      size_t jv = j - fix;

      // Outside of the span?
      if (span - p > j || j > span)
        continue;

      auto vj = directions[j];
      for (size_t k = 0; k < n; ++k) // coefficient of w_k
        if (span - p <= k && k <= span) {
          if (k >= fix && k < n - fix)
            A(i * nv + jv, k - fix) = der[2][k+p-span] * der[2][j+p-span];
          b(i * nv + jv) -= curve.cp[k] * vj * der[2][k+p-span] * der[2][j+p-span];
        }
    }
  }

  // Solve the system
  Eigen::VectorXd x = A.colPivHouseholderQr().solve(b);
  for (size_t i = fix; i < n - fix; ++i)
    curve.cp[i] += directions[i] * x(i - fix);
}
