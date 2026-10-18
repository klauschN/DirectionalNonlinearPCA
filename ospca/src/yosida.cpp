// Yosida (Moreau-envelope) regularization of the piecewise-constant cyclically monotone
// interpolation (used by onlpca(); kappa = 0 gives the piecewise-constant interpolation).
//
// phi(z) = max_j (<z, y_j> - psi_j).  For kappa > 0,
//   Qbar(z) = grad phi_kappa(z) = sum_j lambda_j y_j,
//   lambda  = argmin_{lambda in simplex} (kappa/2) ||sum_j lambda_j y_j||^2 - sum_j lambda_j c_j,
//   c_j     = <z, y_j> - psi_j.
// Only pieces with c_j >= max_j c_j - 2 kappa R^2 (R = max ||y_j||) can be active (exact bound).
// The QP is solved by a primal active-set method: at the solution at most d+1 pieces are
// active, so each step only solves a (|S|+1) x (|S|+1) KKT system. Parallel over points.

// OpenMP flags are set in Makevars / Makevars.win.
#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

// Gaussian elimination with partial pivoting; A is n x n row-major, b length n (overwritten)
static bool solve_dense(std::vector<double>& A, std::vector<double>& b, int n) {
  for (int col = 0; col < n; ++col) {
    int piv = col; double best = std::fabs(A[col * n + col]);
    for (int r = col + 1; r < n; ++r) {
      double v = std::fabs(A[r * n + col]);
      if (v > best) { best = v; piv = r; }
    }
    if (best < 1e-300) return false;
    if (piv != col) {
      for (int k = 0; k < n; ++k) std::swap(A[col * n + k], A[piv * n + k]);
      std::swap(b[col], b[piv]);
    }
    for (int r = col + 1; r < n; ++r) {
      double f = A[r * n + col] / A[col * n + col];
      if (f == 0) continue;
      for (int k = col; k < n; ++k) A[r * n + k] -= f * A[col * n + k];
      b[r] -= f * b[col];
    }
  }
  for (int r = n - 1; r >= 0; --r) {
    double s = b[r];
    for (int k = r + 1; k < n; ++k) s -= A[r * n + k] * b[k];
    b[r] = s / A[r * n + r];
  }
  return true;
}

// [[Rcpp::export]]
NumericMatrix yosida_cpp(NumericMatrix Z, NumericMatrix Y, NumericVector psi, double kappa,
                         int nthreads = 1, int max_iter = 500) {
  const int N = Z.nrow(), n = Y.nrow(), d = Y.ncol();
  if (Z.ncol() != d || psi.size() != n) stop("dimension mismatch");
  const double INF = std::numeric_limits<double>::infinity();
  std::vector<double> y((size_t)n * d), z((size_t)N * d), ps(n);
  double R2 = 0;
  for (int j = 0; j < n; ++j) {
    double r = 0;
    for (int m = 0; m < d; ++m) { y[(size_t)j * d + m] = Y(j, m); r += Y(j, m) * Y(j, m); }
    R2 = std::max(R2, r); ps[j] = psi[j];
  }
  for (int i = 0; i < N; ++i) for (int m = 0; m < d; ++m) z[(size_t)i * d + m] = Z(i, m);
  std::vector<double> out((size_t)N * d);
  int n_fail = 0;                                             // points where max_iter was reached

#ifdef _OPENMP
#pragma omp parallel num_threads(nthreads)
#endif
  {
    std::vector<double> s(n);
    std::vector<int> cand; std::vector<double> c, g, lam, w(d);
    std::vector<int> S; std::vector<char> inS;
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 64) reduction(+:n_fail)
#endif
    for (int i = 0; i < N; ++i) {
      const double* zi = &z[(size_t)i * d];
      // scores and candidate set
      double smax = -INF; int jmax = 0;
      for (int j = 0; j < n; ++j) {
        double v = -ps[j];
        for (int m = 0; m < d; ++m) v += zi[m] * y[(size_t)j * d + m];
        s[j] = v; if (v > smax) { smax = v; jmax = j; }
      }
      double* oi = &out[(size_t)i * d];
      if (kappa <= 0) { for (int m = 0; m < d; ++m) oi[m] = y[(size_t)jmax * d + m]; continue; }
      const double thr = smax - 2.0 * kappa * R2;
      cand.clear(); c.clear();
      int start = 0;
      for (int j = 0; j < n; ++j) if (s[j] >= thr) {
        if (j == jmax) start = (int)cand.size();
        cand.push_back(j); c.push_back(s[j]);
      }
      const int K = (int)cand.size();
      if (K == 1) { for (int m = 0; m < d; ++m) oi[m] = y[(size_t)jmax * d + m]; continue; }
      auto yc = [&](int a, int m) { return y[(size_t)cand[a] * d + m]; };

      // active-set method, starting from the unsmoothed answer
      lam.assign(K, 0.0); inS.assign(K, 0); S.clear();
      S.push_back(start); inS[start] = 1; lam[start] = 1.0;
      g.resize(K);
      double cscale = 1.0; for (int a = 0; a < K; ++a) cscale = std::max(cscale, std::fabs(c[a]));
      const double tol = 1e-12 * cscale;
      bool conv = false;
      for (int it = 0; it < max_iter; ++it) {
        std::fill(w.begin(), w.end(), 0.0);
        for (int a : S) for (int m = 0; m < d; ++m) w[m] += lam[a] * yc(a, m);
        double nu = 0;
        for (int a = 0; a < K; ++a) {
          double v = 0;
          for (int m = 0; m < d; ++m) v += yc(a, m) * w[m];
          g[a] = kappa * v - c[a];
        }
        for (int a : S) nu += lam[a] * g[a];
        int jn = -1; double gmin = INF;
        for (int a = 0; a < K; ++a) if (!inS[a] && g[a] < gmin) { gmin = g[a]; jn = a; }
        if (jn < 0 || gmin >= nu - tol) { conv = true; break; }    // optimal
        S.push_back(jn); inS[jn] = 1; lam[jn] = 0.0;
        for (int inner = 0; inner < 50; ++inner) {
          const int p = (int)S.size(), q = p + 1;
          std::vector<double> A((size_t)q * q, 0.0), b(q, 0.0);
          double hmax = 0;
          for (int r = 0; r < p; ++r) for (int t = 0; t < p; ++t) {
            double v = 0;
            for (int m = 0; m < d; ++m) v += yc(S[r], m) * yc(S[t], m);
            A[r * q + t] = kappa * v; if (r == t) hmax = std::max(hmax, kappa * v);
          }
          for (int r = 0; r < p; ++r) { A[r * q + r] += 1e-12 * std::max(1.0, hmax);
                                        A[r * q + p] = 1.0; A[p * q + r] = 1.0; b[r] = c[S[r]]; }
          b[p] = 1.0;
          if (!solve_dense(A, b, q)) break;
          bool feasible = true;
          for (int r = 0; r < p; ++r) if (b[r] < -1e-14) { feasible = false; break; }
          if (feasible) {
            double sum = 0;
            for (int r = 0; r < p; ++r) { lam[S[r]] = std::max(b[r], 0.0); sum += lam[S[r]]; }
            for (int r = 0; r < p; ++r) lam[S[r]] /= sum;
            break;
          }
          // step towards the subproblem solution until a weight hits zero, then drop it
          double tstep = 1.0;
          for (int r = 0; r < p; ++r) {
            double cur = lam[S[r]], nw = b[r];
            if (nw < cur && nw < 0) tstep = std::min(tstep, cur / (cur - nw));
          }
          std::vector<int> keep;
          for (int r = 0; r < p; ++r) {
            int a = S[r];
            lam[a] = lam[a] + tstep * (b[r] - lam[a]);
            if (lam[a] > 1e-15) keep.push_back(a); else { lam[a] = 0.0; inS[a] = 0; }
          }
          if (keep.empty()) { keep.push_back(start); inS[start] = 1; lam[start] = 1.0; }
          S = keep;
          double sum = 0; for (int a : S) sum += lam[a];
          for (int a : S) lam[a] /= sum;
        }
      }
      std::fill(w.begin(), w.end(), 0.0);
      for (int a : S) for (int m = 0; m < d; ++m) w[m] += lam[a] * yc(a, m);
      for (int m = 0; m < d; ++m) oi[m] = w[m];
      if (!conv) ++n_fail;
    }
  }
  NumericMatrix res(N, d);
  for (int i = 0; i < N; ++i) for (int m = 0; m < d; ++m) res(i, m) = out[(size_t)i * d + m];
  res.attr("not_converged") = n_fail;
  return res;
}
