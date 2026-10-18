// Karp's minimum-mean-cycle algorithm for the potentials of the cyclically monotone
// interpolation (used by onlpca()).
//
// Complete graph on the N grid points, edge weights w(u,v) = <g_u, y_u - y_v> (u != v):
//   D_0(source) = 0,  D_k(v) = min_{u != v} [ D_{k-1}(u) + w(u,v) ],   k = 1..N
//   mu*   = min_v max_{0 <= k < N} ( D_N(v) - D_k(v) ) / (N - k)       (maximal margin)
//   psi_v = -min_k [ D_k(v) - k mu* ] + const                          (potentials)
// Weights are computed on the fly (no N x N matrix). With store = true, all D_k are kept
// (N*(N+1) doubles); otherwise the recursion is recomputed in three passes (O(N) memory).
// The loop over v is parallelized with OpenMP.

// OpenMP flags are set in Makevars / Makevars.win.
#include <Rcpp.h>
#include <vector>
#include <limits>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

static void karp_step(const std::vector<double>& prev, std::vector<double>& cur,
                      const std::vector<double>& g, const std::vector<double>& y,
                      const std::vector<double>& a, int N, int d, int nthreads) {
  const double INF = std::numeric_limits<double>::infinity();
  std::vector<double> b(N);
  for (int u = 0; u < N; ++u) b[u] = prev[u] + a[u];          // Inf stays Inf
#ifdef _OPENMP
#pragma omp parallel num_threads(nthreads)
#endif
  {
    std::vector<double> t(N);                                 // per-thread work array
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (int v = 0; v < N; ++v) {
      // t[u] = b[u] - <g_u, y_v>, computed coordinatewise (vectorizable loops)
      for (int u = 0; u < N; ++u) t[u] = b[u];
      for (int m = 0; m < d; ++m) {
        const double yv = y[(size_t)m * N + v];
        const double* gm = &g[(size_t)m * N];
        double* tp = t.data();
#ifdef _OPENMP
#pragma omp simd
#endif
        for (int u = 0; u < N; ++u) tp[u] -= gm[u] * yv;
      }
      double best = INF;
      const double* tp = t.data();
#ifdef _OPENMP
#pragma omp simd reduction(min:best)
#endif
      for (int u = 0; u < v; ++u) best = tp[u] < best ? tp[u] : best;
#ifdef _OPENMP
#pragma omp simd reduction(min:best)
#endif
      for (int u = v + 1; u < N; ++u) best = tp[u] < best ? tp[u] : best;
      cur[v] = best;
    }
  }
}

// [[Rcpp::export]]
List karp_potentials_cpp(NumericMatrix G, NumericMatrix Y, int nthreads = 1, bool store = true) {
  const int N = G.nrow(), d = G.ncol();
  if (Y.nrow() != N || Y.ncol() != d) stop("G and Y must have the same dimensions");
  const double INF = std::numeric_limits<double>::infinity();

  // scaling (as in the R version): divide grid and data by max ||y_i||
  double s = 0;
  for (int i = 0; i < N; ++i) {
    double r = 0;
    for (int m = 0; m < d; ++m) r += Y(i, m) * Y(i, m);
    s = std::max(s, std::sqrt(r));
  }
  if (s <= 0) s = 1;
  std::vector<double> g((size_t)N * d), y((size_t)N * d), a(N, 0.0);
  for (int m = 0; m < d; ++m)
    for (int i = 0; i < N; ++i) {
      g[(size_t)m * N + i] = G(i, m) / s;
      y[(size_t)m * N + i] = Y(i, m) / s;
    }
  for (int i = 0; i < N; ++i)
    for (int m = 0; m < d; ++m) a[i] += g[(size_t)m * N + i] * y[(size_t)m * N + i];

  std::vector<double> prev(N, INF), cur(N, INF);
  std::vector<double> DN(N), maxr(N, -INF), sd(N, INF);
  double mu;

  if (store) {
    std::vector<double> D((size_t)N * (N + 1));                  // D[k*N + v]
    prev[0] = 0;
    for (int v = 0; v < N; ++v) D[v] = prev[v];
    for (int k = 1; k <= N; ++k) {
      karp_step(prev, cur, g, y, a, N, d, nthreads);
      for (int v = 0; v < N; ++v) D[(size_t)k * N + v] = cur[v];
      std::swap(prev, cur);
      if (k % 50 == 0) Rcpp::checkUserInterrupt();
    }
    for (int v = 0; v < N; ++v) DN[v] = D[(size_t)N * N + v];
    for (int k = 0; k < N; ++k)
      for (int v = 0; v < N; ++v) {
        double Dk = D[(size_t)k * N + v];
        if (Dk < INF && DN[v] < INF) maxr[v] = std::max(maxr[v], (DN[v] - Dk) / (N - k));
      }
    mu = INF;
    for (int v = 0; v < N; ++v) if (maxr[v] > -INF) mu = std::min(mu, maxr[v]);
    for (int k = 0; k <= N; ++k)
      for (int v = 0; v < N; ++v) {
        double Dk = D[(size_t)k * N + v];
        if (Dk < INF) sd[v] = std::min(sd[v], Dk - k * mu);
      }
  } else {
    // pass 1: D_N
    std::fill(prev.begin(), prev.end(), INF); prev[0] = 0;
    for (int k = 1; k <= N; ++k) {
      karp_step(prev, cur, g, y, a, N, d, nthreads); std::swap(prev, cur);
      if (k % 50 == 0) Rcpp::checkUserInterrupt();
    }
    DN = prev;
    // pass 2: max_k (D_N - D_k)/(N - k)
    std::fill(prev.begin(), prev.end(), INF); prev[0] = 0;
    for (int k = 0; k < N; ++k) {
      for (int v = 0; v < N; ++v)
        if (prev[v] < INF && DN[v] < INF) maxr[v] = std::max(maxr[v], (DN[v] - prev[v]) / (N - k));
      karp_step(prev, cur, g, y, a, N, d, nthreads); std::swap(prev, cur);
      if (k % 50 == 0) Rcpp::checkUserInterrupt();
    }
    mu = INF;
    for (int v = 0; v < N; ++v) if (maxr[v] > -INF) mu = std::min(mu, maxr[v]);
    // pass 3: shortest distances with weights w - mu
    std::fill(prev.begin(), prev.end(), INF); prev[0] = 0;
    for (int k = 0; k <= N; ++k) {
      for (int v = 0; v < N; ++v) if (prev[v] < INF) sd[v] = std::min(sd[v], prev[v] - k * mu);
      if (k < N) { karp_step(prev, cur, g, y, a, N, d, nthreads); std::swap(prev, cur); }
      if (k % 50 == 0) Rcpp::checkUserInterrupt();
    }
  }

  NumericVector psi(N);
  for (int v = 0; v < N; ++v) psi[v] = (-sd[v] + sd[0]) * s * s;
  return List::create(_["psi"] = psi, _["e0"] = mu * s * s);
}
