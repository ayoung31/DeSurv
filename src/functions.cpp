#include <RcppArmadillo.h>
#include <limits>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

struct CoxStatsFull {
  arma::vec g_z;   // gradient wrt eta: d - haz_c ∘ h_c  (n×1)
  arma::vec w;     // diag Hessian approx: haz_c ∘ h_c   (n×1)
  double loglik;
};

// eta, d must be ordered by nondecreasing time (risk sets are suffixes)
inline CoxStatsFull cox_suffix_stats_full(const arma::vec& eta,
                                          const arma::vec& d) {
  const arma::uword n = eta.n_rows;
  arma::vec m(n), rsk_c(n);

  // suffix max
  m[n-1] = eta[n-1];
  for (arma::sword i = (arma::sword)n-2; i >= 0; --i) m[i] = std::max(m[i+1], eta[i]);

  // centered cumulative risk sums: rsk_c[i] = Σ_{j>=i} exp(eta[j]-m[i])
  rsk_c[n-1] = std::exp(eta[n-1] - m[n-1]);
  for (arma::sword i = (arma::sword)n-2; i >= 0; --i) {
    const double scale = std::exp(m[i+1] - m[i]);
    rsk_c[i] = std::exp(eta[i]-m[i]) + rsk_c[i+1]*scale;
  }

  arma::vec logS0 = arma::log(rsk_c) + m;
  // cum sum of d / rsk_c (centered domain)
  arma::vec h(n);
  h[0] = d[0] * std::exp(-logS0[0]);
  for (arma::uword i = 1; i < n; ++i)
    h[i] = h[i-1] + d[i] * std::exp(-logS0[i]);

  // haz_c = exp(eta - m)
  arma::vec eeta = arma::exp(eta);

  // gradient and diagonal weight
  arma::vec g_z = d - eeta % h;        // n×1
  arma::vec w   = eeta % h;            // n×1  (diag Hessian approx)

  // partial log-likelihood
  double ll = 0.0;
  for (arma::uword i = 0; i < n; ++i)
    ll += d[i] * (eta[i] - logS0[i]);

  return CoxStatsFull{std::move(g_z), std::move(w), ll};
}


// [[Rcpp::export]]
void update_H_cpp(const arma::mat& X,
                  const arma::colvec& y,
                  const arma::colvec& delta,
                  const arma::mat& W,
                  arma::mat& H,
                  int n,
                  int p,
                  int k,
                  double alpha,
                  double lambdaH) {

  double eps=1e-12;

  arma::mat Wt = W.t();
  arma::mat num = Wt * X;
  arma::mat denom = (Wt * W * H) + ((lambdaH*p / k) * H);

  H = H % num / (denom + eps);

  H.transform([&](double v){ return (v < eps ? eps : v); });

  return;
}


 // [[Rcpp::export]]
List calc_loss_cpp(const arma::mat& X,
                   const arma::vec& delta,
                   const arma::mat& W,
                   const arma::mat& H,
                   const arma::vec& beta,
                   int n,
                   int p,
                   int k,
                   int n_event,
                   double alpha, double lambda, double nu,
                   double lambdaW, double lambdaH, arma::rowvec sdZ) {

  double penalty_H = arma::accu(arma::square(H))/(2.0 * n * k);
  double nmf_loss = arma::accu(arma::square(X - W * H))/(2.0 * n * p);

  double penalty_beta = (1 - nu) * arma::accu(arma::square(beta%sdZ.t())) / 2.0 +
                        nu * arma::accu(arma::abs(beta%sdZ.t()));

  arma::vec lp  = X.t() * W * beta;          // linear predictor
  auto cs = cox_suffix_stats_full(lp, delta);

  double surv_loss = cs.loglik * 2.0 / n_event;
  double penalty_W = arma::accu(arma::square(W))/(2.0 * p * k);

  double loss = (1-alpha)*(nmf_loss + lambdaH*penalty_H) -
                alpha * (surv_loss - lambda*penalty_beta) +
                lambdaW*penalty_W;

   // Rcout << "survloss " << surv_loss<< "\n";
   // Rcout << "pen beta " << lambda*penalty_beta << "\n";

  return List::create(
    Named("loss") = loss,
    Named("nmf_loss") = nmf_loss,
    Named("surv_loss") = surv_loss,
    Named("penalty_beta") = penalty_beta,
    Named("penalty_W") = penalty_W,
    Named("penalty_H") = penalty_H
  );
}

struct WUpdateResult {
  bool   accepted;
  double new_loss;
  double theta_used;
  int    backtracks;
};

void update_W_damped_backtrack(const arma::mat& X,
                              const arma::vec& delta,          // n×1, rows time-sorted upstream
                              arma::mat& W,                    // in/out  <-- not const
                              arma::mat& H,
                              arma::vec& beta,
                              int n,
                              int p,
                              int k,
                              int n_event,
                              double alpha,
                              double lambda,
                              double nu,   // nu = elastic-net mixing
                              double lambdaW,
                              double lambdaH,
                              arma::rowvec sdZ,
                              // tuning
                              double theta_init = 1.0,
                              double rho = 0.5,
                              int max_backtracks = 10)
{
  double eps = 1e-12;
  // Cache
  arma::mat Ht = H.t();      // N×k

  // Current loss L0
  List L0list = calc_loss_cpp(X, delta, W, H, beta, n, p, k, n_event,
                              alpha, lambda, nu, lambdaW, lambdaH,sdZ);
  double L0 = as<double>(L0list["loss"]);

  // NMF parts
  arma::mat num   = ((1.0 - alpha) / (n * p)) * (X * Ht);// P×k
  arma::mat denom = ((1.0 - alpha) / (n * p)) * (W * H) * Ht        // P×k
  + (lambdaW / (static_cast<double>(p) * k)) * W;   // L2(W) match loss



  if(alpha > 0.0){

    // Cox pieces
    arma::vec lp = (X.t() * W) * beta;                 // N×1 linear predictor
    auto cs = cox_suffix_stats_full(lp, delta);        // g_z (N×1), loglik

    // Cox gradient wrt W: (Gᵗ * g_z) * βᵗ → P×k
    arma::mat grad_nmf = ((1.0-alpha)/(n*p))* (W*H - X) * H.t();         // P×k
    arma::mat dW_cox = (2.0*alpha/n_event)*(X * cs.g_z) * beta.t(); // P×k

    double gn = arma::norm(grad_nmf, "fro") + 1e-12;
    double gc = arma::norm(dW_cox, "fro") + 1e-12;

    double cox_scale = alpha * std::min(gn / gc, 1e6);
    num = num + cox_scale * dW_cox;

  }

  arma::mat R = num / (denom + eps);
  // R.transform([&](double v){ return (v < eps ? eps : v); });

  // --- Backtracking on θ ---
  const arma::mat W_old = W;
  double theta = theta_init;
  int bt = 0;
  bool ok = false;
  double Lnew = std::numeric_limits<double>::infinity();

  while (bt <= max_backtracks) {
    arma::mat W_trial = W_old % ((1.0 - theta) + theta * R);
    W_trial.transform([&](double v){ return (v < eps ? eps : v); });

    arma::mat  Wn = W_trial;               // normalized W
    arma::mat  Hn = H;                     // H to rescale inversely
    arma::vec  betan = beta;               // beta to rescale to keep eta
    arma::rowvec sdZn = sdZ;

    arma::rowvec cn = arma::sqrt(arma::sum(arma::square(Wn), 0));
    cn.transform([&](double c){ return (c < eps ? 1.0 : c); });
    Wn.each_row() /= cn;
    Hn.each_col() %= cn.t();
    betan %= cn.t();                       // <-- keep eta invariant
    sdZn /= cn;

    List Llist = calc_loss_cpp(X, delta, Wn, Hn, betan, n, p, k, n_event,
                               alpha, lambda, nu, lambdaW, lambdaH, sdZn);
    Lnew = as<double>(Llist["loss"]);

    if (std::isfinite(Lnew) && Lnew <= L0) {
      W = std::move(Wn);
      H = std::move(Hn);
      beta = std::move(betan);
      ok = true;
      break;
    }
    theta *= rho;
    ++bt;
  }

}




// Everything for updating beta below
double lasso(double z, double l1, double l2, double v) {
  double s = 0;

  // Determine the sign of z
  if (z > 0)
    s = 1;
  else if (z < 0)
    s = -1;

  // Apply lasso penalty based on the value of z
  if (std::abs(z) <= l1)
    return 0; // No penalty
  else
    return s * (std::abs(z) - l1) / (v  + l2); // Lasso penalty
}


// X: n×p (rows already sorted by time); d: events (n×1, 0/1)
// lambda: overall penalty; a: starting beta (p×1)
// m: per-feature penalty multipliers (p×1), nu: elastic-net mixing in your API
arma::vec cdfit_cox_dh_one_lambda(const arma::mat& X,
                                  const arma::vec& d,
                                  arma::vec a,
                                  int n,
                                  int p,
                                  int n_event,
                                  double lambda,
                                  double eps,
                                  int max_iter,
                                  const arma::vec& m,
                                  double nu)
{

  arma::vec beta = a;                 // current coefficients
  arma::vec eta  = X * beta;          // linear predictor

  int tot_iter = 0;
  while (tot_iter < max_iter) {
    ++tot_iter;

    // ---- Cox stats at current eta (stable) ----
    auto cs = cox_suffix_stats_full(eta, d);
    const arma::vec& g = cs.g_z;      // gradient wrt eta
    const arma::vec& w = cs.w;        // diag Hessian approx

    double maxChange = 0.0;

    // ---- Coordinate updates (quadratic surrogate around current eta) ----
    int n_beta = beta.n_elem;
    for (int j = 0; j < n_beta; ++j) {
      const arma::vec xj = X.col(j);

      // sufficient stats for j using current (g,w)
      const double xwr = arma::dot(xj, g);                 // = sum_i x_ij * g_i
      const double xwx = arma::dot(w, xj % xj);            // = sum_i w_i * x_ij^2

      // guard against zero curvature
      const double v = std::max(xwx / n_event, 1e-12);           // curvature (per glmnet-style scaling)
      const double u = (xwr / n_event) + v * beta[j];            // unpenalized update point

      // penalties (your API: l1 = λ * m_j * nu, l2 = λ * m_j * (1 - nu))
      const double l1 = lambda * m[j] * nu;
      const double l2 = lambda * m[j] * (1.0 - nu);

      // proximal update (you already have this)
      const double bj_old = beta[j];
      const double bj_new = lasso(u, l1, l2, v);           // soft-thresh with ridge

      beta[j] = bj_new;

      // (Optional) keep eta in sync for diagnostics; g,w remain fixed this outer iter
      const double shift = bj_new - bj_old;
      if (shift != 0.0) {
        eta += shift * xj;                                  // O(n)
        maxChange = std::max(maxChange, std::abs(shift) * std::sqrt(v));
      }
    }

    if (maxChange < eps) break;
    // next loop: recompute cs (g,w) at updated eta
  }

  return beta;
}

// [[Rcpp::export]]
arma::vec update_beta_cpp(const arma::mat& Z,     // Z = X^T W
                          const arma::vec delta,
                          int n, int p, int n_event,
                          double alpha,
                          double nu,
                          double lambda,
                          arma::vec beta0,
                          arma::rowvec& meanZ,
                          arma::rowvec& sdZ,
                          double eps,
                          int max_iter,
                          bool& flag_nan,
                          bool flag_beta
                          ) {

  if(alpha==0.0 && !flag_beta) return beta0;

  // get mean and sd of each covariate
  meanZ = arma::mean(Z, 0);
  sdZ = arma::stddev(Z, 0, 0);

  // find columns with 0 variance
  arma::uvec ns = arma::find(sdZ > 1e-12);
  if (ns.n_elem == 0) {
    // Nothing usable: return zeros or beta0; zeros is safer
    return arma::zeros<arma::vec>(Z.n_cols);
  }

  // subset Z, mean, and sd to those with nonzero variance
  arma::mat Z_keep = Z.cols(ns);
  arma::rowvec meanZ_keep = meanZ.cols(ns);
  arma::rowvec sd_keep = sdZ.cols(ns);

  // center and scale new subsetted Z
  Z_keep.each_row() -= meanZ_keep;
  Z_keep.each_row() /= sd_keep;

  // get beta0 elements with nonzero variance
  arma::vec beta0_keep = beta0.elem(ns);

  // penalty factor of size number of new Z cols
  arma::vec penalty_factor = arma::ones<arma::vec>(Z_keep.n_cols);

  // cdfit_cox_dh_one_lambda MUST be the stabilized version (log-sum-exp), as discussed.
  arma::vec beta_cd = cdfit_cox_dh_one_lambda(
    Z_keep, delta, beta0_keep,
    n, p, n_event,
    lambda, eps, max_iter,
    penalty_factor, nu);

  arma::vec beta_full = arma::zeros<arma::vec>(Z.n_cols);
  arma::vec bb = beta_cd / sd_keep.t();
  beta_full.elem(ns) = bb;

  // sanity checks
  if (!beta_full.is_finite()) {
    flag_nan = true;
    // return beta0 if we get nonfinite solutions
    return beta0;
  }

  return beta_full;  // accept full step

}


// [[Rcpp::export]]
List optimize_loss_cpp(const arma::mat& X_in,
                       const arma::colvec& y_in,
                       const arma::colvec& delta_in,
                       int k,
                       int n,
                       int p,
                       int n_event,
                       const arma::mat& W0,
                       const arma::mat& H0,
                       const arma::colvec& beta0,
                       double alpha,
                       double lambda,
                       double nu,      // nu = elastic-net mixing
                       double lambdaW,
                       double lambdaH,
                       double tol,
                       int maxit,
                       bool verbose,
                       double theta_init,
                       double rho,
                       int max_backtracks,
                       double eps_beta,
                       int max_iter_beta,
                       int max_iter_beta_final)
{

   // ---- Local copies we can reorder ----
   arma::mat   X = X_in;
   arma::colvec y = y_in;
   arma::colvec d = delta_in;

   // Initialize parameters
   arma::mat   W    = W0;
   arma::mat   H    = H0;
   arma::colvec beta = beta0;

   // Sort by increasing time once (apply same permutation to everything)
   arma::uvec tOrder = arma::sort_index(y);
   y = y.elem(tOrder);
   d = d.elem(tOrder);
   X = X.cols(tOrder);
   H = H.cols(tOrder);

   //l2 normalize W and adjust H to correspond
   arma::rowvec cn = arma::sqrt(arma::sum(arma::square(W), 0));
   cn.transform([&](double c){ return (c < 1e-12 ? 1.0 : c); });
   W.each_row() /= cn;
   H.each_col() %= cn.t();

   arma::mat Z = X.t() * W;
   arma::rowvec meanZ = arma::mean(Z, 0);
   arma::rowvec sdZ = arma::stddev(Z, 0, 0);


   // Iteration bookkeeping
   double loss      = 1e-6;
   double loss_prev = loss;
   double eps       = 1.0;
   int    it        = 0;
   bool   flag_nan  = false;

   arma::vec lossit   = arma::zeros<arma::vec>(maxit);
   arma::vec slossit  = arma::zeros<arma::vec>(maxit);
   arma::vec nlossit  = arma::zeros<arma::vec>(maxit);
   arma::vec pblossit = arma::zeros<arma::vec>(maxit);
   arma::vec pwlossit = arma::zeros<arma::vec>(maxit);
   arma::vec phlossit = arma::zeros<arma::vec>(maxit);


   while (eps > tol && it < maxit) {
     loss_prev = loss;
     // ---- H update (your function) ----
     update_H_cpp(X, y, d, W, H, n, p, k, alpha, lambdaH);

     // ---- W update with damping + backtracking ----
     update_W_damped_backtrack(
       X, d, W, H, beta, n, p, k, n_event,
       alpha, lambda, nu, lambdaW, lambdaH,
       sdZ, theta_init, rho,
       max_backtracks
     );

     // ---- β update on Z = (M%X)^T W ----
     {
       Z = X.t() * W;              // N×k (rows already sorted)
       beta = update_beta_cpp(Z, d, n, p, n_event, alpha, nu, lambda, beta, meanZ, sdZ, eps_beta, max_iter_beta, flag_nan, false);
       if (flag_nan) break;
     }

     // ---- Loss & diagnostics ----
     List l = calc_loss_cpp(X, d, W, H, beta, n, p, k, n_event,
                            alpha, lambda, nu, lambdaW, lambdaH, sdZ);
     double new_loss = as<double>(l["loss"]);
     if (!std::isfinite(new_loss)) {
       flag_nan = true;
       //Rcout << "loss is non-finite\n";
       break;
     }
     loss = new_loss;

     double survloss   = as<double>(l["surv_loss"]);
     double nmfloss    = as<double>(l["nmf_loss"]);
     double penaltyW   = as<double>(l["penalty_W"]);
     double penaltyH   = as<double>(l["penalty_H"]);
     double penaltybet = as<double>(l["penalty_beta"]);

     eps = std::abs(loss - loss_prev) / std::max(1e-12, std::abs(loss_prev));

     lossit[it]   = loss;
     slossit[it]  = survloss;
     nlossit[it]  = nmfloss;
     pblossit[it] = penaltybet;
     pwlossit[it] = penaltyW;
     phlossit[it] = penaltyH;

      if (verbose) {
        Rprintf("iter: %d  eps: %.8e  loss: %.8e  nmf: %.8e  surv: %.8e\n",
                it, eps, loss, nmfloss, survloss);
      }
      ++it;

      if (it == maxit) {
        warning("DeSurv fit hit max iterations without convergence");
      }
   }

   // refit beta on final W

   Z = X.t() * W;              // N×k (rows already sorted)
   beta = update_beta_cpp(Z, d, n, p, n_event, alpha, nu, lambda, beta, meanZ, sdZ, eps_beta, max_iter_beta_final, flag_nan, true);

   List loss_final = calc_loss_cpp(X, d, W, H, beta, n, p, k, n_event,
                             alpha, lambda, nu, lambdaW, lambdaH, sdZ);

   return List::create(
     Named("W")           = W,
     Named("H")           = H,
     Named("beta")        = beta,
     Named("meanZ")       = meanZ,
     Named("sdZ")         = sdZ,
     Named("loss")        = loss_final,
     Named("iter")        = it,
     Named("lossit")      = lossit.head(it),
     Named("slossit")     = slossit.head(it),
     Named("nlossit")     = nlossit.head(it),
     Named("pblossit")    = pblossit.head(it),
     Named("pwlossit")    = pwlossit.head(it),
     Named("phlossit")    = phlossit.head(it),
     Named("convergence") = (eps <= tol),
     Named("nan_flag")    = flag_nan
   );
}
