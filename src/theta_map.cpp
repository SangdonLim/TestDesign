#include "item_probability.h"
#include "item_information.h"
#include "jacobian.h"
#include "hessian.h"

//' @rdname map
//' @export
// [[Rcpp::export]]
List estimate_theta_map(
  const arma::mat& ipar,
  const arma::irowvec& item_model,
  const int& nd,
  const arma::rowvec& response,
  const arma::rowvec& start_theta,
  const arma::mat& sigma,
  const int& max_iteration = 30,
  const double& convergence_criterion = 0.001,
  const bool& use_Fisher = true
) {

  // score input

  int ni = ipar.n_rows;
  arma::rowvec new_estimate = start_theta;
  arma::rowvec old_estimate;
  int iteration = 0;
  bool converged = false;
  arma::mat delta;
  arma::mat abs_delta;
  arma::rowvec theta;
  arma::vec se;
  arma::vec d2_diag;
  arma::mat d2_inv;

  // dLL input

  arma::mat sigma_inv(nd, nd);
  arma::rowvec dll(nd, fill::zeros);
  arma::rowvec deriv1(nd);
  arma::rowvec j(nd);
  arma::mat a = ipar.cols(0, nd - 1);
  arma::colvec d = ipar.col(nd);
  arma::colvec c = ipar.col(nd + 1);
  arma::mat dd = ipar.cols(nd, nd + 1);
  arma::mat P;
  arma::rowvec w(nd, fill::zeros);
  bool Bayesian = true;

  // makeFI input

  arma::mat FI(nd, nd, fill::zeros);
  arma::mat deriv2;
  bool addsigma = true;

  // makeHessian input

  arma::mat H(nd, nd, fill::zeros);

  while ((iteration < max_iteration) & (converged == false)) {

    iteration++;
    old_estimate = new_estimate;

    {

      // calculate first derivative

      dll.zeros();

      if (nd == 1) { sigma_inv = 1; }
      else { sigma_inv = inv(sigma); }

      for (int i=0; i<ni; i++) {
        switch (item_model(i)) {
          case 102: {
            dll += j_m_2pl(old_estimate, a.row(i), d(i), response(i));
          }
          break;
          case 103: {
            dll += j_m_3pl(old_estimate, a.row(i), d(i), c(i), response(i));
          }
          break;
          case 105: {
            dll += j_m_gpc(old_estimate, a.row(i), dd.row(i), response(i));
          }
          break;
          case 106: {
            dll += j_m_gr(old_estimate, a.row(i), dd.row(i), response(i));
          }
          break;
        }
      }

      if (Bayesian == true) {
        arma::rowvec w(nd, fill::zeros);
        for (int h=0; h<nd; h++) {
          w.zeros();
          w(h) = 1;
          dll(h) = dll(h) - arma::as_scalar(w * sigma_inv * old_estimate.t());
        }
      }

      deriv1 = dll;

    }

    // calculate second derivative

    if (use_Fisher == true) {

      FI.zeros();

      for (int i=0; i<ni; i++) {
        switch (item_model(i)) {
          case 102: {
            FI = FI + info_m_2pl(old_estimate, a.row(i), d(i));
          }
          break;
          case 103: {
            FI = FI + info_m_3pl(old_estimate, a.row(i), d(i), c(i));
          }
          break;
          case 105: {
            FI = FI + info_m_gpc(old_estimate, a.row(i), dd.row(i));
          }
          break;
          case 106: {
            FI = FI + info_m_gr(old_estimate, a.row(i), dd.row(i));
          }
          break;
        }
      }

      if (addsigma == true) {
        FI = FI + inv(sigma);
      }

      deriv2 = -FI;

    } else {

      H.zeros();

      for (int i=0; i<ni; i++) {
        switch (item_model(i)) {
          case 102: {
            H = H + h_m_2pl(old_estimate, a.row(i), d(i), response(i));
          }
          break;
          case 103: {
            H = H + h_m_3pl(old_estimate, a.row(i), d(i), c(i), response(i));
          }
          break;
          case 105: {
            H = H + h_m_gpc(old_estimate, a.row(i), dd.row(i), response(i));
          }
          break;
          case 106: {
            H = H + h_m_gr(old_estimate, a.row(i), dd.row(i), response(i));
          }
          break;
        }
      }

      if (addsigma == true) {
        H = H - inv(sigma);
      }

      deriv2 = H;

    }

    // calculate delta

    delta = (inv(deriv2) * (deriv1.t())).t();
    new_estimate = old_estimate - delta;
    abs_delta = abs(delta);
    converged = all(vectorise(abs_delta) < convergence_criterion);

  }

  // makeHessian for deriv2

  H.zeros();

  for (int i=0; i<ni; i++) {
    switch (item_model(i)) {
      case 102: {
        H = H + h_m_2pl(new_estimate, a.row(i), d(i), response(i));
      }
      break;
      case 103: {
        H = H + h_m_3pl(new_estimate, a.row(i), d(i), c(i), response(i));
      }
      break;
      case 105: {
        H = H + h_m_gpc(new_estimate, a.row(i), dd.row(i), response(i));
      }
      break;
      case 106: {
        H = H + h_m_gr(new_estimate, a.row(i), dd.row(i), response(i));
      }
      break;
    }
  }

  if (addsigma == true) {
    H = H - inv(sigma);
  }

  deriv2 = H;

  theta = new_estimate;
  d2_inv = inv(deriv2);
  d2_diag = d2_inv.diag();
  se = pow(abs(d2_diag),0.5);

  return List::create(
    Named("theta")=theta,
    Named("se")=se,
    Named("hessian")=deriv2,
    Named("iteration")=iteration
  );

}
