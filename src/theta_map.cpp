List score_cpp(
  arma::mat ipar,
  arma::rowvec resp,
  arma::colvec th,
  const int nd,
  arma::mat sigma,
  const int maxIter = 30,
  const double conv = 0.001,
  bool Fisher = true
) {

  // score input

  int ni = ipar.n_rows;
  arma::colvec new_estimate = th;
  arma::colvec old_estimate;
  int iter = 0;
  bool converged = false;
  arma::mat delta;
  arma::mat abs_delta;
  arma::vec theta;
  arma::vec SE;
  arma::vec d2_diag;
  arma::mat d2_inv;
  bool conv_status;

  // dLL input

  arma::mat sigma_inv(nd, nd);
  arma::rowvec dll(nd, fill::zeros);
  arma::rowvec deriv1(nd);
  arma::mat a = ipar.cols(0, nd - 1);
  arma::colvec d = ipar.col(nd);
  arma::colvec c = ipar.col(nd + 1);
  arma::mat P;
  double num;
  double u;
  arma::rowvec w(nd, fill::zeros);
  bool Bayesian = true;

  // makeFI input

  arma::mat FI(nd, nd, fill::zeros);
  arma::mat deriv2;
  arma::mat FI_temp(nd, nd, fill::zeros);
  arma::mat cf;
  bool addsigma = true;

  // makeHessian input

  arma::mat H(nd, nd, fill::zeros);
  arma::mat H_temp(nd, nd, fill::zeros);

  while ((iter < maxIter) & (converged == false)) {

    iter++;
    old_estimate = new_estimate;

    {

      // calculate first derivative

      dll.zeros();

      if (nd == 1) { sigma_inv = 1; }
      else { sigma_inv = inv(sigma); }

      for (int i=0; i<ni; i++) {
        u = arma::as_scalar(resp.col(i));
        if ((u==1) | (u==0)) {
          P = p_m_3pl(old_estimate, a.row(i), d.row(i), c.row(i));
          num = arma::as_scalar((P-c.row(i))*(u-P)/((1-c.row(i))*P));
          dll = dll + a.row(i)*num;
        }
      }

      if (Bayesian == true) {
        arma::rowvec w(nd, fill::zeros);
        for (int h=0; h<nd; h++) {
          w.zeros();
          w.col(h) = 1;
          dll.col(h) = dll.col(h) - w*sigma_inv*old_estimate;
        }
      }

      deriv1 = dll;

    }

    // calculate second derivative

    if (Fisher == true) {

      FI.zeros();
      FI_temp.zeros();

      for (int i=0; i<ni; i++) {

        P = p_m_3pl(old_estimate, a.row(i), d.row(i), c.row(i));
        cf = (1-P)*pow(P-c.row(i),2.0)/(P*pow(1-c.row(i),2.0));
        FI_temp = trans(a.row(i))*a.row(i);
        num = arma::as_scalar(cf);
        FI_temp = num*FI_temp;
        FI = FI + FI_temp;

      }

      if (addsigma == true) {
        FI = FI + inv(sigma);
      }

      deriv2 = -FI;

    } else {

      H.zeros();
      H_temp.zeros();

      for (int i=0; i<ni; i++) {
        u = arma::as_scalar(resp.col(i));
        if ((u==1) | (u==0)) {
          P = p_m_3pl(old_estimate, a.row(i), d.row(i), c.row(i));
          cf = arma::as_scalar((1-P)*(P-c.row(i))*(c.row(i)*u-pow(P,2.0))/(pow(P,2.0)*pow(1-c.row(i),2.0)));
          H_temp = trans(a.row(i))*a.row(i);
          H_temp = arma::as_scalar(cf)*H_temp;
          H = H + H_temp;
        }
      }

      if (addsigma == true) {
        H = H - inv(sigma);
      }

      deriv2 = H;

    }

    // calculate delta

    delta = inv(deriv2)*trans(deriv1);
    new_estimate = old_estimate - delta;
    abs_delta = abs(delta);
    conv_status = all(vectorise(abs_delta) < conv );
    if (conv_status == true) { converged = true; }

  }

  // makeHessian for deriv2

  H.zeros();
  H_temp.zeros();

  for (int i=0; i<ni; i++) {
    u = arma::as_scalar(resp.col(i));
    if ((u==1) | (u==0)) {
      P = p_m_3pl(new_estimate, a.row(i), d.row(i), c.row(i));
      cf = arma::as_scalar((1-P)*(P-c.row(i))*(c.row(i)*u-pow(P,2.0))/(pow(P,2.0)*pow(1-c.row(i),2.0)));
      H_temp = trans(a.row(i))*a.row(i);
      H_temp = arma::as_scalar(cf)*H_temp;
      H = H + H_temp;
    }
  }

  if (addsigma == true) {
    H = H - inv(sigma);
  }

  deriv2 = H;

  theta = new_estimate;
  d2_inv = inv(deriv2);
  d2_diag = d2_inv.diag();
  SE = pow(abs(d2_diag),0.5);

  return List::create(
    Named("theta")=theta,
    Named("SE")=SE,
    Named("Hessian")=deriv2,
    Named("niter")=iter
  );

}
