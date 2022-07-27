#ifndef _INFO_ARRAY_FUNCTIONS_H
#define _INFO_ARRAY_FUNCTIONS_H

arma::colvec calc_info(
  const arma::rowvec&,
  const arma::mat&,
  const arma::icolvec&,
  const arma::icolvec&);

arma::colvec calc_thisdirinfo(
  const arma::rowvec&,
  const arma::mat&,
  const int&,
  const arma::icolvec&,
  const arma::icolvec&,
  const arma::rowvec&);

arma::mat calc_info_array(
  const arma::mat&,
  const arma::mat&,
  const arma::icolvec&,
  const arma::icolvec&);

arma::colvec calc_info_EB (
  const arma::mat&,
  const arma::mat&,
  const arma::icolvec&,
  const arma::icolvec&);

arma::colvec calc_info_FB (
  const arma::mat&,
  const List&,
  const arma::icolvec&,
  const arma::icolvec&,
  const bool&);

arma::colvec calc_MI_FB (
  const arma::rowvec&,
  const List&,
  const arma::icolvec&,
  const arma::icolvec&);

#endif
