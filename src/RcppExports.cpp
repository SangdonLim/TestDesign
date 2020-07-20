// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// calc_info
arma::colvec calc_info(const arma::rowvec& x, const arma::mat& item_parm, const arma::icolvec& ncat, const arma::icolvec& model);
RcppExport SEXP _TestDesign_calc_info(SEXP xSEXP, SEXP item_parmSEXP, SEXP ncatSEXP, SEXP modelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::rowvec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type item_parm(item_parmSEXP);
    Rcpp::traits::input_parameter< const arma::icolvec& >::type ncat(ncatSEXP);
    Rcpp::traits::input_parameter< const arma::icolvec& >::type model(modelSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_info(x, item_parm, ncat, model));
    return rcpp_result_gen;
END_RCPP
}
// calc_info_matrix
arma::mat calc_info_matrix(const arma::mat& x, const arma::mat& item_parm, const arma::icolvec& ncat, const arma::icolvec& model);
RcppExport SEXP _TestDesign_calc_info_matrix(SEXP xSEXP, SEXP item_parmSEXP, SEXP ncatSEXP, SEXP modelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type item_parm(item_parmSEXP);
    Rcpp::traits::input_parameter< const arma::icolvec& >::type ncat(ncatSEXP);
    Rcpp::traits::input_parameter< const arma::icolvec& >::type model(modelSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_info_matrix(x, item_parm, ncat, model));
    return rcpp_result_gen;
END_RCPP
}
// calc_info_EB
arma::colvec calc_info_EB(const arma::mat& x, const arma::mat& item_parm, const arma::icolvec& ncat, const arma::icolvec& model);
RcppExport SEXP _TestDesign_calc_info_EB(SEXP xSEXP, SEXP item_parmSEXP, SEXP ncatSEXP, SEXP modelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type item_parm(item_parmSEXP);
    Rcpp::traits::input_parameter< const arma::icolvec& >::type ncat(ncatSEXP);
    Rcpp::traits::input_parameter< const arma::icolvec& >::type model(modelSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_info_EB(x, item_parm, ncat, model));
    return rcpp_result_gen;
END_RCPP
}
// calc_info_FB
arma::colvec calc_info_FB(const arma::mat& x, const List& items_list, const arma::icolvec& ncat, const arma::icolvec& model, const bool& useEAP);
RcppExport SEXP _TestDesign_calc_info_FB(SEXP xSEXP, SEXP items_listSEXP, SEXP ncatSEXP, SEXP modelSEXP, SEXP useEAPSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const List& >::type items_list(items_listSEXP);
    Rcpp::traits::input_parameter< const arma::icolvec& >::type ncat(ncatSEXP);
    Rcpp::traits::input_parameter< const arma::icolvec& >::type model(modelSEXP);
    Rcpp::traits::input_parameter< const bool& >::type useEAP(useEAPSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_info_FB(x, items_list, ncat, model, useEAP));
    return rcpp_result_gen;
END_RCPP
}
// calc_MI_FB
arma::colvec calc_MI_FB(const arma::rowvec& x, const List items_list, const arma::icolvec& ncat, const arma::icolvec& model);
RcppExport SEXP _TestDesign_calc_MI_FB(SEXP xSEXP, SEXP items_listSEXP, SEXP ncatSEXP, SEXP modelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::rowvec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const List >::type items_list(items_listSEXP);
    Rcpp::traits::input_parameter< const arma::icolvec& >::type ncat(ncatSEXP);
    Rcpp::traits::input_parameter< const arma::icolvec& >::type model(modelSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_MI_FB(x, items_list, ncat, model));
    return rcpp_result_gen;
END_RCPP
}
// p_1pl
double p_1pl(const arma::rowvec& x, const double& b);
RcppExport SEXP _TestDesign_p_1pl(SEXP xSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::rowvec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(p_1pl(x, b));
    return rcpp_result_gen;
END_RCPP
}
// p_2pl
double p_2pl(const arma::rowvec& x, const double& a, const double& b);
RcppExport SEXP _TestDesign_p_2pl(SEXP xSEXP, SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::rowvec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double& >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(p_2pl(x, a, b));
    return rcpp_result_gen;
END_RCPP
}
// p_3pl
double p_3pl(const arma::rowvec& x, const double& a, const double& b, const double& c);
RcppExport SEXP _TestDesign_p_3pl(SEXP xSEXP, SEXP aSEXP, SEXP bSEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::rowvec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const double& >::type c(cSEXP);
    rcpp_result_gen = Rcpp::wrap(p_3pl(x, a, b, c));
    return rcpp_result_gen;
END_RCPP
}
// p_pc
arma::rowvec p_pc(const arma::rowvec& x, const arma::rowvec& b);
RcppExport SEXP _TestDesign_p_pc(SEXP xSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::rowvec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(p_pc(x, b));
    return rcpp_result_gen;
END_RCPP
}
// p_gpc
arma::rowvec p_gpc(const arma::rowvec& x, const double& a, const arma::rowvec& b);
RcppExport SEXP _TestDesign_p_gpc(SEXP xSEXP, SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::rowvec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(p_gpc(x, a, b));
    return rcpp_result_gen;
END_RCPP
}
// p_gr
arma::rowvec p_gr(const arma::rowvec& x, const double& a, const arma::rowvec& b);
RcppExport SEXP _TestDesign_p_gr(SEXP xSEXP, SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::rowvec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(p_gr(x, a, b));
    return rcpp_result_gen;
END_RCPP
}
// array_p_1pl
arma::colvec array_p_1pl(const arma::mat& x, const double& b);
RcppExport SEXP _TestDesign_array_p_1pl(SEXP xSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(array_p_1pl(x, b));
    return rcpp_result_gen;
END_RCPP
}
// array_p_2pl
arma::colvec array_p_2pl(const arma::mat& x, const double& a, const double& b);
RcppExport SEXP _TestDesign_array_p_2pl(SEXP xSEXP, SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double& >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(array_p_2pl(x, a, b));
    return rcpp_result_gen;
END_RCPP
}
// array_p_3pl
arma::colvec array_p_3pl(const arma::mat& x, const double& a, const double& b, const double& c);
RcppExport SEXP _TestDesign_array_p_3pl(SEXP xSEXP, SEXP aSEXP, SEXP bSEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const double& >::type c(cSEXP);
    rcpp_result_gen = Rcpp::wrap(array_p_3pl(x, a, b, c));
    return rcpp_result_gen;
END_RCPP
}
// array_p_pc
arma::mat array_p_pc(const arma::mat& x, const arma::rowvec& b);
RcppExport SEXP _TestDesign_array_p_pc(SEXP xSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(array_p_pc(x, b));
    return rcpp_result_gen;
END_RCPP
}
// array_p_gpc
arma::mat array_p_gpc(const arma::mat& x, const double& a, const arma::rowvec& b);
RcppExport SEXP _TestDesign_array_p_gpc(SEXP xSEXP, SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(array_p_gpc(x, a, b));
    return rcpp_result_gen;
END_RCPP
}
// array_p_gr
arma::mat array_p_gr(const arma::mat& x, const double& a, const arma::rowvec& b);
RcppExport SEXP _TestDesign_array_p_gr(SEXP xSEXP, SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(array_p_gr(x, a, b));
    return rcpp_result_gen;
END_RCPP
}
// info_1pl
double info_1pl(const arma::rowvec& x, const double& b);
RcppExport SEXP _TestDesign_info_1pl(SEXP xSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::rowvec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(info_1pl(x, b));
    return rcpp_result_gen;
END_RCPP
}
// info_2pl
double info_2pl(const arma::rowvec& x, const double& a, const double& b);
RcppExport SEXP _TestDesign_info_2pl(SEXP xSEXP, SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::rowvec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double& >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(info_2pl(x, a, b));
    return rcpp_result_gen;
END_RCPP
}
// info_3pl
double info_3pl(const arma::rowvec& x, const double& a, const double& b, const double& c);
RcppExport SEXP _TestDesign_info_3pl(SEXP xSEXP, SEXP aSEXP, SEXP bSEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::rowvec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const double& >::type c(cSEXP);
    rcpp_result_gen = Rcpp::wrap(info_3pl(x, a, b, c));
    return rcpp_result_gen;
END_RCPP
}
// info_pc
double info_pc(const arma::rowvec& x, const arma::rowvec& b);
RcppExport SEXP _TestDesign_info_pc(SEXP xSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::rowvec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(info_pc(x, b));
    return rcpp_result_gen;
END_RCPP
}
// info_gpc
double info_gpc(const arma::rowvec& x, const double& a, const arma::rowvec& b);
RcppExport SEXP _TestDesign_info_gpc(SEXP xSEXP, SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::rowvec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(info_gpc(x, a, b));
    return rcpp_result_gen;
END_RCPP
}
// info_gr
double info_gr(const arma::rowvec& x, const double& a, const arma::rowvec& b);
RcppExport SEXP _TestDesign_info_gr(SEXP xSEXP, SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::rowvec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(info_gr(x, a, b));
    return rcpp_result_gen;
END_RCPP
}
// array_info_1pl
arma::colvec array_info_1pl(const arma::mat& x, const double& b);
RcppExport SEXP _TestDesign_array_info_1pl(SEXP xSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(array_info_1pl(x, b));
    return rcpp_result_gen;
END_RCPP
}
// array_info_2pl
arma::colvec array_info_2pl(const arma::mat& x, const double& a, const double& b);
RcppExport SEXP _TestDesign_array_info_2pl(SEXP xSEXP, SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double& >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(array_info_2pl(x, a, b));
    return rcpp_result_gen;
END_RCPP
}
// array_info_3pl
arma::colvec array_info_3pl(const arma::mat& x, const double& a, const double& b, const double& c);
RcppExport SEXP _TestDesign_array_info_3pl(SEXP xSEXP, SEXP aSEXP, SEXP bSEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const double& >::type c(cSEXP);
    rcpp_result_gen = Rcpp::wrap(array_info_3pl(x, a, b, c));
    return rcpp_result_gen;
END_RCPP
}
// array_info_pc
arma::colvec array_info_pc(const arma::mat& x, const arma::rowvec& b);
RcppExport SEXP _TestDesign_array_info_pc(SEXP xSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(array_info_pc(x, b));
    return rcpp_result_gen;
END_RCPP
}
// array_info_gpc
arma::colvec array_info_gpc(const arma::mat& x, const double& a, const arma::rowvec& b);
RcppExport SEXP _TestDesign_array_info_gpc(SEXP xSEXP, SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(array_info_gpc(x, a, b));
    return rcpp_result_gen;
END_RCPP
}
// array_info_gr
arma::colvec array_info_gr(const arma::mat& x, const double& a, const arma::rowvec& b);
RcppExport SEXP _TestDesign_array_info_gr(SEXP xSEXP, SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(array_info_gr(x, a, b));
    return rcpp_result_gen;
END_RCPP
}
// find_segment
IntegerVector find_segment(NumericVector x, NumericVector segment);
RcppExport SEXP _TestDesign_find_segment(SEXP xSEXP, SEXP segmentSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type segment(segmentSEXP);
    rcpp_result_gen = Rcpp::wrap(find_segment(x, segment));
    return rcpp_result_gen;
END_RCPP
}
// calc_likelihood
double calc_likelihood(const arma::rowvec& x, const arma::mat& item_parm, const arma::icolvec& resp, const arma::icolvec& ncat, const arma::icolvec& model);
RcppExport SEXP _TestDesign_calc_likelihood(SEXP xSEXP, SEXP item_parmSEXP, SEXP respSEXP, SEXP ncatSEXP, SEXP modelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::rowvec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type item_parm(item_parmSEXP);
    Rcpp::traits::input_parameter< const arma::icolvec& >::type resp(respSEXP);
    Rcpp::traits::input_parameter< const arma::icolvec& >::type ncat(ncatSEXP);
    Rcpp::traits::input_parameter< const arma::icolvec& >::type model(modelSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_likelihood(x, item_parm, resp, ncat, model));
    return rcpp_result_gen;
END_RCPP
}
// calc_likelihood_function
arma::colvec calc_likelihood_function(const arma::mat& theta_grid, const arma::mat& item_parm, const arma::icolvec& resp, const arma::icolvec& ncat, const arma::icolvec& model);
RcppExport SEXP _TestDesign_calc_likelihood_function(SEXP theta_gridSEXP, SEXP item_parmSEXP, SEXP respSEXP, SEXP ncatSEXP, SEXP modelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type theta_grid(theta_gridSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type item_parm(item_parmSEXP);
    Rcpp::traits::input_parameter< const arma::icolvec& >::type resp(respSEXP);
    Rcpp::traits::input_parameter< const arma::icolvec& >::type ncat(ncatSEXP);
    Rcpp::traits::input_parameter< const arma::icolvec& >::type model(modelSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_likelihood_function(theta_grid, item_parm, resp, ncat, model));
    return rcpp_result_gen;
END_RCPP
}
// calc_log_likelihood
double calc_log_likelihood(const arma::rowvec& x, const arma::mat& item_parm, const arma::icolvec& resp, const arma::icolvec& ncat, const arma::icolvec& model, const int& prior, const arma::rowvec& prior_parm);
RcppExport SEXP _TestDesign_calc_log_likelihood(SEXP xSEXP, SEXP item_parmSEXP, SEXP respSEXP, SEXP ncatSEXP, SEXP modelSEXP, SEXP priorSEXP, SEXP prior_parmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::rowvec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type item_parm(item_parmSEXP);
    Rcpp::traits::input_parameter< const arma::icolvec& >::type resp(respSEXP);
    Rcpp::traits::input_parameter< const arma::icolvec& >::type ncat(ncatSEXP);
    Rcpp::traits::input_parameter< const arma::icolvec& >::type model(modelSEXP);
    Rcpp::traits::input_parameter< const int& >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type prior_parm(prior_parmSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_log_likelihood(x, item_parm, resp, ncat, model, prior, prior_parm));
    return rcpp_result_gen;
END_RCPP
}
// calc_log_likelihood_function
arma::colvec calc_log_likelihood_function(const arma::mat& theta_grid, const arma::mat& item_parm, const arma::icolvec& resp, const arma::icolvec& ncat, const arma::icolvec& model, const int& prior, const arma::rowvec& prior_parm);
RcppExport SEXP _TestDesign_calc_log_likelihood_function(SEXP theta_gridSEXP, SEXP item_parmSEXP, SEXP respSEXP, SEXP ncatSEXP, SEXP modelSEXP, SEXP priorSEXP, SEXP prior_parmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type theta_grid(theta_gridSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type item_parm(item_parmSEXP);
    Rcpp::traits::input_parameter< const arma::icolvec& >::type resp(respSEXP);
    Rcpp::traits::input_parameter< const arma::icolvec& >::type ncat(ncatSEXP);
    Rcpp::traits::input_parameter< const arma::icolvec& >::type model(modelSEXP);
    Rcpp::traits::input_parameter< const int& >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type prior_parm(prior_parmSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_log_likelihood_function(theta_grid, item_parm, resp, ncat, model, prior, prior_parm));
    return rcpp_result_gen;
END_RCPP
}
// calc_posterior
double calc_posterior(const arma::rowvec& x, const arma::mat& item_parm, const arma::icolvec& resp, const arma::icolvec& ncat, const arma::icolvec& model, const int& prior, const arma::rowvec& prior_parm);
RcppExport SEXP _TestDesign_calc_posterior(SEXP xSEXP, SEXP item_parmSEXP, SEXP respSEXP, SEXP ncatSEXP, SEXP modelSEXP, SEXP priorSEXP, SEXP prior_parmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::rowvec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type item_parm(item_parmSEXP);
    Rcpp::traits::input_parameter< const arma::icolvec& >::type resp(respSEXP);
    Rcpp::traits::input_parameter< const arma::icolvec& >::type ncat(ncatSEXP);
    Rcpp::traits::input_parameter< const arma::icolvec& >::type model(modelSEXP);
    Rcpp::traits::input_parameter< const int& >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type prior_parm(prior_parmSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_posterior(x, item_parm, resp, ncat, model, prior, prior_parm));
    return rcpp_result_gen;
END_RCPP
}
// calc_posterior_function
arma::colvec calc_posterior_function(const arma::mat& theta_grid, const arma::mat& item_parm, const arma::icolvec& resp, const arma::icolvec& ncat, const arma::icolvec& model, const int& prior, const arma::rowvec& prior_parm);
RcppExport SEXP _TestDesign_calc_posterior_function(SEXP theta_gridSEXP, SEXP item_parmSEXP, SEXP respSEXP, SEXP ncatSEXP, SEXP modelSEXP, SEXP priorSEXP, SEXP prior_parmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type theta_grid(theta_gridSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type item_parm(item_parmSEXP);
    Rcpp::traits::input_parameter< const arma::icolvec& >::type resp(respSEXP);
    Rcpp::traits::input_parameter< const arma::icolvec& >::type ncat(ncatSEXP);
    Rcpp::traits::input_parameter< const arma::icolvec& >::type model(modelSEXP);
    Rcpp::traits::input_parameter< const int& >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type prior_parm(prior_parmSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_posterior_function(theta_grid, item_parm, resp, ncat, model, prior, prior_parm));
    return rcpp_result_gen;
END_RCPP
}
// calc_posterior_single
double calc_posterior_single(const arma::rowvec& x, const arma::rowvec& item_parm, const int& resp, const int& ncat, const int& model, const int& prior, const arma::rowvec& prior_parm);
RcppExport SEXP _TestDesign_calc_posterior_single(SEXP xSEXP, SEXP item_parmSEXP, SEXP respSEXP, SEXP ncatSEXP, SEXP modelSEXP, SEXP priorSEXP, SEXP prior_parmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::rowvec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type item_parm(item_parmSEXP);
    Rcpp::traits::input_parameter< const int& >::type resp(respSEXP);
    Rcpp::traits::input_parameter< const int& >::type ncat(ncatSEXP);
    Rcpp::traits::input_parameter< const int& >::type model(modelSEXP);
    Rcpp::traits::input_parameter< const int& >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type prior_parm(prior_parmSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_posterior_single(x, item_parm, resp, ncat, model, prior, prior_parm));
    return rcpp_result_gen;
END_RCPP
}
// theta_EAP
arma::colvec theta_EAP(const arma::mat& theta_grid, const arma::mat& item_parm, const arma::icolvec& resp, const arma::icolvec& ncat, const arma::icolvec& model, const int& prior, const arma::rowvec& prior_parm);
RcppExport SEXP _TestDesign_theta_EAP(SEXP theta_gridSEXP, SEXP item_parmSEXP, SEXP respSEXP, SEXP ncatSEXP, SEXP modelSEXP, SEXP priorSEXP, SEXP prior_parmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type theta_grid(theta_gridSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type item_parm(item_parmSEXP);
    Rcpp::traits::input_parameter< const arma::icolvec& >::type resp(respSEXP);
    Rcpp::traits::input_parameter< const arma::icolvec& >::type ncat(ncatSEXP);
    Rcpp::traits::input_parameter< const arma::icolvec& >::type model(modelSEXP);
    Rcpp::traits::input_parameter< const int& >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type prior_parm(prior_parmSEXP);
    rcpp_result_gen = Rcpp::wrap(theta_EAP(theta_grid, item_parm, resp, ncat, model, prior, prior_parm));
    return rcpp_result_gen;
END_RCPP
}
// theta_EAP_matrix
arma::mat theta_EAP_matrix(const arma::mat& theta_grid, const arma::mat& item_parm, const arma::imat& resp, const arma::icolvec& ncat, const arma::icolvec& model, const int& prior, const arma::rowvec& prior_parm);
RcppExport SEXP _TestDesign_theta_EAP_matrix(SEXP theta_gridSEXP, SEXP item_parmSEXP, SEXP respSEXP, SEXP ncatSEXP, SEXP modelSEXP, SEXP priorSEXP, SEXP prior_parmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type theta_grid(theta_gridSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type item_parm(item_parmSEXP);
    Rcpp::traits::input_parameter< const arma::imat& >::type resp(respSEXP);
    Rcpp::traits::input_parameter< const arma::icolvec& >::type ncat(ncatSEXP);
    Rcpp::traits::input_parameter< const arma::icolvec& >::type model(modelSEXP);
    Rcpp::traits::input_parameter< const int& >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type prior_parm(prior_parmSEXP);
    rcpp_result_gen = Rcpp::wrap(theta_EAP_matrix(theta_grid, item_parm, resp, ncat, model, prior, prior_parm));
    return rcpp_result_gen;
END_RCPP
}
// theta_EB
arma::mat theta_EB(const int& nx, const arma::rowvec& theta_init, const double& theta_prop, const arma::mat& item_parm, const arma::icolvec& resp, const arma::icolvec& ncat, const arma::icolvec& model, const int& prior, const arma::rowvec& prior_parm);
RcppExport SEXP _TestDesign_theta_EB(SEXP nxSEXP, SEXP theta_initSEXP, SEXP theta_propSEXP, SEXP item_parmSEXP, SEXP respSEXP, SEXP ncatSEXP, SEXP modelSEXP, SEXP priorSEXP, SEXP prior_parmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type nx(nxSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type theta_init(theta_initSEXP);
    Rcpp::traits::input_parameter< const double& >::type theta_prop(theta_propSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type item_parm(item_parmSEXP);
    Rcpp::traits::input_parameter< const arma::icolvec& >::type resp(respSEXP);
    Rcpp::traits::input_parameter< const arma::icolvec& >::type ncat(ncatSEXP);
    Rcpp::traits::input_parameter< const arma::icolvec& >::type model(modelSEXP);
    Rcpp::traits::input_parameter< const int& >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type prior_parm(prior_parmSEXP);
    rcpp_result_gen = Rcpp::wrap(theta_EB(nx, theta_init, theta_prop, item_parm, resp, ncat, model, prior, prior_parm));
    return rcpp_result_gen;
END_RCPP
}
// theta_EB_single
arma::mat theta_EB_single(const int& nx, const arma::rowvec& theta_init, const double& theta_prop, const arma::rowvec& item_parm, const int& resp, const int& ncat, const int& model, const int& prior, const arma::rowvec& prior_parm);
RcppExport SEXP _TestDesign_theta_EB_single(SEXP nxSEXP, SEXP theta_initSEXP, SEXP theta_propSEXP, SEXP item_parmSEXP, SEXP respSEXP, SEXP ncatSEXP, SEXP modelSEXP, SEXP priorSEXP, SEXP prior_parmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type nx(nxSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type theta_init(theta_initSEXP);
    Rcpp::traits::input_parameter< const double& >::type theta_prop(theta_propSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type item_parm(item_parmSEXP);
    Rcpp::traits::input_parameter< const int& >::type resp(respSEXP);
    Rcpp::traits::input_parameter< const int& >::type ncat(ncatSEXP);
    Rcpp::traits::input_parameter< const int& >::type model(modelSEXP);
    Rcpp::traits::input_parameter< const int& >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type prior_parm(prior_parmSEXP);
    rcpp_result_gen = Rcpp::wrap(theta_EB_single(nx, theta_init, theta_prop, item_parm, resp, ncat, model, prior, prior_parm));
    return rcpp_result_gen;
END_RCPP
}
// theta_FB
arma::mat theta_FB(const int& nx, const arma::rowvec& theta_init, const double& theta_prop, const List& items_list, const arma::mat& item_init, const arma::icolvec& resp, const arma::icolvec& ncat, const arma::icolvec& model, const int& prior, const arma::rowvec& prior_parm);
RcppExport SEXP _TestDesign_theta_FB(SEXP nxSEXP, SEXP theta_initSEXP, SEXP theta_propSEXP, SEXP items_listSEXP, SEXP item_initSEXP, SEXP respSEXP, SEXP ncatSEXP, SEXP modelSEXP, SEXP priorSEXP, SEXP prior_parmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type nx(nxSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type theta_init(theta_initSEXP);
    Rcpp::traits::input_parameter< const double& >::type theta_prop(theta_propSEXP);
    Rcpp::traits::input_parameter< const List& >::type items_list(items_listSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type item_init(item_initSEXP);
    Rcpp::traits::input_parameter< const arma::icolvec& >::type resp(respSEXP);
    Rcpp::traits::input_parameter< const arma::icolvec& >::type ncat(ncatSEXP);
    Rcpp::traits::input_parameter< const arma::icolvec& >::type model(modelSEXP);
    Rcpp::traits::input_parameter< const int& >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type prior_parm(prior_parmSEXP);
    rcpp_result_gen = Rcpp::wrap(theta_FB(nx, theta_init, theta_prop, items_list, item_init, resp, ncat, model, prior, prior_parm));
    return rcpp_result_gen;
END_RCPP
}
// theta_FB_single
arma::mat theta_FB_single(const int& nx, const arma::rowvec& theta_init, const double& theta_prop, const arma::mat& item_mcmc, const arma::rowvec& item_init, const int& resp, const int& ncat, const int& model, const int& prior, const arma::rowvec& prior_parm);
RcppExport SEXP _TestDesign_theta_FB_single(SEXP nxSEXP, SEXP theta_initSEXP, SEXP theta_propSEXP, SEXP item_mcmcSEXP, SEXP item_initSEXP, SEXP respSEXP, SEXP ncatSEXP, SEXP modelSEXP, SEXP priorSEXP, SEXP prior_parmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type nx(nxSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type theta_init(theta_initSEXP);
    Rcpp::traits::input_parameter< const double& >::type theta_prop(theta_propSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type item_mcmc(item_mcmcSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type item_init(item_initSEXP);
    Rcpp::traits::input_parameter< const int& >::type resp(respSEXP);
    Rcpp::traits::input_parameter< const int& >::type ncat(ncatSEXP);
    Rcpp::traits::input_parameter< const int& >::type model(modelSEXP);
    Rcpp::traits::input_parameter< const int& >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type prior_parm(prior_parmSEXP);
    rcpp_result_gen = Rcpp::wrap(theta_FB_single(nx, theta_init, theta_prop, item_mcmc, item_init, resp, ncat, model, prior, prior_parm));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_TestDesign_calc_info", (DL_FUNC) &_TestDesign_calc_info, 4},
    {"_TestDesign_calc_info_matrix", (DL_FUNC) &_TestDesign_calc_info_matrix, 4},
    {"_TestDesign_calc_info_EB", (DL_FUNC) &_TestDesign_calc_info_EB, 4},
    {"_TestDesign_calc_info_FB", (DL_FUNC) &_TestDesign_calc_info_FB, 5},
    {"_TestDesign_calc_MI_FB", (DL_FUNC) &_TestDesign_calc_MI_FB, 4},
    {"_TestDesign_p_1pl", (DL_FUNC) &_TestDesign_p_1pl, 2},
    {"_TestDesign_p_2pl", (DL_FUNC) &_TestDesign_p_2pl, 3},
    {"_TestDesign_p_3pl", (DL_FUNC) &_TestDesign_p_3pl, 4},
    {"_TestDesign_p_pc", (DL_FUNC) &_TestDesign_p_pc, 2},
    {"_TestDesign_p_gpc", (DL_FUNC) &_TestDesign_p_gpc, 3},
    {"_TestDesign_p_gr", (DL_FUNC) &_TestDesign_p_gr, 3},
    {"_TestDesign_array_p_1pl", (DL_FUNC) &_TestDesign_array_p_1pl, 2},
    {"_TestDesign_array_p_2pl", (DL_FUNC) &_TestDesign_array_p_2pl, 3},
    {"_TestDesign_array_p_3pl", (DL_FUNC) &_TestDesign_array_p_3pl, 4},
    {"_TestDesign_array_p_pc", (DL_FUNC) &_TestDesign_array_p_pc, 2},
    {"_TestDesign_array_p_gpc", (DL_FUNC) &_TestDesign_array_p_gpc, 3},
    {"_TestDesign_array_p_gr", (DL_FUNC) &_TestDesign_array_p_gr, 3},
    {"_TestDesign_info_1pl", (DL_FUNC) &_TestDesign_info_1pl, 2},
    {"_TestDesign_info_2pl", (DL_FUNC) &_TestDesign_info_2pl, 3},
    {"_TestDesign_info_3pl", (DL_FUNC) &_TestDesign_info_3pl, 4},
    {"_TestDesign_info_pc", (DL_FUNC) &_TestDesign_info_pc, 2},
    {"_TestDesign_info_gpc", (DL_FUNC) &_TestDesign_info_gpc, 3},
    {"_TestDesign_info_gr", (DL_FUNC) &_TestDesign_info_gr, 3},
    {"_TestDesign_array_info_1pl", (DL_FUNC) &_TestDesign_array_info_1pl, 2},
    {"_TestDesign_array_info_2pl", (DL_FUNC) &_TestDesign_array_info_2pl, 3},
    {"_TestDesign_array_info_3pl", (DL_FUNC) &_TestDesign_array_info_3pl, 4},
    {"_TestDesign_array_info_pc", (DL_FUNC) &_TestDesign_array_info_pc, 2},
    {"_TestDesign_array_info_gpc", (DL_FUNC) &_TestDesign_array_info_gpc, 3},
    {"_TestDesign_array_info_gr", (DL_FUNC) &_TestDesign_array_info_gr, 3},
    {"_TestDesign_find_segment", (DL_FUNC) &_TestDesign_find_segment, 2},
    {"_TestDesign_calc_likelihood", (DL_FUNC) &_TestDesign_calc_likelihood, 5},
    {"_TestDesign_calc_likelihood_function", (DL_FUNC) &_TestDesign_calc_likelihood_function, 5},
    {"_TestDesign_calc_log_likelihood", (DL_FUNC) &_TestDesign_calc_log_likelihood, 7},
    {"_TestDesign_calc_log_likelihood_function", (DL_FUNC) &_TestDesign_calc_log_likelihood_function, 7},
    {"_TestDesign_calc_posterior", (DL_FUNC) &_TestDesign_calc_posterior, 7},
    {"_TestDesign_calc_posterior_function", (DL_FUNC) &_TestDesign_calc_posterior_function, 7},
    {"_TestDesign_calc_posterior_single", (DL_FUNC) &_TestDesign_calc_posterior_single, 7},
    {"_TestDesign_theta_EAP", (DL_FUNC) &_TestDesign_theta_EAP, 7},
    {"_TestDesign_theta_EAP_matrix", (DL_FUNC) &_TestDesign_theta_EAP_matrix, 7},
    {"_TestDesign_theta_EB", (DL_FUNC) &_TestDesign_theta_EB, 9},
    {"_TestDesign_theta_EB_single", (DL_FUNC) &_TestDesign_theta_EB_single, 9},
    {"_TestDesign_theta_FB", (DL_FUNC) &_TestDesign_theta_FB, 10},
    {"_TestDesign_theta_FB_single", (DL_FUNC) &_TestDesign_theta_FB_single, 10},
    {NULL, NULL, 0}
};

RcppExport void R_init_TestDesign(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
