// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// CoxFitCpp
Rcpp::List CoxFitCpp(const arma::rowvec& lp, const arma::irowvec& status, const arma::umat& riskset);
RcppExport SEXP penalizedcpp_CoxFitCpp(SEXP lpSEXP, SEXP statusSEXP, SEXP risksetSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::rowvec& >::type lp(lpSEXP);
    Rcpp::traits::input_parameter< const arma::irowvec& >::type status(statusSEXP);
    Rcpp::traits::input_parameter< const arma::umat& >::type riskset(risksetSEXP);
    __result = Rcpp::wrap(CoxFitCpp(lp, status, riskset));
    return __result;
END_RCPP
}
// Lasso
Rcpp::List Lasso(arma::vec beta, const arma::vec& lambda, const arma::vec& lambda2, const arma::uvec& positive, const arma::mat& X, const Rcpp::Function& fit, const bool trace, const double epsilon, const double maxiter);
RcppExport SEXP penalizedcpp_Lasso(SEXP betaSEXP, SEXP lambdaSEXP, SEXP lambda2SEXP, SEXP positiveSEXP, SEXP XSEXP, SEXP fitSEXP, SEXP traceSEXP, SEXP epsilonSEXP, SEXP maxiterSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type lambda2(lambda2SEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type positive(positiveSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Function& >::type fit(fitSEXP);
    Rcpp::traits::input_parameter< const bool >::type trace(traceSEXP);
    Rcpp::traits::input_parameter< const double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< const double >::type maxiter(maxiterSEXP);
    __result = Rcpp::wrap(Lasso(beta, lambda, lambda2, positive, X, fit, trace, epsilon, maxiter));
    return __result;
END_RCPP
}
// StepLasso
Rcpp::List StepLasso(arma::vec beta, const arma::vec& lambda, const arma::vec& lambda2, const arma::uvec& positive, const arma::mat& X, const Rcpp::Function& fit, const bool trace, const double epsilon, const double maxiter);
RcppExport SEXP penalizedcpp_StepLasso(SEXP betaSEXP, SEXP lambdaSEXP, SEXP lambda2SEXP, SEXP positiveSEXP, SEXP XSEXP, SEXP fitSEXP, SEXP traceSEXP, SEXP epsilonSEXP, SEXP maxiterSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type lambda2(lambda2SEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type positive(positiveSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Function& >::type fit(fitSEXP);
    Rcpp::traits::input_parameter< const bool >::type trace(traceSEXP);
    Rcpp::traits::input_parameter< const double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< const double >::type maxiter(maxiterSEXP);
    __result = Rcpp::wrap(StepLasso(beta, lambda, lambda2, positive, X, fit, trace, epsilon, maxiter));
    return __result;
END_RCPP
}
// Ridge
Rcpp::List Ridge(arma::vec beta, arma::vec eta, const arma::mat& Lambda, const arma::mat& X, const Rcpp::Function& fit, const bool trace, const double epsilon, const double maxiter, const Rcpp::List& fitInput);
RcppExport SEXP penalizedcpp_Ridge(SEXP betaSEXP, SEXP etaSEXP, SEXP LambdaSEXP, SEXP XSEXP, SEXP fitSEXP, SEXP traceSEXP, SEXP epsilonSEXP, SEXP maxiterSEXP, SEXP fitInputSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Lambda(LambdaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Function& >::type fit(fitSEXP);
    Rcpp::traits::input_parameter< const bool >::type trace(traceSEXP);
    Rcpp::traits::input_parameter< const double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< const double >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type fitInput(fitInputSEXP);
    __result = Rcpp::wrap(Ridge(beta, eta, Lambda, X, fit, trace, epsilon, maxiter, fitInput));
    return __result;
END_RCPP
}
