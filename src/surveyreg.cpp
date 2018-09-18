#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec SGD_lm(arma::mat X, arma::vec y, bool intercept=true, int iter=1000, double learning_rate=0.001){
  if(X.n_rows!=y.n_rows)
    stop("Wrong dimensions!");

  arma::vec beta = arma::randn(X.n_cols);
  if(intercept){
    X.insert_cols(0, arma::ones(y.n_rows));
    beta.insert_rows(0, arma::ones(1));
  }

  arma::vec gradient=2*X.t()*(X*beta - y);
  for(int i=0; i<iter; i++){
    beta -= gradient * learning_rate;
    gradient = 2*X.t()*(X*beta - y);
  }
  return beta;
}

arma::vec OLS_lm(arma::mat X, arma::vec y, bool intercept=true){
  if(X.n_rows!=y.n_rows)
    stop("Wrong dimensions!");
  if(intercept)
    X.insert_cols(0, arma::ones(y.n_rows));
  return arma::solve(X, y);
}
