#ifndef S_H
#define S_H
#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

double jk_var_w_mean(NumericVector jk_rep, double mean);
double jk_se_w_mean(NumericVector jk_rep, double mean);
double brr_var_w_mean(NumericVector brr_rep, double mean, double BRRfay);
double brr_se_w_mean(NumericVector brr_rep, double mean, double BRRfay);

#endif
