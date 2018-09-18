#include "scripts.h"
#include <Rcpp.h>
#include <map>

using namespace Rcpp;
using namespace std;

// freq tab
template <int A>
NumericVector freq_tab_help(Vector<A> arr, NumericVector wt, bool includeNA, double NA_wt){
  int l=arr.size();
  if(l!=wt.size())
    stop("Lengths of array and weight are different!\n");

  typedef typename traits::storage_type<A>::type map_type;
  typedef std::map<map_type, double, internal::NAComparator<map_type> > SORTED_MAP;
  SORTED_MAP m; map_type s;
  double tot=0;
  for(int i=0; i<l; i++){
    if(Vector<A>::is_na(arr[i]) && !includeNA)
      continue;
    s=arr[i];
    if(!NumericVector::is_na(wt[i])){
      m[s]+=wt[i]; tot+=wt[i];
    } else {
      m[s]+=NA_wt; tot+=NA_wt;
    }
  }

  // Convert to DataFrame output
  NumericVector out=wrap(m);
  StringVector var=out.attr("names");
  out.push_back(tot);
  var.push_back("Total");

  out.attr("names")=var;
  return out;
}

template <int A>
NumericMatrix freq_tab(Vector<A> arr, Nullable<NumericVector> wt, bool includeNA, double NA_wt) {
  NumericMatrix out;
  if(wt.isNull()) {
    NumericVector tmp = freq_tab_help(arr, NumericVector(arr.size(), 1.0), includeNA, NA_wt);
    out=NumericMatrix(tmp.size(), 1);
    out(_, 0) = tmp;
    out.attr("dimnames") = List::create(tmp.attr("names"), StringVector::create("Count"));
  } else{
    NumericVector tmp = freq_tab_help(arr, as<NumericVector>(wt), includeNA, NA_wt);
    int l=tmp.size();
    out=NumericMatrix(l, 2);
    out(_, 0) = tmp;
    for(int i=0; i<l; i++)
      out(_, 1)=out(_, 0)/out(l-1, 0);
    out.attr("dimnames") = List::create(tmp.attr("names"), StringVector::create("Freq", "Pct"));
  }
  return out;
}



template <int A>
NumericMatrix rep_freq_tab(Vector<A> arr,
                           Nullable<NumericVector> based_wt,
                           Nullable<NumericMatrix> rep_wts,
                           bool includeNA,
                           double NA_wt,
                           String varmethod,
                           double BRRfay) {
  int l=arr.size();
  NumericVector bwt(l);

  if(rep_wts.isNotNull()) {
    NumericMatrix rpwt=as<NumericMatrix>(rep_wts);
    if(based_wt.isNotNull()){
      Rprintf("based_wt: Provided, rep_wts: Provided\n");
      bwt=as<NumericVector>(based_wt);
    } else{
      Rprintf("based_wt: Null, rep_wts: Provided\nAverage of rep_wts is used as based_wt.\n");
      NumericMatrix rpwt=as<NumericMatrix>(rep_wts);
      for(int i=0; i<l; i++)
        bwt[i]=mean(rpwt(i, _));
    }

    if(rpwt.nrow()!=bwt.size())
      stop("Unequal dimensions of based and replicated wgts!!!\n");

    int n_col = rpwt.ncol();
    StringVector col_name(n_col);
    NumericMatrix based_out=freq_tab(arr, bwt, includeNA, NA_wt);
    int n_row = based_out.nrow();
    NumericMatrix rep_freq(n_row, n_col), rep_pct(n_row, n_col), tmp(n_row, 2), out(n_row, 5);
    out(_, 0)=freq_tab(arr, R_NilValue, includeNA, NA_wt);
    out(_, 1)=based_out(_, 0);
    out(_, 3)=based_out(_, 1);

    NumericVector tmp_rpwt(l);
    for(int i=0; i<n_col; i++) {
      tmp_rpwt=rpwt(_, i);
      tmp=freq_tab(arr, tmp_rpwt, includeNA, NA_wt);
      rep_freq(_, i)=tmp(_, 0); rep_pct(_, i)=tmp(_, 1);
    }

    if(varmethod=="jackknife") {
      for(int i=0; i<n_row; i++) {
        out(i, 2)=jk_se_w_mean(rep_freq(i, _), out(i, 1));
        if(i!=n_row-1)
          out(i ,4)=jk_se_w_mean(rep_pct(i, _), out(i, 3));
        else
          out(i, 4)=NA_REAL;
      }

    } else if(varmethod=="BRR"){
      if(BRRfay<0 || BRRfay>1)
        stop("BRRfay must be between 0 and 1.\n");
      for(int i=0; i<n_row; i++) {
        out(i, 2)=brr_se_w_mean(rep_freq(i, _), out(i, 1), BRRfay);
        if(i!=n_row-1)
          out(i ,4)=brr_se_w_mean(rep_pct(i, _), out(i, 3), BRRfay);
        else
          out(i, 4)=NA_REAL;
      }

    } else
      stop("Only jackknife and BRR are supported currently");

    out.attr("dimnames")=List::create(rownames(based_out), StringVector::create("Count", "Freq", "Freq_SE", "FreqPct", "FreqPct_SE"));
    return out;

  } else{
    Rprintf("No \"rep_wts\" found, parameter \"varmethod\" is ignored.\n");
    if(based_wt.isNotNull())
      Rprintf("based_wt: Provided, rep_wts: Null\nOne-way freq table only.\n");
    else
      Rprintf("based_wt: Null, rep_wts: Null\nStraight 1's are used as based_wt.\n");
    return freq_tab(arr, based_wt, includeNA, NA_wt);
  }
}

//' Generate the weighted frequency table
//'
//' Generate the weighted frequency table with jackknife or BRR similar with PROC SURVEYFREQ procedure.
//' It takes roughly 10+ seconds to process 10
//' @param arr A \code{vector} (integer, double, character, logical or complex) of dimension M
//' @param based_wt A \code{vector} of dimension of M
//' @param rep_wts A \code{matrix} of M x N
//' @param includeNA Logical. Should missing values (including NaN) be removed?
//' @param NA_wt double Assign a constance value if weight is missing
//' @param varmethod Character, "jackknife" or "BRR"
//' @param BRRfay 0 - 1, only available when varmethod=="BRR"
//' @author Wenke Huang
//' @export
//' @examples
//' \dontrun{
//' # Generate a vector, based weight and replicated weights
//' arr=sample(1:10, 1000, replace=TRUE)
//' based_wt=runif(1000)
//' rep_wts=matrix(runif(1000*80), ncol=80)
//'
//' # Generate table
//' surveyfreq2(arr, based_wt, rep_wts)
//' }
// [[Rcpp::export]]
NumericMatrix surveyfreq(SEXP arr,
                         Nullable<NumericVector> based_wt=R_NilValue,
                         Nullable<NumericMatrix> rep_wts=R_NilValue,
                         bool includeNA=false,
                         double NA_wt=1,
                         String varmethod="jackknife",
                         double BRRfay=0){
  switch(TYPEOF(arr)){
  case REALSXP:
    return rep_freq_tab(as<NumericVector>(arr), based_wt, rep_wts, includeNA, NA_wt, varmethod, BRRfay);
  case STRSXP:
    return rep_freq_tab(as<StringVector>(arr), based_wt, rep_wts, includeNA, NA_wt, varmethod, BRRfay);
  case INTSXP:
    return rep_freq_tab(as<IntegerVector>(arr), based_wt, rep_wts, includeNA, NA_wt, varmethod, BRRfay);
  case LGLSXP:
    return rep_freq_tab(as<LogicalVector>(arr), based_wt, rep_wts, includeNA, NA_wt, varmethod, BRRfay);
  case CPLXSXP:
    return rep_freq_tab(as<ComplexVector>(arr), based_wt, rep_wts, includeNA, NA_wt, varmethod, BRRfay);
  default:
    stop("You are not supposed to see this warning!\n");
  }
}
