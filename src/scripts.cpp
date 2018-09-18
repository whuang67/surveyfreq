#include "scripts.h"
#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// O(n)
double jk_var_w_mean(NumericVector jk_rep, double mean) {
  int l=jk_rep.size();
  if(l<2)
    stop("There must be at least 2 replicate weights.");
  double out = pow(jk_rep[0]-mean, 2);
  for(int i=1; i<jk_rep.size(); i++)
    out += pow(jk_rep[i]-mean, 2);
  return out*(l-1)/l;
}
double jk_se_w_mean(NumericVector jk_rep, double mean) {
  return pow(jk_var_w_mean(jk_rep, mean), 0.5);
}


double brr_var_w_mean(NumericVector brr_rep, double mean, double BRRfay) {
  int l=brr_rep.size();
  if(l<2)
    stop("There must be at least 2 replicate weights.");
  double out = pow(brr_rep[0]-mean, 2);
  for(int i=1; i<brr_rep.size(); i++)
    out += pow(brr_rep[i]-mean, 2);
  return out/(l*pow(1-BRRfay, 2));
}
double brr_se_w_mean(NumericVector brr_rep, double mean, double BRRfay) {
  return pow(brr_var_w_mean(brr_rep, mean, BRRfay), 0.5);
}

// [[Rcpp::export]]
int fibonacci(int n){
  if(n<0)
    stop("Not available, input cannot be less than 0.\n");
  else if(n==0)
    return 0;
  else if(n==1||n==2)
    return 1;
  else
    return fibonacci(n-1)+fibonacci(n-2);
}

// ------ Sort Algorithms
// Bubble Sort
template <int A>
void bubble_sort(Vector<A> x) {
  typename traits::storage_type<A>::type tmp;
  for(int i=x.size()-1; i>=0; i--){
    for(int j=0; j<i; j++){
      if(x[i]<x[j]){tmp=x[j]; x[j]=x[i]; x[i]=tmp;}
    }
  }
}
// Selection Sort
template <int A>
void selection_sort(Vector<A> x) {
  int idxOfMax; typename traits::storage_type<A>::type tmp;
  for(int i=x.size()-1; i>=0; i--){
    idxOfMax=i;
    for(int j=0; j<i; j++){
      if(x[idxOfMax]<x[j]){idxOfMax=j;}
    }
    if(idxOfMax!=i){tmp=x[idxOfMax]; x[idxOfMax]=x[i]; x[i]=tmp;}
  }
}
// Insertion Sort
template <int A>
void insertion_sort(Vector<A> x, int start=0, int gap=1) {
  typename traits::storage_type<A>::type currentValue; int j;
  for(int i=start+gap; i<x.size(); i=i+gap){
    currentValue=x[i]; j=i;
    while(j>=gap && x[j-gap]>currentValue){x[j]=x[j-gap]; j=j-gap;}
    x[j]=currentValue;
  }
}
// Shell Sort
template <int A>
void shell_sort(Vector<A> x) {
  int gap=x.size()/2;
  while(gap>0){
    for(int start=0; start<gap; start++){
      insertion_sort(x, start, gap);
    }
  }
}
// Merge Sort
template <int A>
void merge_sort(Vector<A> x) {
  if(x.size()>1){
    int midpoint=x.size()/2;
    Vector<A> lefthalf=x[Range(0, midpoint-1)], righthalf=x[Range(midpoint, x.size()-1)];
    merge_sort(lefthalf); merge_sort(righthalf);
    int i=0,j=0,k=0;
    while(i<lefthalf.size() && j<righthalf.size()){
      if(lefthalf[i]<righthalf[j]){
        x[k]=lefthalf[i]; i++;
      } else{
        x[k]=righthalf[j]; j++;
      }
      k++;
    }
    while(i<lefthalf.size()){
      x[k]=lefthalf[i]; i++; k++;
    }
    while(j<righthalf.size()){
      x[k]=righthalf[j]; j++; k++;
    }
  }
}
// Quick Sort
template <int A>
int partition(Vector<A> x, int first, int last){
  typename traits::storage_type<A>::type pivotvalue=x[first], tmp;
  int leftmark=first+1, rightmark=last; bool done=false;
  while(!done){
    while(leftmark<=rightmark &&x[leftmark]<=pivotvalue){leftmark+=1;}
    while(x[rightmark]>=pivotvalue && rightmark>=leftmark){rightmark-=1;}
    if(leftmark>rightmark){
      done=true;
    } else{
      tmp=x[leftmark]; x[leftmark]=x[rightmark]; x[rightmark]=tmp;
    }
  }
  tmp=x[first]; x[first]=x[rightmark]; x[rightmark]=tmp;
  return rightmark;
}
template <int A>
void quick_sort_help(Vector<A> x, int first, int last){
  int splitpoint;
  if(first<last){
    splitpoint = partition(x, first, last);
    quick_sort_help(x, first, splitpoint);
    quick_sort_help(x, splitpoint+1, last);
  }
}
template <int A>
void quick_sort(Vector<A> x){
  quick_sort_help(x,0,x.size()-1);
}

template <int A>
Vector<A> sort_help(Vector<A> arr, String method, bool overwrite){
  Vector<A> x;
  if(overwrite) x=clone(arr);
  else x=arr;

  if(method=="bubble") bubble_sort(x);
  else if(method=="selection") selection_sort(x);
  else if(method=="insertion") insertion_sort(x);
  else if(method=="shell") shell_sort(x);
  else if(method=="merge") merge_sort(x);
  else if(method=="quick") quick_sort(x);
  else{
    Rprintf("Unknown method, MERGE sort is used.\n");
    merge_sort(x);
  }
  return x;
}

// [[Rcpp::export]]
SEXP sort_cpp(SEXP arr, String method="merge", bool overwrite=true){
  switch(TYPEOF(arr)){
  case INTSXP: return sort_help(as<IntegerVector>(arr), method, overwrite);
  case REALSXP: return sort_help(as<NumericVector>(arr), method, overwrite);
  case STRSXP: return sort_help(as<StringVector>(arr), method, overwrite);
  default: stop("Only support INT, NUM and CHAR datatype.\n");
  }
}
