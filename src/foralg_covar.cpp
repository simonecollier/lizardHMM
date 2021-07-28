// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double foralg_covar(int n, int N, arma::mat foo, arma::vec gamma, arma::mat allprobs) {

  arma::cube gammaArray(gamma.begin(), N, N, n, false);

  double lscale=0;
  double sumfoo=0;

  foo = foo%allprobs.row(1);

  for(int j=0; j<N; j++){
    sumfoo += foo(0,j);
  }

  lscale+=log(sumfoo);
  foo = foo/sumfoo;
  sumfoo = 0;

  for (int i=1; i < n; i++){
    /*Rcout << foo<< std::endl;
     Rcout << allprobs.row(i)<< std::endl;*/

    foo = foo*gammaArray.slice(i)%allprobs.row(i);

    /*Rcout<<foo<<std::endl;
     Rcout<<gammaArray.slice(i)<<std::endl;*/

    for(int j=0; j<N; j++){
      sumfoo += foo(0,j);
    }

    /*Rcout<<sumfoo<<std::endl;*/

    lscale+=log(sumfoo);
    foo = foo/sumfoo;
    sumfoo = 0;
  }

  return(lscale);

}
