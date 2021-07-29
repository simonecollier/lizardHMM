/* Forward algorithm for an N-state
 * hidden Markov model with no stationary distribution
 */


#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//' Execute the forward algorithm
//'
//' @param n The length of the timeseries.
//' @param N The number of states in the HMM.
//' @param foo The initial state distributions (delta).
//' @param gamma The transition probability matrix.
//' @param allprobs As outputted by the funtions norm_allporbs
//' @export
// [[Rcpp::export]]
double foralg(int n, int N, arma::mat foo, arma::mat gamma, arma::mat allprobs) {

  double lscale=0;
  double sumfoo=0;

  foo = foo%allprobs.row(0);

  for(int j=0; j<N; j++){
    sumfoo += foo(0,j);
  }

  lscale+=log(sumfoo);
  foo = foo/sumfoo;
  sumfoo = 0;

  for (int i=1; i < n; i++){
    foo = foo*gamma%allprobs.row(i);

    for(int j=0; j<N; j++){
      sumfoo += foo(0,j);
    }

    lscale+=log(sumfoo);
    foo = foo/sumfoo;
    sumfoo = 0;
  }

  return(lscale);

}
