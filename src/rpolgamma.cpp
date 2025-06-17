#include <stdio.h>
#include <cmath>
#include <R.h>
#include <Rmath.h>
#include <Rcpp.h>
#include <RcppEigen.h>

using std::pow;
using std::fabs;
using std::sqrt;
using std::log;
using std::exp;
using namespace Eigen;
using namespace Rcpp;

/* 
 * 
 * see also: https://cran.r-project.org/web/packages/pgdraw/index.html
 * see also: http://dirk.eddelbuettel.com/papers/rcpp_uzuerich_2015_part5_packaging.pdf
 *
 */
 
double anx(int n, double x, double t){
  // evaluates the piecewise coefficients in Polson et al. (2013, p.13)
  double z;
  if(x <= 0) {
    return(0);
  } else if(x <= t) {
    z = M_PI*(n+0.5)*std::pow(2/(M_PI*x),3/2)*std::exp(-2*(n+0.5)*(n+0.5)/x);
  } else {
    z = M_PI*(n+0.5)*std::exp(-x*M_PI*M_PI*(n+0.5)*(n+0.5)/2);
  }
  return(z);
}

double rinvtnorm(const double m, const double t){
  // samples from truncated normal distribution, mean mu and sd 1, truncated at t
  // shouldn't be visible from R, so I'm not checking for mu, t > 0
  // these first lines will be removed depending on how it's called in rpolgamma
  double ret = 0;
  Rcpp::NumericVector e(1);
  Rcpp::NumericVector e2(1);
  Rcpp::NumericVector X(1);
  Rcpp::NumericVector Y(1);
  Rcpp::NumericVector nu(1);
  Rcpp::NumericVector alpha(1);
  Rcpp::LogicalVector test(1);
  Rcpp::LogicalVector test2(1);
  if(m > t){ // m must be length one, from code
    do{
      do{
        e = rexp(1);
        e2 = rexp(1);
        test = (e*e <= 2*m*e2);
      } while (test);
      X = (1/m)*(1 + e/m)*(1 + e/m);
      alpha = Rcpp::exp(-0.5*X/(m*m));
      test2 = (Rcpp::runif(1) <= alpha);
    } while (test2);
    ret = as<double>(X);
  } else { // mu(i) <= t
    // Devroye's book is online
    // http://luc.devroye.org/chapter_four.pdf, p.149
    // This is the only part that runs in the code, actually... 
    // rinvtnorm is only called if the first if() fails
    nu = Rcpp::rnorm(1);
    Y = nu*nu;
    X = m + 0.5*m*m*Y - 0.5*m*sqrt(4*m*Y + m*m*Y*Y);
    alpha = m/(m+X);
    test2 = (Rcpp::runif(1) > alpha);
    if(test2(0)){
      X = m*m/X;
    } 
    ret = as<double>(X);
  }
  // }
  return(ret);
}

Rcpp::NumericVector pigauss(const double t, const Rcpp::NumericVector mu){
  // inverse gaussian cdf, but with lambda = 1
  // Rpigauss <- function(t,mu) pnorm((t/mu-1)/sqrt(t)) + exp(2/mu)*pnorm(-(t/mu+1)/sqrt(t))
  int L = mu.size();
  Rcpp::NumericVector retPiGauss(L);
  for(int i = 0; i < L; i++){
    retPiGauss(i) = R::pnorm((t/mu(i) - 1)/std::sqrt(t), 0.0, 1.0, 1, 0) + std::exp(2/mu(i))*R::pnorm(-(t/mu(i) + 1)/std::sqrt(t), 0.0, 1.0, 1, 0);
  }
  return(retPiGauss);
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Rcpp::NumericVector rpolgamma(const int n, const double b, const Rcpp::NumericVector mu) {
  // Part 1: change mu to z, so that it matches dimensions with n
  Rcpp::NumericVector m(n);
  if(mu.size() == 1){ // constant mu parameter
    for(int i = 0; i < n; i++){
      m(i) = mu(0);
    }
  } else if (mu.size() == n){ // variable mu parameter
    m = mu;
  } else { // raise error
    stop("mu parameter must have length 1 or n.");
  }
  // From here onward, must call m, not mu
  // Part 2: use additivity of Polya-Gamma for parameter b, 
  //         and consider fractional b values as well.
  if(b <= 0) stop("Invalid value for b in rpolgamma.");
  if((b > 0) & (b < 1)) {
    // Rcpp::NumericVector temp(n);
    // // ... TODO
    // return(temp);
    stop("Non-integer b in rpolgamma not yet implemented.");
  }
  if((b > 1) & (b - floor(b) != 0)) {
    // Rcpp::NumericVector temp(n);
    // temp += rpolgamma(n, b - floor(b), m); // start with fractional b
    // for(int i = 0; i < floor(b); i++){ // use additivity of Polya-Gamma
    //   temp += rpolgamma(n, 1, m);
    // }
    // return(temp);
    stop("Non-integer b in rpolgamma not yet implemented.");
  }
  if((b > 1) & (b - floor(b) == 0)) {
    Rcpp::NumericVector temp(n);
    for(int i = 0; i < floor(b); i++){ // use additivity of Polya-Gamma
      temp += rpolgamma(n, 1, m);
    }
    return(temp);
  }
  // else b == 1
  
  // "Devroye shows that the best choice of t is near 0.64" (Polson et al., 2013, p.13)
  Rcpp::NumericVector z = abs(m)/2;
  double t = 0.64;
  Rcpp::NumericVector K = M_PI*M_PI/8 + z*z/2;
  Rcpp::NumericVector p = M_PI*Rcpp::exp(-K*t)/(2*K);
  Rcpp::NumericVector newMu = 1/z;
  Rcpp::NumericVector q = 2*Rcpp::exp(-abs(z))*pigauss(t, newMu);
  // Rcpp::Rcout << z << std::endl;
  // Rcpp::Rcout << q << std::endl;
  Rcpp::NumericVector ret(n);
  Rcpp::NumericVector S(n);
  Rcpp::NumericVector Y(1);
  int k = 0;
  int i = 0;
  
  while(i < n){
    
    Rcpp::NumericVector U = Rcpp::runif(1);
    Rcpp::NumericVector V = Rcpp::runif(1);
    // Rcpp::Rcout << U(0) << '\t';
    // Rcpp::Rcout << p(i)/(p(i) + q(i)) << '\n';
    
    if(U(0) < p(i)/(p(i) + q(i))){
      Rcpp::NumericVector E = rexp(1);
      ret(i) = t + E(0)/K(i); 
    } else {
      if(newMu(i) > t) {
        int condition = 1;
        do{
          Rcpp::NumericVector chi = rchisq(1, 1.0);
          while(chi(0) < t){
            chi = Rcpp::rchisq(1, 1.0);
          }
          ret(i) = 1/chi(0);
          double temp = std::exp(-z(i)*z(i)*ret(i)/2);
          Rcpp::NumericVector compare = Rcpp::runif(1);
          if(compare(0) < temp){
            condition = 0;
          }
        } while (condition); // repeat until a == repeat while !a
      } else { // Remember if m is NumericVector, m(i) is double!!
        int condition = 1;
        // Rcpp::Rcout << "caso A" << std::endl;
        do{
          // Rcpp::Rcout << "m(i) = " << m(i) << std::endl;
          // Rcpp::Rcout << "t = " << t << std::endl;
          ret(i) = rinvtnorm(newMu(i), t); // as<double>(rinvtnorm(1, m(i), t));
          // Rcpp::Rcout << ret(i) << std::endl;
          if(ret(i) < t){
            condition = 0;
          }
        } while (condition); // repeat until a == repeat while !a
      }
    }
    S(i) = anx(0, ret(i), t);
    Y = V*S(i);
    k = 0;
    for(;;){ // sample until not rejecting
      k++;
      if(k % 2){
        S(i) = S(i) - anx(k, ret(i), t);
        if(Y(0) < S(i)) {
          ret(i) = ret(i)/4;
          i++;
          break;
        }
      } else {
        S(i) = S(i) + anx(k, ret(i), t);
        if(Y(0) > S(i)) {
          break; // break for(;;), won't increment i so regenerates Y
        }
      }
    }
    
  }
  
  // Rcpp::NumericVector z = wrap(mu);
  // Rcpp::NumericVector z = rnorm(n);
  return(ret);
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
# set.seed(1)
# rpolgamma(100, 1, .2)
# rpolgamma(100, 1, 2.2)
# rpolgamma(100, 5, .2)
# rpolgamma(100, 2, 2.2)
# rpolgamma(20, 1, seq(0.01,10, length = 20)) # fine if mu is constant
*/
