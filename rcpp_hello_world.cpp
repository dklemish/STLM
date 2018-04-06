//#include <Rcpp.h>
#include <R.h>
#include <Rmath.h>
#include <RcppGSL.h>
#include <RcppArmadillo.h>
#include <algorithm>

using namespace Rcpp;
using namespace arma;

const double n_cdf   = 5000;

// Declare dependencies so Rcpp knows to link libraries
// [[Rcpp::depends(RcppGSL)]]
// [[Rcpp::depends(RcppArmadillo)]]

/** Function headers **/
double E1(double x){
  return(R::pgamma(x, 1e-9, 1.0, FALSE, FALSE) / 1e-9);
}

double E1inv(double x){
  return(R::qgamma(1e-9*x, 1e-9, 1.0, FALSE, FALSE));
}

// [[Rcpp::export]]
SEXP E1_wrap(SEXP x_){
  double x;
  
  x = E1(as<double>(x_));
  
  return wrap(x);
}

int find(double x, NumericVector breaks) {
  NumericVector::iterator it;
  NumericVector::iterator pos;
  int out;
  
  pos = std::upper_bound(breaks.begin(), breaks.end(), x);
  out = std::distance(breaks.begin(), pos);
  
  return out;
}

NumericVector drawGamma(double epsilon, double shape, 
                        double rate, double norm_constant){
  NumericVector cdf(n_cdf);
  NumericVector X(n_cdf);
  
  double stepSize = (1-2*epsilon) / n_cdf;
  double currentStep;
  double u;
  
  int i, J;
  
  for(i = 0; i < n_cdf; i++){
    currentStep = (i+1)*stepSize;
    cdf[i] = currentStep;
    X[i] = E1inv(E1(rate*epsilon)*(1-currentStep));
  }
  
  GetRNGstate();
  
  J = R::rpois(norm_constant);
  
  NumericVector output(J);
  
  for(i=0; i < J; i++){
    u = R::runif(epsilon, 1.0-epsilon);
    output[i] = X[find(u, cdf)];
  }
  
  return(output);
}

// [[Rcpp::export]]
SEXP drawRF_Gamma1D(SEXP X_, 
                    SEXP shape_, 
                    SEXP rate_, 
                    SEXP epsilon_, 
                    SEXP radius_)
{
  NumericVector X = as<NumericVector>(X_);
  NumericVector numMassPts = NumericVector(X.nrow());
  NumericVector levyDraws;
  IntegerVector locToStore;
  double shape    = as<double>(shape_);
  double rate     = as<double>(rate_);
  double epsilon  = as<double>(epsilon_);
  double radius   = as<double>(radius_);
  double norm_constant;
  NumericVector draw_radius;
  NumericVector draw_angle;
  int totalMassPts = 0;
  int num_draws;
  int i;
  
  norm_constant = shape*E1(rate*epsilon);
  
  NumericMatrix massPts = mat(10*round(norm_constant*X.nrow()), 2, fill::zeros);
  
  //Draw mass points for first location
  levyDraws = drawGamma(epsilon, shape, rate, norm_constant);
  
  num_draws = levyDraws.size();
  draw_radius = runif(num_draws, 0.0, radius);
  for(i = 0; i < num_draws; i++){
    draw_radius[i] = sqrt(draw_radius[i]);
  }
  draw_angle  = runif(num_draws, 0.0, 2*PI);

  locToStore = seq(totalMassPts + 1, totalMassPts + num_draws);
  massPts[locToStore, 1] = levyDraws;
  massPts[locToStore, 2] = draw.radius*cos(draw.angle) + X[1];

  numAdded[1] <- num.draws
  total.mass.pts <- total.mass.pts + num.draws
    
    for(i in 2:nrow(X)){
# Draw mass points for other locations, but discard points that fall 
# in the regions of previous locations
      levy.draws <- drawGammaLM(1e-5, 10, 1, norm.constant)
      num.draws <- length(levy.draws)
      draw.radius <- sqrt(runif(num.draws, 0, R))
      draw.angle <- runif(num.draws, 0, 2*pi)
      
      x.loc <- X[i,1] + draw.radius*cos(draw.angle)
      y.loc <- X[i,2] + draw.radius*sin(draw.angle)
      
      valid.new.pt <- rep(TRUE, num.draws)
      
      for(j in 1:(i-1)){
        valid.new.pt <- valid.new.pt & !(((x.loc-X[j,1])^2 + (y.loc-X[j,2])^2) < R^2)
      }
      
      num.new <- sum(valid.new.pt)
        num.added[i] <- num.new
        
        if(num.new > 0){
          mass.pts[(total.mass.pts+1):(total.mass.pts+num.new), 1] <- levy.draws[valid.new.pt]
          mass.pts[(total.mass.pts+1):(total.mass.pts+num.new), 2] <- x.loc[valid.new.pt]
          mass.pts[(total.mass.pts+1):(total.mass.pts+num.new), 3] <- y.loc[valid.new.pt]
          total.mass.pts <- total.mass.pts + num.new
        }
    }
    Y <- rep(0, nrow(X))
      for(i in 1:nrow(X)){
        Y[i] <- sum(mass.pts[(mass.pts[,2]-X[i,1])^2 + (mass.pts[,3]-X[i,2])^2 < R^2,1])
      }
      
      return(cbind(X, Y, num.added))
}