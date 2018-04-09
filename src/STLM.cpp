// Declare dependencies so Rcpp knows to link libraries
// [[Rcpp::depends(RcppGSL)]]

#include <R.h>
#include <Rmath.h>
#include <RcppGSL.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_sf_exp.h>
#include <algorithm>
#include <cmath>
#include <limits>

using namespace Rcpp;

/** Global variables **/
double shape, rate; // Gamma shape & rate parameters
double epsilon;     // Epsilon Levy threshold
double radius;      // Spatial radius
double tau;
double u;
int type;         // type indicator for covariance shape
int dimSpatial;   // number of dimensions needed for cov shape
int nLoc;         // number of locations
int dimX;         // number of columns of X
int numMassPts;   // total number of mass points
int numDraws;
int totalDraws;

/** Function headers **/
double E1(double);
double E1inv(double, double, bool);
std::vector<double> drawGamma(int, double, double, double);

SEXP E1(SEXP);
SEXP E1inv(SEXP, SEXP, SEXP);
SEXP drawGamma(SEXP, SEXP, SEXP, SEXP);

// Exponential integral functions
double E1(double x){return(gsl_sf_expint_E1(x));}
double E1inv(double y, double epsilon, bool Newton){
  double x_old = R::qgamma(1.0e-9*y, 1e-9, 1.0, 0, 0);
  if(Newton == FALSE){
    return(x_old);
  }else{
    double x_new = 0;
    double diff = 1.0e9;
    
    while(diff > epsilon){
      x_new = x_old * (1 + gsl_sf_exp(x_old)*(gsl_sf_expint_E1(x_old) - y));
      diff = std::abs(x_new - x_old);
      x_old = x_new;
    }
    
    return(x_new);
  }
}

// [[Rcpp::export]]
SEXP E1(SEXP x){
  return(wrap(E1(as<double>(x))));
}

// [[Rcpp::export]]
SEXP E1inv(SEXP y, SEXP epsilon, SEXP Newton){
  return(wrap(E1inv(as<double>(y), as<double>(epsilon), as<bool>(Newton))));
}

// [[Rcpp::export]]
SEXP drawGamma(SEXP init_elem, SEXP shape, SEXP rate, SEXP epsilon){
  double shape_ = as<double>(shape);
  double rate_ = as<double>(rate);
  double epsilon_ = as<double>(epsilon);
  
  return(wrap(drawGamma(as<int>(init_elem), 
                        shape_, 
                        rate_,
                        epsilon_)));
}

std::vector<double> drawGamma(int init_elem, double shape, double rate, double epsilon){
  // Draw Gamma random variable using the Inverse Levy Measure 
  // approach of Wolpert & Ickstadt
  
  // Returns a vector of mass values such that the sum of all mass values 
  // is a draw from a Gamma(shape, rate) distribution (up to an error of epsilon)
  const double E1inv_epsilon = 1e-15;
  
  std::vector<double> result(init_elem);
  
  numDraws = 0;
  tau = 0;
  u = std::numeric_limits<double>::infinity();
  
  tau = R::rexp(1.0);
  u = E1inv(tau/shape, E1inv_epsilon, FALSE)/rate;
  result[numDraws++] = u;
  
  while(u > epsilon){
    tau += R::rexp(1.0);
    u = E1inv(tau/shape, E1inv_epsilon, FALSE)/rate;
    result[numDraws++] = u;
    
    if(numDraws == init_elem){
      init_elem = 2*init_elem;
      result.resize(init_elem);
    }
  }
  return(result);
}

// [[Rcpp::export]]
SEXP drawGammaRF_1D(SEXP X_,
                    SEXP shape_,
                    SEXP rate_,
                    SEXP epsilon_,
                    SEXP radius_,
                    SEXP type_,
                    SEXP distMatrix_){
  // Variable definitions
  std::vector<double> levyDraws;
  std::vector<std::vector<double> > massPts;
  int i, j, k;
  double tempLoc;
  bool validLoc;
  NumericVector Y;
  
  // Input data
  // Convert SEXP objects to Rcpp objects
  NumericMatrix X = as<NumericMatrix>(X_);
  nLoc = X.nrow();
  dimX = 1;
  NumericMatrix distMatrix = as<NumericMatrix>(distMatrix_);
  
  shape   = as<double>(shape_);
  rate    = as<double>(rate_);
  epsilon = as<double>(epsilon_);
  radius  = as<double>(radius_);
  type    = as<int>(type_);
  
  numMassPts = nLoc*10; // initial # of mass points needed for all locations; grow as necessary
  totalDraws = 0;
  
  Y = NumericVector(nLoc);
  
  dimSpatial = dimX + 2; // Number of spatial dimensions + 2 (for Levy mass value & "height")
  
  massPts = std::vector<std::vector<double> > (numMassPts, std::vector<double>(dimSpatial));
  
  GetRNGstate();
  
  levyDraws = drawGamma(20, shape, rate, epsilon); // Draw mass points for first location
  
  if(numDraws > numMassPts){
    while(numDraws > numMassPts){
      massPts.resize(2*numMassPts);
      numMassPts = 2*numMassPts;
      Rcout << "Resized!" << std::endl;
    }
  }
  
  // Assign draws for first location to point in space
  for(j = 0; j < numDraws; j++){
    massPts[j][0] = levyDraws[j];
    massPts[j][1] = X[0] + R::runif(-1*radius, radius);
  }
  totalDraws += numDraws;
  
  // Loop through remaining locations, drawing mass points and assigning 
  // them to points in space not yet accounted for by previous locations
  for(i = 1; i < nLoc; i++){
    Rcout << "i = " << i-1 << ", " << "total mass points = " << totalDraws << std::endl;
    levyDraws = drawGamma(20, shape, rate, epsilon);

    Rcout << "1";
    if((totalDraws + numDraws) > numMassPts){
      while((totalDraws + numDraws) > numMassPts){
        massPts.resize(2*numMassPts);
        numMassPts = 2*numMassPts;
        Rcout << "Resized!" << std::endl;
      }
    }
    Rcout << "2";    
    // Assign mass point to location in space
    for(j = 0; j < numDraws; j++){
      tempLoc = X[i] + R::runif(-1*radius, radius);
      validLoc = TRUE;
      Rcout << "3";
      // Check previous locations; only need to check if current location is 
      // within 2*radius of previous location
      for(k = 0; (k<i) & validLoc; k++){
        if(distMatrix(i,k) < 2*radius){
          validLoc = validLoc & (std::abs(X[k] - tempLoc) > radius); 
        }
      }
      Rcout << "4";
      // Current draw not in any previous location's "shape"
      if(validLoc){
        massPts[totalDraws][0] = levyDraws[j];
        massPts[totalDraws++][1] = tempLoc;
      }
    }
    Rcout << "5";
    Rcout << "total Draws = " << totalDraws << std::endl;
  }
  
  // Rcout << "total draws = " << totalDraws << std::endl;

  // Calculate value of random variable at each location
  for(i = 0; i < nLoc; i++){
    Rcout << "6";
    for(j = 0; j < totalDraws; j++){
      Rcout << "7";
      if(std::abs(X[i] - massPts[j][1]) < radius){
        Rcout << "8";
        Y(i) += massPts[j][0];
        // Rcout << "X = " << X[i] << ", " << "mass = " <<
        //   massPts[j][0] << ", " << "Y = " << Y(i) << std::endl;
      }
    }    
  }
  
  // Create NumericMatrix holding mass points
  NumericMatrix mp(totalDraws, 2);
  for(i = 0; i < totalDraws; i++){
    mp(i,0) = massPts[i][0];
    mp(i,1) = massPts[i][1];
  }    
  
  //return(wrap(Y));
  return(List::create(Named("massPts") = mp, 
                      Named("Y") = Y));
}

// SEXP drawGammaRF_2D(SEXP X_,
//                     SEXP shape_,
//                     SEXP rate_,
//                     SEXP epsilon_,
//                     SEXP radius_,
//                     SEXP type_,
//                     SEXP distMatrix_){
//   // Variable definitions
//   std::vector<double> levyDraws;
//   std::vector<std::vector<double> > massPts;
//   
//   // Input data
//   // Convert SEXP objects to Rcpp objects
//   NumericMatrix X = as<NumericMatrix>(X_);
//   nLoc = X.nrow();
//   dimX = 1;
//   
//   shape   = as<double>(shape_);
//   rate    = as<double>(rate_);
//   epsilon = as<double>(epsilon_);
//   radius  = as<double>(radius_);
//   type    = as<int>(type_);
//   
//   numMassPts = 1000; // initial # of mass points needed for all locations; grow as necessary
//   
//   dimSpatial = dimX + 1; // Number of spatial dimensions + 1 for Levy measure mass value
//   
//   massPts = std::vector<std::vector<double> > (numMassPts, std::vector<double>(dimSpatial));
//   
//   GetRNGstate();
//   
//   levyDraws = drawGamma(20, epsilon, shape, rate); // Draw mass points for first location
//   numDraws = levyDraws.size();
//   
//   // Assign draws to a location
//   
//   
//   // Store results
//   NumericMatrix result(numMassPts, dimSpatial);
//   for(int i = 0; i < numMassPts; i++){
//     for(int j = 0; j < dimSpatial; j++){
//       result(i,j) = massPts[i][j];  
//     }
//   }
//   return(wrap(X));
// }

// 
// // [[Rcpp::export]]
// SEXP drawRF_Gamma1D(SEXP X_, 
//                     SEXP shape_, 
//                     SEXP rate_, 
//                     SEXP epsilon_, 
//                     SEXP radius_)
// {
//   NumericVector X = as<NumericVector>(X_);
//   NumericVector numMassPts = NumericVector(X.nrow());
//   NumericVector levyDraws;
//   IntegerVector locToStore;
//   double shape    = as<double>(shape_);
//   double rate     = as<double>(rate_);
//   double epsilon  = as<double>(epsilon_);
//   double radius   = as<double>(radius_);
//   double nu_plus;
//   NumericVector draw_radius;
//   NumericVector draw_angle;
//   int totalMassPts = 0;
//   int num_draws;
//   int i;
//   
//   nu_plus = shape*E1(rate*epsilon);
//   
//   NumericMatrix massPts = mat(10*round(nu_plus*X.nrow()), 2, fill::zeros);
//   
//   //Draw mass points for first location
//   levyDraws = drawGamma(epsilon, shape, rate, nu_plus);
//   
//   num_draws = levyDraws.size();
//   draw_radius = runif(num_draws, 0.0, radius);
//   for(i = 0; i < num_draws; i++){
//     draw_radius[i] = sqrt(draw_radius[i]);
//   }
//   draw_angle  = runif(num_draws, 0.0, 2*PI);
// 
//   locToStore = seq(totalMassPts + 1, totalMassPts + num_draws);
//   massPts[locToStore, 1] = levyDraws;
//   massPts[locToStore, 2] = draw.radius*cos(draw.angle) + X[1];
// 
//   numAdded[1] <- num.draws
//   total.mass.pts <- total.mass.pts + num.draws
//     
//     for(i in 2:nrow(X)){
// # Draw mass points for other locations, but discard points that fall 
// # in the regions of previous locations
//       levy.draws <- drawGammaLM(1e-5, 10, 1, norm.constant)
//       num.draws <- length(levy.draws)
//       draw.radius <- sqrt(runif(num.draws, 0, R))
//       draw.angle <- runif(num.draws, 0, 2*pi)
//       
//       x.loc <- X[i,1] + draw.radius*cos(draw.angle)
//       y.loc <- X[i,2] + draw.radius*sin(draw.angle)
//       
//       valid.new.pt <- rep(TRUE, num.draws)
//       
//       for(j in 1:(i-1)){
//         valid.new.pt <- valid.new.pt & !(((x.loc-X[j,1])^2 + (y.loc-X[j,2])^2) < R^2)
//       }
//       
//       num.new <- sum(valid.new.pt)
//         num.added[i] <- num.new
//         
//         if(num.new > 0){
//           mass.pts[(total.mass.pts+1):(total.mass.pts+num.new), 1] <- levy.draws[valid.new.pt]
//           mass.pts[(total.mass.pts+1):(total.mass.pts+num.new), 2] <- x.loc[valid.new.pt]
//           mass.pts[(total.mass.pts+1):(total.mass.pts+num.new), 3] <- y.loc[valid.new.pt]
//           total.mass.pts <- total.mass.pts + num.new
//         }
//     }
//     Y <- rep(0, nrow(X))
//       for(i in 1:nrow(X)){
//         Y[i] <- sum(mass.pts[(mass.pts[,2]-X[i,1])^2 + (mass.pts[,3]-X[i,2])^2 < R^2,1])
//       }
//       
//       return(cbind(X, Y, num.added))
// }