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
double cone_height1D(double, double);
double E1(double);
double E1inv(double, double, bool);
std::vector<double> drawGamma(int, double, double, double);

SEXP E1(SEXP);
SEXP E1inv(SEXP, SEXP, SEXP);
SEXP drawGamma(SEXP, SEXP, SEXP, SEXP);
SEXP cone_height1D(SEXP, SEXP);

double cone_height1D(double d, double r){
  return((r-std::abs(d))/pow(r,2.0));
}

// [[Rcpp::export]]
SEXP cone_height1D(SEXP d, SEXP r){
  return(wrap(cone_height1D(as<double>(d), as<double>(r))));
}

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
  int i, j, k, ind = 0;
  double x_inc, x_star, h_star;
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
      massPts.reserve(2*numMassPts);
      numMassPts = 2*numMassPts;
    }
  }
  
  // Assign Levy draws to each point in space
  if(type==1){
    for(j = 0; j < numDraws; j++){
      massPts.push_back(std::vector<double>(dimSpatial));
      massPts[j][0] = levyDraws[j];
      x_inc = R::runif(-1*radius, radius);
      massPts[j][1] = X[0] + x_inc;
    }
    totalDraws += numDraws;
    
    // Loop through remaining locations, drawing mass points and assigning 
    // them to points in space not yet accounted for by previous locations
    for(i = 1; i < nLoc; i++){
      levyDraws = drawGamma(20, shape, rate, epsilon);
      
      if((totalDraws + numDraws) > numMassPts){
        while((totalDraws + numDraws) > numMassPts){
          massPts.reserve(2*numMassPts);
          numMassPts = 2*numMassPts;
        }
      }
      
      // Assign mass point to location in space
      for(j = 0; j < numDraws; j++){
        x_inc = R::runif(-1*radius, radius);
        x_star = X[i] + x_inc;
        validLoc = TRUE;
        // Check previous locations; only need to check if current location is 
        // within 2*radius of previous location
        for(k = 0; (k<i) & validLoc; k++){
          if(distMatrix(i,k) < 2*radius){
            validLoc = validLoc & (std::abs(X[k] - x_star) > radius); 
          }
        }
        // Current draw not in any previous location's "shape"
        if(validLoc){
          massPts.push_back(std::vector<double>(dimSpatial));
          massPts[totalDraws][0] = levyDraws[j];
          massPts[totalDraws++][1] = x_star;
        }
      }
    }

    // Calculate value of random variable at each location
    for(i = 0; i < nLoc; i++){
      for(j = 0; j < totalDraws; j++){
        if(std::abs(X[i] - massPts[j][1]) < radius){
          Y(i) += massPts[j][0];
        }
      }    
    }
  }else if(type==2){
    for(j = 0; j < numDraws; j++){
      massPts.push_back(std::vector<double>(dimSpatial));
      ind = std::floor(R::rbinom(1.0, 0.5) + 0.5);
      
      if(ind == 0){
        x_star = (X[0] - radius) + radius*R::rbeta(2.0, 1.0);
      }else{
        x_star = X[0] + radius*R::rbeta(1.0, 2.0);
      }
      h_star = R::runif(0, cone_height1D(x_star - X[0], radius));
      
      massPts[j][0] = levyDraws[j];
      massPts[j][1] = x_star;
      massPts[j][2] = h_star;
    }
    totalDraws += numDraws;
    
    // Loop through remaining locations, drawing mass points and assigning 
    // them to points in space not yet accounted for by previous locations
    for(i = 1; i < nLoc; i++){
      levyDraws = drawGamma(20, shape, rate, epsilon);
      
      if((totalDraws + numDraws) > numMassPts){
        while((totalDraws + numDraws) > numMassPts){
          massPts.reserve(2*numMassPts);
          numMassPts = 2*numMassPts;
        }
      }
      
      // Assign mass point to location in space
      for(j = 0; j < numDraws; j++){
        ind = std::floor(R::rbinom(1.0, 0.5) + 0.5);
        
        if(ind == 0){
          x_star = (X[i] - radius) + radius*R::rbeta(2.0, 1.0);
        }else{
          x_star = X[i] + radius*R::rbeta(1.0, 2.0);
        }
        h_star = R::runif(0, cone_height1D(x_star - X[i], radius));
        
        validLoc = TRUE;
        // Check previous locations; only need to check if current location is 
        // within 2*radius of previous location
        for(k = 0; (k<i) & validLoc; k++){
          if(distMatrix(i,k) < 2*radius){
            validLoc = validLoc & 
              (h_star > cone_height1D(x_star - X[k], radius));
          }
        }
        
        // Current draw not in any previous location's "shape"
        if(validLoc){
          massPts.push_back(std::vector<double>(dimSpatial));
          massPts[totalDraws][0] = levyDraws[j];
          massPts[totalDraws][1] = x_star;
          massPts[totalDraws++][2] = h_star;
        }
      }
    }
    
    // Calculate value of random variable at each location
    for(i = 0; i < nLoc; i++){
      for(j = 0; j < totalDraws; j++){
        if(massPts[j][2] < cone_height1D(X[i] - massPts[j][1], radius)){
          Y(i) += massPts[j][0];
        }
      }    
    }
  }
  
  return(wrap(Y));
}

//     // [[Rcpp::export]]
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
//   int i, j, k, ind = 0;
//   double draw_radius, draw_angle, x_star, y_star, h_star;
//   bool validLoc;
//   NumericVector Y;
//   
//   // Input data
//   // Convert SEXP objects to Rcpp objects
//   NumericMatrix X = as<NumericMatrix>(X_);
//   nLoc = X.nrow();
//   dimX = 2;
//   dimSpatial = dimX + 2; // Number of spatial dimensions + 2 (for Levy mass value & "height")
//   
//   NumericMatrix distMatrix = as<NumericMatrix>(distMatrix_);
//   
//   shape   = as<double>(shape_);
//   rate    = as<double>(rate_);
//   epsilon = as<double>(epsilon_);
//   radius  = as<double>(radius_);
//   type    = as<int>(type_);
//   
//   numMassPts = nLoc*100; // initial # of mass points needed for all locations; grow as necessary
//   totalDraws = 0;
//   
//   Y = NumericVector(nLoc);
//   
//   massPts = std::vector<std::vector<double> > (numMassPts, std::vector<double>(dimSpatial));
//   
//   GetRNGstate();
//   
//   levyDraws = drawGamma(20, shape, rate, epsilon); // Draw mass points for first location
//   
//   if(numDraws > numMassPts){
//     while(numDraws > numMassPts){
//       massPts.reserve(2*numMassPts);
//       numMassPts = 2*numMassPts;
//     }
//   }
//   
//   // Assign Levy draws to each point in space
//   if(type==1){
//     for(j = 0; j < numDraws; j++){
//       massPts.push_back(std::vector<double>(dimSpatial));
//       massPts[j][0] = levyDraws[j];
//       
//       draw_radius = std::sqrt(R::runif(0, radius));
//       draw_angle  = R::runif(0, M_2PI);
//       x_star = X(0,0) + draw_radius*cos(draw_angle);
//       y_star = X(0,1) + draw_radius*sin(draw_angle);
//         
//       massPts[j][1] = x_star;
//       massPts[j][2] = y_star;
//     }
//     totalDraws += numDraws;
//     
//     // Loop through remaining locations, drawing mass points and assigning 
//     // them to points in space not yet accounted for by previous locations
//     for(i = 1; i < nLoc; i++){
//       levyDraws = drawGamma(20, shape, rate, epsilon);
//       
//       if((totalDraws + numDraws) > numMassPts){
//         while((totalDraws + numDraws) > numMassPts){
//           massPts.reserve(2*numMassPts);
//           numMassPts = 2*numMassPts;
//         }
//       }
//       
//       // Assign mass point to location in space
//       for(j = 0; j < numDraws; j++){
//         draw_radius = std::sqrt(R::runif(0, radius));
//         draw_angle  = R::runif(0, M_2PI);
//         x_star = X(i,0) + draw_radius*cos(draw_angle);
//         y_star = X(i,1) + draw_radius*sin(draw_angle);
//           
//         validLoc = TRUE;
//         // Check previous locations; only need to check if current location is 
//         // within 2*radius of previous location
//         for(k = 0; (k<i) & validLoc; k++){
//           if(distMatrix(i,k) < 2*radius){
//             // validLoc = validLoc & (std::abs(X[k] - x_star) > radius); 
//             validLoc = validLoc & 
//               (((X(k,0)-x_star)*(X(k,0)-x_star) + (X(k,1)-y_star)*(X(k,1)-y_star)) > 
//               radius*radius);
//           }
//         }
//         // Current draw not in any previous location's "shape"
//         if(validLoc){
//           massPts.push_back(std::vector<double>(dimSpatial));
//           massPts[totalDraws][0] = levyDraws[j];
//           massPts[totalDraws][1] = x_star;
//           massPts[totalDraws++][2] = y_star;
//         }
//       }
//     }
//     
//     // Calculate value of random variable at each location
//     for(i = 0; i < nLoc; i++){
//       for(j = 0; j < totalDraws; j++){
//         // if(std::abs(X[i] - massPts[j][1]) < radius){
//         if((((X(i,0)-massPts[j][1])*(X(i,0)-massPts[j][1]) + 
//            (X(i,1)-massPts[j][2])*(X(i,1)-massPts[j][2]))) < (radius*radius)){
//           Y(i) += massPts[j][0];
//         }
//       }    
//     }
//   }else if(type==2){
//     for(j = 0; j < numDraws; j++){
//       massPts.push_back(std::vector<double>(dimSpatial));
//       ind = std::floor(R::rbinom(1.0, 0.5) + 0.5);
//       
//       if(ind == 0){
//         x_star = (X[0] - radius) + radius*R::rbeta(2.0, 1.0);
//       }else{
//         x_star = X[0] + radius*R::rbeta(1.0, 2.0);
//       }
//       h_star = R::runif(0, cone_height1D(x_star - X[0], radius));
//       
//       massPts[j][0] = levyDraws[j];
//       massPts[j][1] = x_star;
//       massPts[j][2] = h_star;
//     }
//     totalDraws += numDraws;
//     
//     // Loop through remaining locations, drawing mass points and assigning 
//     // them to points in space not yet accounted for by previous locations
//     for(i = 1; i < nLoc; i++){
//       levyDraws = drawGamma(20, shape, rate, epsilon);
//       
//       if((totalDraws + numDraws) > numMassPts){
//         while((totalDraws + numDraws) > numMassPts){
//           massPts.reserve(2*numMassPts);
//           numMassPts = 2*numMassPts;
//         }
//       }
//       
//       // Assign mass point to location in space
//       for(j = 0; j < numDraws; j++){
//         ind = std::floor(R::rbinom(1.0, 0.5) + 0.5);
//         
//         if(ind == 0){
//           x_star = (X[i] - radius) + radius*R::rbeta(2.0, 1.0);
//         }else{
//           x_star = X[i] + radius*R::rbeta(1.0, 2.0);
//         }
//         h_star = R::runif(0, cone_height1D(x_star - X[i], radius));
//         
//         validLoc = TRUE;
//         // Check previous locations; only need to check if current location is 
//         // within 2*radius of previous location
//         for(k = 0; (k<i) & validLoc; k++){
//           if(distMatrix(i,k) < 2*radius){
//             validLoc = validLoc & 
//               (h_star > cone_height1D(x_star - X[k], radius));
//           }
//         }
//         
//         // Current draw not in any previous location's "shape"
//         if(validLoc){
//           massPts.push_back(std::vector<double>(dimSpatial));
//           massPts[totalDraws][0] = levyDraws[j];
//           massPts[totalDraws][1] = x_star;
//           massPts[totalDraws++][2] = h_star;
//         }
//       }
//     }
//     
//     // Calculate value of random variable at each location
//     for(i = 0; i < nLoc; i++){
//       for(j = 0; j < totalDraws; j++){
//         if(massPts[j][2] < cone_height1D(X[i] - massPts[j][1], radius)){
//           Y(i) += massPts[j][0];
//         }
//       }    
//     }
//   }
//   
//   return(wrap(Y));
// }