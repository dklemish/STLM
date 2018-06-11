// // [[Rcpp::plugins(cpp11)]]
// 
// // Declare dependencies so Rcpp knows to link libraries
// // [[Rcpp::depends(RcppGSL)]]
// 
// #include <R.h>
// #include <Rcpp.h>
// #include <Rmath.h>
// #include <RcppGSL.h>
// #include <gsl/gsl_blas.h>
// #include <gsl/gsl_rng.h>
// #include <gsl/gsl_randist.h>
// #include <gsl/gsl_linalg.h>
// #include <chrono>
// 
// using namespace Rcpp;
// 
// /***** Function headers ******/
// // Random number generation
// SEXP rmvt(SEXP, SEXP, SEXP);
// SEXP rmvnorm(SEXP, SEXP);
// SEXP rmvt_chol(SEXP, SEXP, SEXP);
// SEXP rmvnorm_chol(SEXP, SEXP);
// // SEXP rmvt_chol_R(SEXP, SEXP, SEXP);
// // 
// RcppGSL::vector<double> rmvt(RcppGSL::vector<double>,
//                              RcppGSL::matrix<double>,
//                              double,
//                              gsl_rng *);
// 
// RcppGSL::vector<double> rmvnorm(RcppGSL::vector<double>,
//                                 RcppGSL::matrix<double>,
//                                 gsl_rng *);
// 
// RcppGSL::vector<double> rmvt_chol(RcppGSL::vector<double>, 
//                                   RcppGSL::matrix<double>, 
//                                   double,
//                                   gsl_rng *);
// 
// RcppGSL::vector<double> rmvnorm_chol(RcppGSL::vector<double>,
//                                      RcppGSL::matrix<double>,
//                                      gsl_rng *);
// 
// // RcppGSL::vector<double> rmvt_chol_R(RcppGSL::vector<double>, 
// //                                     RcppGSL::matrix<double>, 
// //                                     double);
// 
// // Density functions
// // SEXP dmvt(SEXP, SEXP, SEXP, SEXP);
// // SEXP dmvnorm(SEXP, SEXP, SEXP, SEXP);
// // SEXP dmvt_chol(SEXP, SEXP, SEXP, SEXP);
// // SEXP dmvnorm_chol(SEXP, SEXP, SEXP, SEXP);
// // 
// // double dmvt(RcppGSL::vector<double>, 
// //             RcppGSL::vector<double>, 
// //             RcppGSL::matrix<double>, 
// //             double);
// // 
// // double dmvnorm(RcppGSL::vector<double>, 
// //                RcppGSL::vector<double>, 
// //                RcppGSL::matrix<double>);
// // 
// // double dmvt_chol(RcppGSL::vector<double>, 
// //                  RcppGSL::vector<double>, 
// //                  RcppGSL::matrix<double>, 
// //                  double);
// // 
// // double dmvnorm_chol(RcppGSL::vector<double>, 
// //                     RcppGSL::vector<double>, 
// //                     RcppGSL::matrix<double>);
// 
// /***** End function headers *****/
// 
// 
// /***** Function definitions *****/
// gsl_rng * set_seed(void){
//   const gsl_rng_type * T;
//   gsl_rng * r;
//   
//   gsl_rng_env_setup();
//   T = gsl_rng_default;
//   r = gsl_rng_alloc(T);
//   
//   // Set seed to random time
//   int time = static_cast<long int>(std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count());
//   gsl_rng_set(r, time);
//   
//   return(r);
// }
// 
// // [[Rcpp::export]]
// SEXP rmvt_c(SEXP mu_, SEXP Sigma_, SEXP df){
//   gsl_rng * r = set_seed();
//   RcppGSL::vector<double> mu = mu_;
//   RcppGSL::matrix<double> Sigma = Sigma_;
//   
//   return(wrap(rmvt(mu, Sigma, as<double>(df), r)));
// }
// 
// // [[Rcpp::export]]
// SEXP rmvnorm_c(SEXP mu_, SEXP Sigma_){
//   gsl_rng * r = set_seed();
//   RcppGSL::vector<double> mu = mu_;
//   RcppGSL::matrix<double> Sigma = Sigma_;
//   
//   return(wrap(rmvnorm(mu, Sigma, r)));
// }
// 
// // [[Rcpp::export]]
// SEXP rmvt_chol(SEXP mu_, SEXP L_, SEXP df){
//   gsl_rng * r = set_seed();
//   RcppGSL::vector<double> mu = mu_;
//   RcppGSL::matrix<double> L  = L_;
//   
//   return(wrap(rmvt_chol(mu, L, as<double>(df), r)));
// }
// 
// // [[Rcpp::export]]
// SEXP rmvnorm_chol(SEXP mu_, SEXP L_){
//   gsl_rng * r = set_seed();
//   RcppGSL::vector<double> mu = mu_;
//   RcppGSL::matrix<double> L  = L_;
//   
//   return(wrap(rmvnorm_chol(mu, L, r)));
// }
// 
// // // [[Rcpp::export]]
// // SEXP dmvt(SEXP x_, SEXP mu_, SEXP Sigma_, SEXP df_){
// //   RcppGSL::vector<double> x     = x_;
// //   RcppGSL::vector<double> mu    = mu_;
// //   RcppGSL::matrix<double> Sigma = Sigma_;
// //   
// //   return(wrap(rmvt(x, mu, Sigma, as<double>(df))));
// // }
// // 
// // // [[Rcpp::export]]
// // SEXP dmvnorm(SEXP x_, SEXP mu_, SEXP Sigma_){
// //   RcppGSL::vector<double> x     = x_;
// //   RcppGSL::vector<double> mu    = mu_;
// //   RcppGSL::matrix<double> Sigma = Sigma_;
// //   
// //   return(wrap(rmvnorm(x, mu, Sigma)));
// // }
// // 
// // // [[Rcpp::export]]
// // SEXP dmvt_chol(SEXP x_, SEXP mu_, SEXP L_, SEXP df){
// //   RcppGSL::vector<double> x  = x_;
// //   RcppGSL::vector<double> mu = mu_;
// //   RcppGSL::matrix<double> L  = L_;
// // 
// //   return(wrap(rmvt_chol(x, mu, L, as<double>(df))));
// // }
// // 
// // // [[Rcpp::export]]
// // SEXP dmvnorm_chol(SEXP x_, SEXP mu_, SEXP L_){
// //   RcppGSL::vector<double> x  = x_;
// //   RcppGSL::vector<double> mu = mu_;
// //   RcppGSL::matrix<double> L  = L_;
// //   
// //   return(wrap(rmvnorm_chol(x, mum L)));
// // }
// RcppGSL::vector<double> 
//   rmvt(RcppGSL::vector<double> mu, 
//        RcppGSL::matrix<double> Sigma, 
//        double df,
//        gsl_rng *r){
//     
//     int i;
//     int d = Sigma.ncol();
//     double u;
//     RcppGSL::vector<double> normal(d);
//     RcppGSL::vector<double> result(d);
//     
//     for(i = 0; i < d; i++){
//       normal[i] = gsl_ran_gaussian(r, 1.0);
//     }
//     u = gsl_ran_chisq(r, df);
//     
//     gsl_linalg_cholesky_decomp(Sigma);
//     gsl_blas_dsymv(CblasLower, 1.0, Sigma, normal, 0.0, result);
//     gsl_vector_scale(result, 1/sqrt(u/df));
//     gsl_vector_add(result, mu);
//     
//     return(result);
//   }
// 
// RcppGSL::vector<double> 
//   rmvnorm(RcppGSL::vector<double> mu, 
//           RcppGSL::matrix<double> Sigma, 
//           gsl_rng *r){
//     
//     int i;
//     int d = Sigma.ncol();
//     RcppGSL::vector<double> normal(d);
//     RcppGSL::vector<double> result(d);
//     
//     for(i = 0; i < d; i++){
//       normal[i] = gsl_ran_gaussian(r, 1.0);
//     }
//     
//     gsl_linalg_cholesky_decomp(Sigma);
//     gsl_blas_dsymv(CblasLower, 1.0, Sigma, normal, 0.0, result);
//     gsl_vector_add(result, mu);
//     
//     return(result);
//   }
// 
// RcppGSL::vector<double> 
//   rmvt_chol(RcppGSL::vector<double> mu, 
//             RcppGSL::matrix<double> L, 
//             double df,
//             gsl_rng *r){
//     
//     int i;
//     int d = L.ncol();
//     double u;
//     RcppGSL::vector<double> normal(d);
//     RcppGSL::vector<double> result(d);
//     
//     for(i = 0; i < d; i++){
//       normal[i] = gsl_ran_gaussian(r, 1.0);
//     }
//     u = gsl_ran_chisq(r, df);
//     
//     gsl_blas_dsymv(CblasLower, 1.0, L, normal, 0.0, result);
//     gsl_vector_scale(result, 1/sqrt(u/df));
//     gsl_vector_add(result, mu);
//     
//     return(result);
//   }
// 
// RcppGSL::vector<double> 
//   rmvnorm_chol(RcppGSL::vector<double> mu, 
//                RcppGSL::matrix<double> L, 
//                gsl_rng *r){
//     
//     int i;
//     int d = L.ncol();
//     RcppGSL::vector<double> normal(d);
//     RcppGSL::vector<double> result(d);
//     
//     for(i = 0; i < d; i++){
//       normal[i] = gsl_ran_gaussian(r, 1.0);
//     }
//     
//     gsl_blas_dsymv(CblasLower, 1.0, L, normal, 0.0, result);
//     gsl_vector_add(result, mu);
//     
//     return(result);
//   }
// 
// // // [[Rcpp::export]]
// // SEXP rmvt_chol_R(SEXP mu_, SEXP L_, SEXP df){
// //   RcppGSL::vector<double> mu = mu_;
// //   RcppGSL::matrix<double> L  = L_;
// //   
// //   return(wrap(rmvt_chol_R(mu, L, as<double>(df))));
// // }
// 
// // RcppGSL::vector<double> 
// //   rmvt_chol_R(RcppGSL::vector<double> mu, 
// //             RcppGSL::matrix<double> L, 
// //             double df){
// //     int i;
// //     int d = L.ncol();
// //     double u;
// //     RcppGSL::vector<double> normal(d);
// //     RcppGSL::vector<double> result(d);
// //     
// //     for(i = 0; i < d; i++){
// //       normal[i] = R::rnorm(0.0, 1.0);
// //     }
// //     u = R::rchisq(df);
// //     
// //     gsl_blas_dsymv(CblasUpper, 1.0, L, normal, 0.0, result);
// //     gsl_vector_scale(result, 1/sqrt(u/df));
// //     gsl_vector_add(result, mu);
// //     
// //     return(result);
// //   }