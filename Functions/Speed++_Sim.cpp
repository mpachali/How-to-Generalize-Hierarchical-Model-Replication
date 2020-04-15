// [[Rcpp::depends(RcppArmadillo)]]
#include "bayesm.h"
#include "RcppArmadillo.h"
#include <stdio.h>
#include <time.h>
using namespace arma; // use the Armadillo library for matrix computations
using namespace Rcpp;


// [[Rcpp::export]]
mat startobeta_HP_LLMns(mat const& betastar){
  
  // by Max Pachali & Thomas Otter (2016-06-11)
  
  // converts normally distributed betastars to constrained betas
  
  int nvar = betastar.n_cols;
  int draws = betastar.n_rows;
  
  mat beta = zeros(draws,nvar);
  //First Att
  //high rel. to low
  beta(span(),0) = exp(betastar(span(),0));
  //Price
  beta(span(),1) = -exp(betastar(span(),1));
  //Unrestricted Atts.
  beta(span(),span(2,3)) = betastar(span(),span(2,3)); 
  return beta;
}

// [[Rcpp::export]]
cube startobeta_cube(cube const& betastar){
  
  // by Max Pachali & Thomas Otter (2016-06-11)
  
  // converts normally distributed betastars to constrained betas
  
  int N = betastar.n_rows;
  int nvar = betastar.n_cols;
  int draws = betastar.n_slices;
  cube beta = zeros(N,nvar,draws);
  
  //First Attribute
  //high rel. to low
  beta(span(),span(0,0),span()) = exp(betastar(span(),span(0,0),span()));
  //Price
  beta(span(),span(1,1),span()) = -exp(betastar(span(),span(1,1),span()));
  //Unrestricted Atts. 
  beta(span(),span(2,3),span()) = betastar(span(),span(2,3),span());
  return beta;
}

//Posterior Means
// [[Rcpp::export]]
Mat<double> Product_Probs_PM_cpp(mat const& beta, List const& Hdata, int N, int H, int tasks, int p){
  //Convert List to std::vector of struct
  List Hdatai;
  std::vector<moments_H> Hdata_vector;
  moments_H Hdatai_struct;
  for (int h = 0; h<H; h++){
    Hdatai = Hdata[h];
    Hdatai_struct.y = as<vec>(Hdatai["y"]);
    Hdatai_struct.X = as<mat>(Hdatai["X"]);
    Hdata_vector.push_back(Hdatai_struct);    
  }
  Cube<double> Xbeta_H = zeros(N,tasks*p,H);
  //Compute Xbeta now...
  for (int h = 0; h<H; h++){
    Xbeta_H.slice(h) = beta * trans(Hdata_vector[h].X); 
  }
  //Compute Prob_chosen now...
  Cube<double> Prob_chosen = zeros(N,tasks,H);
  vec choices_y_H = zeros(tasks);
  //Main loop now...
  for(int h = 0; h<H; h++) {
    for(int k = 0; k<tasks; k++) {
      for(int r = 0; r<N; r++) {
        Cube<double> Xbeta_temp_cube = zeros(1,p,1);
        Xbeta_temp_cube(span(0,0),span(),span(0,0)) = Xbeta_H(span(r,r),span((k*p),(k*p+p-1)),span(h,h));
        vec Xbeta_temp = zeros(p);
        Xbeta_temp(span()) = trans(Xbeta_temp_cube.slice(0).cols(0,p-1));
        vec Xbeta_temp_stab = zeros(p);
        vec max_temp = zeros(p);
        max_temp.fill(max(Xbeta_temp)); 
        Xbeta_temp_stab = Xbeta_temp - max_temp;
        choices_y_H = Hdata_vector[h].y;
        double numerator = exp(Xbeta_temp_stab(choices_y_H(k)-1));
        double denominator = sum(exp(Xbeta_temp_stab));
        Prob_chosen(span(r,r),span(k,k),span(h,h)) = numerator/denominator;
      }
    }
  }
  //Compute Prob_prod now...
  Mat<double> Prob_prod = zeros(N,H);
  for(int r=0; r<N; r++){
    for(int h=0; h<H; h++){
      Prob_prod(r,h) = prod(Prob_chosen.slice(h).row(r));
    } 
  }
  return Prob_prod;
}

//Stabilized Version on log-scale
// [[Rcpp::export]]
Mat<double> Product_Probs_PM_stab_cpp(mat const& beta, List const& Hdata, int N, int H, int tasks, int p){
  //Convert List to std::vector of struct
  List Hdatai;
  std::vector<moments_H> Hdata_vector;
  moments_H Hdatai_struct;
  for (int h = 0; h<H; h++){
    Hdatai = Hdata[h];
    Hdatai_struct.y = as<vec>(Hdatai["y"]);
    Hdatai_struct.X = as<mat>(Hdatai["X"]);
    Hdata_vector.push_back(Hdatai_struct);    
  }
  Cube<double> Xbeta_H = zeros(N,tasks*p,H);
  //Compute Xbeta now...
  for (int h = 0; h<H; h++){
    Xbeta_H.slice(h) = beta * trans(Hdata_vector[h].X); 
  }
  //Compute Prob_chosen on the log-space now...
  Cube<double> Prob_chosen = zeros(N,tasks,H);
  vec choices_y_H = zeros(tasks);
  //Main loop now...
  for(int h = 0; h<H; h++) {
    for(int k = 0; k<tasks; k++) {
      for(int r = 0; r<N; r++) {
        Cube<double> Xbeta_temp_cube = zeros(1,p,1);
        Xbeta_temp_cube(span(0,0),span(),span(0,0)) = Xbeta_H(span(r,r),span((k*p),(k*p+p-1)),span(h,h));
        vec Xbeta_temp = zeros(p);
        Xbeta_temp(span()) = trans(Xbeta_temp_cube.slice(0).cols(0,p-1));
        vec Xbeta_temp_stab = zeros(p);
        vec max_temp = zeros(p);
        max_temp.fill(max(Xbeta_temp)); 
        Xbeta_temp_stab = Xbeta_temp - max_temp;
        choices_y_H = Hdata_vector[h].y;
        double numerator = Xbeta_temp_stab(choices_y_H(k)-1);
        double denominator = log(sum(exp(Xbeta_temp_stab)));
        Prob_chosen(span(r,r),span(k,k),span(h,h)) = numerator-denominator;
      }
    }
  }
  //Compute sum of log probs now...
  Mat<double> Prob_Sum_log = zeros(N,H);
  int c = 1; 
  for(int r=0; r<N; r++){
    for(int h=0; h<H; h++){
      double log_prob_sum = sum(Prob_chosen.slice(h).row(r));
      Prob_Sum_log(r,h) = log_prob_sum + c;
    } 
  }
  return Prob_Sum_log;
}

//Hierarchical Prior & Lower Level Model n.s.
//[[Rcpp::export]]
Mat<double> Product_Probs_draws_cpp(mat const& beta, int R, int l, List const& Hdata, int H, int tasks, int p){
  //Convert List to std::vector of struct
  List Hdatai;
  std::vector<moments_H> Hdata_vector;
  moments_H Hdatai_struct;
  for (int h = 0; h<H; h++){
    Hdatai = Hdata[h];
    Hdatai_struct.y = as<vec>(Hdatai["y"]);
    Hdatai_struct.X = as<mat>(Hdatai["X"]);
    Hdata_vector.push_back(Hdatai_struct);    
  }
  int draws = R*l;
  Cube<double> Xbeta_H = zeros(draws,tasks*p,H);
  //Compute Xbeta now...
  for (int h = 0; h<H; h++){
    Xbeta_H.slice(h) = beta * trans(Hdata_vector[h].X); 
  }
  //Compute Prob_chosen now...
  Cube<double> Prob_chosen = zeros(draws,tasks,H);
  vec choices_y_H = zeros(tasks);
  //Main loop now...
  for(int h = 0; h<H; h++) {
    for(int k = 0; k<tasks; k++) {
      for(int r = 0; r<draws; r++) {
        Cube<double> Xbeta_temp_cube = zeros(1,p,1);
        Xbeta_temp_cube(span(0,0),span(),span(0,0)) = Xbeta_H(span(r,r),span((k*p),(k*p+p-1)),span(h,h));
        vec Xbeta_temp = zeros(p);
        Xbeta_temp(span()) = trans(Xbeta_temp_cube.slice(0).cols(0,p-1));
        vec Xbeta_temp_stab = zeros(p);
        vec max_temp = zeros(p);
        max_temp.fill(max(Xbeta_temp)); 
        Xbeta_temp_stab = Xbeta_temp - max_temp;
        choices_y_H = Hdata_vector[h].y;
        double numerator = exp(Xbeta_temp_stab(choices_y_H(k)-1));
        double denominator = sum(exp(Xbeta_temp_stab));
        Prob_chosen(span(r,r),span(k,k),span(h,h)) = numerator/denominator;
      }
    }
  }
  //Compute Prob_prod now...
  Mat<double> Prob_prod = zeros(draws,H);
  for(int r=0; r<draws; r++){
    for(int h=0; h<H; h++){
      Prob_prod(r,h) = prod(Prob_chosen.slice(h).row(r));
    } 
  }
  return Prob_prod;
}
    
//Log-stabilized version
//[[Rcpp::export]]
Mat<double> Product_Probs_draws_stab_cpp(mat const& beta, int R, int l, List const& Hdata, int H, int tasks, int p){
  //Convert List to std::vector of struct
  List Hdatai;
  std::vector<moments_H> Hdata_vector;
  moments_H Hdatai_struct;
  for (int h = 0; h<H; h++){
    Hdatai = Hdata[h];
    Hdatai_struct.y = as<vec>(Hdatai["y"]);
    Hdatai_struct.X = as<mat>(Hdatai["X"]);
    Hdata_vector.push_back(Hdatai_struct);    
  }
  int draws = R*l;
  Cube<double> Xbeta_H = zeros(draws,tasks*p,H);
  //Compute Xbeta now...
  for (int h = 0; h<H; h++){
    Xbeta_H.slice(h) = beta * trans(Hdata_vector[h].X); 
  }
  //Compute Prob_chosen now...
  Cube<double> Prob_chosen = zeros(draws,tasks,H);
  vec choices_y_H = zeros(tasks);
  //Main loop now...
  for(int h = 0; h<H; h++) {
    for(int k = 0; k<tasks; k++) {
      for(int r = 0; r<draws; r++) {
        Cube<double> Xbeta_temp_cube = zeros(1,p,1);
        Xbeta_temp_cube(span(0,0),span(),span(0,0)) = Xbeta_H(span(r,r),span((k*p),(k*p+p-1)),span(h,h));
        vec Xbeta_temp = zeros(p);
        Xbeta_temp(span()) = trans(Xbeta_temp_cube.slice(0).cols(0,p-1));
        vec Xbeta_temp_stab = zeros(p);
        vec max_temp = zeros(p);
        max_temp.fill(max(Xbeta_temp)); 
        Xbeta_temp_stab = Xbeta_temp - max_temp;
        choices_y_H = Hdata_vector[h].y;
        double numerator = Xbeta_temp_stab(choices_y_H(k)-1);
        double denominator = log(sum(exp(Xbeta_temp_stab)));
        Prob_chosen(span(r,r),span(k,k),span(h,h)) = numerator-denominator;
      }
    }
  }
  //Compute sum of log probs now...
  Mat<double> Prob_Sum_log = zeros(draws,H);
  int c = 1; 
  for(int r=0; r<draws; r++){
    for(int h=0; h<H; h++){
      double log_prob_sum = sum(Prob_chosen.slice(h).row(r));
      Prob_Sum_log(r,h) = log_prob_sum + c;
    } 
  }
  return Prob_Sum_log;
}
    
//Market Share Computations
//[[Rcpp::export]]
vec exp_market_share_cpp(mat const& beta, cube const& design_joint){
  int draws = beta.n_rows;
  int pa = design_joint.n_slices;
  int nplay = design_joint.n_rows;
  Cube<double> Xbeta = zeros(draws,nplay,pa);
  //Xbeta of all product possibilities
  for(int p=0; p<pa; p++){
    mat sub_design = design_joint.slice(p);   
    Xbeta.slice(p) = beta * trans(sub_design);
  }
  Mat<double> market_share_draws = zeros(draws,pa);
  //Compute stabilized choice probs now...
  for(int p = 0; p<pa; p++){
    for(int r = 0; r<draws; r++){
      vec Xbeta_temp = zeros(nplay);
      Xbeta_temp = trans(Xbeta.slice(p).row(r));
      vec Xbeta_temp_stab = zeros(nplay);
      vec max_temp = zeros(nplay);
      max_temp.fill(max(Xbeta_temp)); 
      Xbeta_temp_stab = Xbeta_temp - max_temp;
      double numerator = exp(Xbeta_temp_stab(0));
      double denominator = sum(exp(Xbeta_temp_stab));
      market_share_draws(r,p) = numerator/denominator;
    }
  }
  vec exp_ms = zeros(pa);
  for(int p=0; p<pa; p++){
    vec pa_draws = market_share_draws(span(),p);
    exp_ms(p) = mean(pa_draws);
  }
  return exp_ms;
}
    
    
