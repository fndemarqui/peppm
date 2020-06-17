#include <Rcpp.h>
#include <stdio.h>
#include <math.h>
#include <Rmath.h>
#include <stdlib.h>
#include <R_ext/Random.h>
#include <R.h>
#include <Rdefines.h>
#include <vector>


using namespace Rcpp;



Rcpp::IntegerVector idblocks(Rcpp::NumericVector time,
                        Rcpp::NumericVector grid){
  int n = time.size();
  int j=0;
  int i;
  Rcpp::IntegerVector id(n);
  for (i = 0; i < n; i++){
    if(time[i] <= grid[j+1])
      id[i] = j;
    else{
      j++;
      id[i] = j;
    }
  }
  return id;
}

Rcpp::IntegerVector tableRcpp(Rcpp::IntegerVector x){
  int n = x.size();
  int b = max(x) + 1;
  Rcpp::IntegerVector nj(b);
  for (int i = 0; i < n; i++) {
    nj[x[i]]++;
  }
  return nj;
}


Rcpp::IntegerVector idrates(Rcpp::IntegerVector U){
  int m = U.size() + 1;
  int i;
  int j=0;
  
  Rcpp::IntegerVector idr(m);
  for(i=0;i<m-1; i++){
    if(U[i]==1)
      idr[i+1]=j;
    else{
      j++;
      idr[i+1]=j;
    }
  }
  return idr;
}


//' Computes the time grid from the auxiliary vector U.
//' @aliases getGrid
//' @param U vector of change point indicators
//' @param ftgrid vector with the finest time grid (distinct observed failure times)
//' @return the time grid associated with the auxiliary vector U.
//'
// [[Rcpp::export]]
Rcpp::NumericVector getGrid(Rcpp::IntegerVector U, Rcpp::NumericVector ftgrid){
//  int getGrid(Rcpp::IntegerVector U, Rcpp::NumericVector ftgrid){
  int m = U.size();
  int j=0;
  int i;
  int b = m - sum(U) + 1;
  NumericVector grid(b+1);

  j=1;
  grid[0] = 0;
  for(i=0;i<m;i++){
    if(U[i]==0){
      grid[j] = ftgrid[i];				
      j++;
    }
  }
  grid[b] = R_PosInf;
  return grid;
}



Rcpp::List suffstats(Rcpp::NumericVector time, Rcpp::IntegerVector status,
                     Rcpp::NumericVector grid, Rcpp::IntegerVector idb){

  int n = time.size();
  int b = grid.size() - 1;
  int j=0;
  int i;
  
  Rcpp::IntegerVector nu(b);
  Rcpp::NumericVector xi(b);
  
  for (j = 0; j < b; j++){
    nu[j] = 0; 
    xi[j] = 0;
    for (i = 0; i < n; i++){ 
      nu[j] += status[i]*(idb[i]==j);
      xi[j] += ( std::min(time[i], grid[j+1]) - grid[j] )*(time[i] >= grid[j]);
    }
  }
  return Rcpp :: List :: create ( Rcpp::Named("nu") = nu , Rcpp::Named("xi") = xi );
}


double loglik(Rcpp::IntegerVector U, Rcpp::NumericVector ftgrid,
                Rcpp::NumericVector time, Rcpp::IntegerVector status,
                double a_rates, double b_rates){

  int j;
  double logfatj = 0;
  Rcpp::NumericVector grid = getGrid(U, ftgrid);
  Rcpp::IntegerVector id = idblocks(time, grid);
  Rcpp::List stats = suffstats(time, status, grid, id);
  
  
  IntegerVector nu = stats["nu"];
  NumericVector xi = stats["xi"];
  
  int b = nu.size();
  for(j=0; j<b; j++){
    logfatj += a_rates*log(b_rates) - (a_rates + nu[j])*log(b_rates + xi[j])
            + lgamma(a_rates) - lgamma(a_rates + nu[j]);
  }
   return logfatj;
}


double get_logpred(Rcpp::IntegerVector U, Rcpp::NumericVector ftgrid,
                   Rcpp::NumericVector time, Rcpp::IntegerVector status,
                   double a_rates, double b_rates, int cohesion,
                   double a_beta, double b_beta){
  
  int j;
  int m = U.size() + 1;
  double logfatj = 0;
  Rcpp::NumericVector grid = getGrid(U, ftgrid);
  Rcpp::IntegerVector id = idblocks(time, grid);
  Rcpp::List stats = suffstats(time, status, grid, id);
  //Rcpp::IntegerVector nj = tableRcpp(id);
  
  IntegerVector nu = stats["nu"];
  NumericVector xi = stats["xi"];
  IntegerVector nj = stats["nu"];
  
  int b = nu.size();
  if(cohesion==1){
    for(j=0; j<b; j++){
      logfatj += a_rates*log(b_rates) - (a_rates + nu[j])*log(b_rates + xi[j]) + lgamma(a_rates) - lgamma(a_rates + nu[j]);
    }
  }else if(cohesion==2){
    for(j=0; j<b; j++){
      logfatj += a_rates*log(b_rates) - (a_rates + nu[j])*log(b_rates + xi[j]) + lgamma(a_rates) - lgamma(a_rates + nu[j]) + nj[j];
    }
  }else if(cohesion==3){
    for(j=0; j<b; j++){
      logfatj += a_rates*log(b_rates) - (a_rates + nu[j])*log(b_rates + xi[j]) + lgamma(a_rates) - lgamma(a_rates + nu[j]) +  1/nj[j];
    }
  }else{
    for(j=0; j<b; j++){
      logfatj += a_rates*log(b_rates) - (a_rates + nu[j])*log(b_rates + xi[j]) + lgamma(a_rates) - lgamma(a_rates + nu[j]) +
        lgamma(m + b_beta - b + 1) - lgamma(m + b_beta - b) + lgamma(b + a_beta - 2)  - lgamma(b + a_beta - 1);
    }
  }
  
  return logfatj;
}


Rcpp::IntegerVector samplerU(Rcpp::IntegerVector U, Rcpp::NumericVector ftgrid,
                Rcpp::NumericVector time, Rcpp::IntegerVector status,
                double a_rates, double b_rates, int cohesion,
                double a_beta, double b_beta){
  
  int r;
  int m = U.size() + 1;
  double logfat0;
  double logfat1;
  double R;
  Rcpp::NumericVector w(m-1);
  GetRNGstate();
  w = Rcpp::runif(m-1 , 0, 1 );
  for(r=0; r<m-1; r++){
    U[r]=1;
    logfat1 = get_logpred(U, ftgrid, time, status, a_rates, b_rates, cohesion, a_beta, b_beta);
    U[r]=0;
    logfat0 = get_logpred(U, ftgrid, time, status, a_rates, b_rates, cohesion, a_beta, b_beta);
    R = exp(logfat0 - logfat1);
    if(R > (1-w[r])/w[r]){
      U[r] = 1;
    }else{
      U[r] = 0;
    }
  }
  return U;
}

Rcpp::NumericVector samplerRates(Rcpp::IntegerVector U, Rcpp::NumericVector ftgrid,
                             Rcpp::NumericVector time, Rcpp::IntegerVector status,
                             double a_rates, double b_rates){
  
  RNGScope scope;
  int j, k;
  int m = U.size() + 1;
  Rcpp::NumericVector grid = getGrid(U, ftgrid);
  Rcpp::IntegerVector idb = idblocks(time, grid);
  Rcpp::IntegerVector idr = idrates(U);
  Rcpp::List stats = suffstats(time, status, grid, idb);
  IntegerVector nu = stats["nu"];
  NumericVector xi = stats["xi"];
  
  int b = nu.size();
  NumericVector rates(m);
  NumericVector rates_aux(b);
  for(j=0; j<b; j++){
    rates_aux[j] = Rcpp::rgamma(1, a_rates + nu[j], 1/(b_rates + xi[j]) )[0];
  }
  
  j = 0;
  rates[0] = rates_aux[j];
  for(k=1;k<m;k++){
    if(idr[k-1]==idr[k]){
      rates[k]=rates_aux[j];
    }else{
      j++;
      rates[k] = rates_aux[j];
    }
  }
  return rates;
}

//' Runs the Gibbs sampler
//' @aliases gibbs
//' @param U0 vector of change point indicators
//' @param ftgrid vector of indexes of distinct failure times
//' @param time vector of observed failure times.
//' @param status vector of failure indicators
//' @param a_rates shape parameter of the gamma distribution (prior for failure rates).
//' @param b_rates scale parameter of the gamma distribution (prior for failure rates).
//' @param cohesion type of prior cohesion (1 to 4).
//' @param a_beta shape1 parameter of the beta distribution (prior for p - cohesion 4).
//' @param b_beta shape2 parameter of the beta distribution (prior for p - cohesion 4).
//' @param nburnin number of iterations to be discarded.
//' @param npost desired posterior sample size
//' @param nlag number of jumps to eliminate autocorrelation of the chain.
//' @return posterior sample
//' 
// [[Rcpp::export]]
Rcpp::List gibbs(Rcpp::IntegerVector U0, Rcpp::NumericVector ftgrid,
                             Rcpp::NumericVector time, Rcpp::IntegerVector status,
                             double a_rates, double b_rates,
                             int cohesion, double a_beta, double b_beta,
                             int npost, int nburnin, int nlag){
  
  int nsim = nburnin + npost*nlag;
  int k;
  int s = 0;
  int r = 0;
  int m  = U0.size() + 1;
  int b0;
  Rcpp::NumericMatrix rates(npost, m);
  Rcpp::NumericMatrix U(npost, m-1);
  IntegerVector b(npost);
  NumericVector lpred(npost);

  //IntegerVector b(npost);

  for(s=0; s<nburnin; s++){
    U0 = samplerU(U0, ftgrid, time, status, a_rates, b_rates, cohesion, a_beta, b_beta);
    b0 = m - sum(U0); 
  }  

  for(s=0; s<npost; s++){
    U0 = samplerU(U0, ftgrid, time, status, a_rates, b_rates, cohesion, a_beta, b_beta);
    b[s] = m - sum(U0) ; 
    Rcpp::NumericVector rates_aux = samplerRates(U0, ftgrid, time, status, a_rates, b_rates);
    for(k=0; k< m-1; k++){
      rates(s, k) = rates_aux[k];
      U(s, k) = U0[k];
    }
    rates(s, m-1) = rates_aux[m-1];
    lpred[s] = loglik(U0, ftgrid, time, status, a_rates, b_rates);
  }
  for(r=0; r<nlag; r++){
    U0 = samplerU(U0, ftgrid, time, status, a_rates, b_rates, cohesion, a_beta, b_beta);      
  }
   
  return Rcpp::List::create(        
    Rcpp::Named("b", wrap(b)),
    Rcpp::Named("rates", wrap(rates)),
    Rcpp::Named("U", wrap(U)),
    Rcpp::Named("lpred", wrap(lpred))
  );             // Return to R
                              
}
