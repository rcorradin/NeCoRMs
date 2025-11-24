#include <RcppArmadillo.h>
#include <RcppDist.h>
using namespace arma;
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

//-------------------------------------------
//-------------------------------------------
//-------------------------------------------

int rint_not_unif(vec weights){
  double u = randu();
  vec probs = weights / sum(weights);
  probs = cumsum(probs);
  
  for(uword k = 0; k < probs.n_elem; k++) {
    if(u <= probs(k)) {
      return k;
    }
  }
  return -1;
}

int rint_not_unif_log(vec lweights){
  
  double u = randu();
  vec probs(lweights);
  // double mw = max(lweights);
  // 
  // probs = exp(lweights - mw);
  for(uword j = 0; j < lweights.n_elem; j++){
    probs(j) = 1.0 / sum(exp(lweights - lweights(j)));
  }
  probs /= sum(probs);
  probs = cumsum(probs);
  
  for(uword k = 0; k < probs.n_elem; k++) {
    if(u <= probs[k]) {
      return k;
    }
  }
  return -1;
}

//-----------------------------------------------------
// Update freqs
// compute the frequencies for each group and each global active cluster
//-----------------------------------------------------

void update_freqs(int d, 
                  vec group,  
                  vec clust,
                  mat &freqs){
  int k = max(clust) + 1;
  freqs.resize(d,k);
  
  for(uword l = 0; l < k; l++){
    for(uword j = 0; j < d; j++){
      freqs(j,l) = accu(clust == l && group == j);
    }
  }
}


void update_freqs_nested(int d, 
                         vec group,  
                         vec nested_c, 
                         vec clust,
                         mat &freqs_nested){
  int k = max(clust) + 1;
  freqs_nested.resize(d,k);
  freqs_nested.fill(0.0);
  vec temp_alloc(group);
  
  for(uword i = 0; i < group.n_elem; i++){
    temp_alloc(i) = nested_c(group(i));
  }
  
  for(uword l = 0; l < k; l++){
    for(uword j = 0; j < d; j++){
      freqs_nested(j,l) = accu(clust == l && temp_alloc == j);
    }
  }
}

//-----------------------------------------------------
// Update Tlj's
// sample the group-and-cluster specific weights for the existent atoms
// the full conditional is a Gamma(n_lj + phi, sigma_l u_j + 1)
//-----------------------------------------------------

void update_Tlj_nested(mat freqs_nested,      
                       int d,
                       double phi, 
                       vec nested_c, 
                       vec u, 
                       vec sigmal,     
                       mat &Tlj){      
  
  int k = freqs_nested.n_cols;
  Tlj.resize(d, k);
  Tlj.fill(0.0);
  
  for(uword l = 0; l < k; l++){
    for(uword j = 0; j < d; j++){
      if(accu(nested_c == j) > 0){
        Tlj(j,l) = randg(distr_param(freqs_nested(j,l) + phi, 1.0 / (sigmal(l) * u(j) + 1.0)));
      } else {
        Tlj(j,l) = randg(distr_param(phi, 1.0));
      }
    }
  }
}

//-----------------------------------------------------
// Update sigmal
//-----------------------------------------------------

void update_sigmal_nested(mat freqs_nested,
                          mat freqs,     
                          vec nested_c,
                          double sig,    
                          double phi,
                          vec u, 
                          mat Tlj,
                          vec &sigmal, 
                          int IS_nval = 100){      
  
  int k = freqs_nested.n_cols;
  // int d = freqs.n_rows; 
  sigmal.resize(k);
  sigmal.fill(0.0);
  
  
  vec temp_location, temp_weight;
  int temp_val;

  for(uword l = 0; l < k; l++){
    temp_location = randg(IS_nval, distr_param(1.0, 0.1));
    temp_weight.resize(IS_nval);

    for(uword i = 0; i < IS_nval; i++){
      temp_weight(i) = (accu(freqs_nested.col(l)) - sig - 1.0) * log(temp_location(i)) - 0.1 * temp_location(i) + log(0.1);
      for(uword j = 0; j < freqs.n_rows; j++){
        temp_weight(i) -= (freqs_nested(j,l) + phi) * log(temp_location(i) * u(nested_c(j)) + 1.0);
      }
    }
    temp_val = rint_not_unif_log(temp_weight);
    sigmal(l) = temp_location(temp_val);
  }
}

//-----------------------------------------------------
// Update mji's
// sample the weights 
//-----------------------------------------------------

void update_mji_nested(vec eta_weights, 
                       vec u, 
                       double phi,
                       vec nested_c, 
                       int d,
                       mat &mji){
  
  int max_idx_eps = eta_weights.n_elem;
  mji.resize(d,max_idx_eps);
  
  for(uword j = 0; j < d; j++){
    for(uword i = 0; i < max_idx_eps; i++){
      if(accu(nested_c == j) > 0){
        mji(j,i) = randg(distr_param(phi, 1.0 / (u(j) * eta_weights(i) + 1.0)));
      } else {
        mji(j,i) = randg(distr_param(phi, 1.0));
      }
    }
  }
}

//-----------------------------------------------------
// function integrand 
// transformed integrand which map in (0,1)
//-----------------------------------------------------

double integrand_levy_inversion(double y,
                                vec u,
                                double J,
                                double sigma,
                                double phi){
  int d = u.n_elem;
  double logint = -(1.0 + sigma) * log(atan(y) * 2.0 / M_PI) + log(2.0 / (M_PI * (1.0 + pow(y, 2.0))));
  // for(uword i = 0; i < d; i++){
  //   logint -= phi * log((u(i) * atan(y) * 2.0 / M_PI) + 1.0);
  // }
  logint -= accu(phi * log((u * atan(y) * 2.0 / M_PI) + 1.0));
  return exp(logint);
}

//-----------------------------------------------------
// function Simpson's rule integration 
// approximate the integral for the tranformed integrand
//-----------------------------------------------------

double integrate_levy_inversion(double temp_S,
                                vec u,
                                double J,
                                double sigma,
                                double phi,
                                int M = 100){
  
  vec xval = linspace(atan(J) * 2.0 / M_PI, 1.0, M + 1);
  double output = 0.0;
  
  for(uword i = 1; i < M/2 + 1; i++){
    output += integrand_levy_inversion(xval(2 * i - 2), u, J, sigma, phi) +
      4.0 * integrand_levy_inversion(xval(2 * i - 1), u, J, sigma, phi) +
      integrand_levy_inversion(xval(2 * i), u, J, sigma, phi);
  }
  output *= (1.0 - atan(J) * 2.0 / M_PI) / (3.0 * M);
  
  // correct the integral
  output *= sigma * exp(lgamma(phi) - lgamma(sigma + phi) - lgamma(1 - sigma));
  output -= temp_S;
  // return the value
  return output;
}

//-----------------------------------------------------
// bisection to find the weight
//-----------------------------------------------------

double bisection_rule(double S, 
                      vec u, 
                      double sigma, 
                      double phi, 
                      int M = 100, 
                      double tol = 0.00001, 
                      int max_iter = 1000){
  int iter = 0;
  double a = 0.0;
  double b = 1.0;
  double m, fm;
  
  while(integrate_levy_inversion(S,u,b,sigma,phi,M) > 0){
    b += 1.0;
  }
  
  while(iter < max_iter){
    
    m = (a + b) / 2.0;
    fm = integrate_levy_inversion(S,u,m,sigma,phi,M);
    
    if(sqrt(pow(fm, 2.0)) < tol){
      break;
    } else {
      if(fm * integrate_levy_inversion(S,u,a,sigma,phi,M) > 0.0){
        a = m;
      } else {
        b = m;
      }
    }
    iter += 1;
  }
  return m;
}

//-----------------------------------------------------
// update CRM eta
//-----------------------------------------------------

void update_eta(vec u, 
                double sigma, 
                double phi, 
                double m0, 
                double k0, 
                double a0, 
                double b0,
                vec &eta_weights, 
                vec &eta_location,
                vec &eta_scale,
                int M = 100, 
                double tol = 0.00001, 
                int max_iter = 1000, 
                double epsilon = 0.001){
  
  eta_weights.resize(0);
  bool cond = true; 
  double S = 0.0; 
  int index = 0;
  
  while(cond){
    S += randg(distr_param(1, 1));
    eta_weights.resize(index + 1);
    eta_weights(index) = bisection_rule(S, u, sigma, phi, M, tol, max_iter);
    
    if(eta_weights(index) < epsilon){
      cond = false;
    }
    index += 1; 
  }
  
  eta_location.resize(eta_weights.n_elem);
  eta_scale.resize(eta_weights.n_elem);
  
  for(uword j = 0; j < eta_weights.n_elem; j++){
    eta_scale(j) = 1 / randg(distr_param(a0, 1.0 / b0)); 
    eta_location(j) = randn() * (eta_scale(j) / k0) + m0;
  }
}

//-----------------------------------------------------
// update U
//-----------------------------------------------------

void update_u(vec &u, 
              int d,
              vec eta_weights, 
              mat mij,
              mat Tlj, 
              vec sigmal, 
              mat freqs, 
              vec nested_c){
  double mass; 
  for(uword j = 0; j < d; j++){
    if(accu(nested_c == j) > 0){
      mass = dot(mij.row(j), eta_weights) + dot(Tlj.row(j), sigmal);
      u(j) = randg(distr_param(accu(freqs.row(j)), 1.0 / mass));
    } else {
      u(j) = 0.0;
    }
  }
}

//-----------------------------------------------------
// update measures allocation
//-----------------------------------------------------

void update_nested_distribution(vec &nested_weigths, 
                                double beta, 
                                vec nested_c){
  
  for(uword i = 0; i < nested_weigths.n_elem; i++){
    nested_weigths(i) = randg(distr_param(beta + sum(nested_c == i), 1.0 / nested_weigths.n_elem));
  }
  nested_weigths /= accu(nested_weigths);
}

void update_nested_allocation(vec nested_weigths, 
                              vec &nested_c, 
                              vec Y, 
                              vec group, 
                              vec eta_weights, 
                              vec eta_location,
                              vec eta_scale,
                              mat mij,
                              mat Tlj, 
                              vec sigmal,
                              vec active_location, 
                              vec active_scale,
                              int d){
  
  vec temp_prob(nested_weigths);
  double accu_mass, temp_val;
  
  for(uword l = 0; l < max(group) + 1; l++){
    temp_prob = nested_weigths;
    
    for(uword i = 0; i < nested_weigths.n_elem; i++){
      for(uword n = 0; n < Y.n_elem; n++){
        
        // for obs n
        accu_mass = 0.0;
        temp_val = 0.0;
        
        if(group(n) == l){
          for(uword j = 0; j < sigmal.n_elem; j++){
            accu_mass += sigmal(j) * Tlj(i, j);
            temp_val += sigmal(j) * Tlj(i, j) * normpdf(Y(n), active_location(j), sqrt(active_scale(j)));
          }
          for(uword j = 0; j < eta_weights.n_elem; j++){
            accu_mass += eta_weights(j) * mij(i, j);
            temp_val += eta_weights(j) * mij(i, j) * normpdf(Y(n), eta_location(j), sqrt(eta_scale(j)));
          }
          temp_prob(i) += log(temp_val) - log(accu_mass);
        }
      }
    }
    nested_c(l) = rint_not_unif_log(temp_prob);
  }
}

//-----------------------------------------------------
// update cluster
//-----------------------------------------------------

void update_cluster_nested(vec Y, 
                           vec &clust, 
                           vec group, 
                           vec nested_c,
                           vec eta_weights, 
                           vec eta_location,
                           vec eta_scale,
                           mat mij,
                           mat Tlj, 
                           vec sigmal,
                           vec active_location, 
                           vec active_scale){
  
  vec prob(eta_weights.n_elem + sigmal.n_elem);
  for(uword i = 0; i < clust.n_elem; i++){
    
    for(uword j = 0; j < sigmal.n_elem; j++){
      prob(j) = sigmal(j) * Tlj(nested_c(group(i)), j) * normpdf(Y(i), active_location(j), sqrt(active_scale(j)));
    }
    
    for(uword j = 0; j < eta_weights.n_elem; j++){
      prob(sigmal.n_elem + j) = eta_weights(j) * mij(nested_c(group(i)), j) * normpdf(Y(i), eta_location(j), sqrt(eta_scale(j)));
    }
    clust(i) = rint_not_unif(prob);
  }
}

//-----------------------------------------------------
// clean parameters
//-----------------------------------------------------

void para_clean(vec &mu,
                vec &s2,
                vec &clust) {
  int k = max(clust) + 1;
  double tmu, ts2;
  for(uword i = 0; i < k; i++){
    if((int) sum(clust == i) == 0){
      for(uword j = k; j > i; j--){
        if((int) sum(clust == j) != 0){
          clust( find(clust == j) ).fill(i);
          break;
        }
      }
    }
  }
  // reduce dimensions
  int u_bound = 0;
  for(uword i = 0; i < k; i++){
    if(accu(clust == i) > 0){
      u_bound += 1;
    }
  }
  // resize object to the correct dimension
  mu.resize(u_bound);
  s2.resize(u_bound);
  mu.fill(0.0);
  s2.fill(0.0);
}

//-----------------------------------------------------
// acceleration step
//-----------------------------------------------------


void accelerate_location_scale(vec Y, 
                               vec clust, 
                               vec &active_location, 
                               vec &active_scale, 
                               double m0, 
                               double k0, 
                               double a0, 
                               double b0){
  double k_n, m_n, a_n, b_n, data_m;
  int nj;
  vec tY;
  
  for(uword j = 0; j < active_location.n_elem; j++){
    nj = sum(clust == j);
    tY = Y.elem(find(clust == j));
    
    data_m = accu(tY) / nj;
    k_n = (k0 + nj);
    m_n = ((m0 * k0) + nj * data_m) / k_n;
    a_n = a0 + (nj / 2.0);
    b_n = b0 + (accu(pow(tY - data_m, 2)) + (nj * k0 * pow(data_m - m0, 2)) / (k_n)) / 2;
    
    active_scale(j) = 1.0 / randg(distr_param(a_n, 1 / b_n));
    active_location(j) = randn() * sqrt(active_scale(j) / k_n) + m_n;
  }
}


//-----------------------------------------------------
// hyper-acceleration step
//-----------------------------------------------------

void hyper_accelerate(vec active_location,
                      vec active_scale,
                      double &m0,
                      double &k0,
                      double a0,
                      double &b0,
                      double m1,
                      double s21,
                      double tau1,
                      double tau2,
                      double a1,
                      double b1){
  
  double m_temp, s2_temp, tau1_temp, tau2_temp, a_temp, b_temp, mu_m;
  int k = active_location.n_elem;
  
  tau1_temp = tau1 + k / 2;
  tau2_temp = tau2 + accu(pow(active_location - m0, 2) / active_scale) / 2;
  k0 = randg(distr_param(tau1_temp, 1 / tau2_temp));
  
  s2_temp = 1 / ( 1 / s21 + k0 * accu(1 / active_scale) );
  m_temp  = s2_temp * ( m1 / s21 + k0 * accu(active_location / active_scale) );
  m0  = randn() * sqrt(s2_temp) + m_temp;
  
  a_temp = a1 + k * a0;
  b_temp = b1 + accu(active_scale);
  b0 = randg(distr_param(a_temp, 1 / b_temp));
}

//-----------------------------------------------------
// evaluate the density for each group
//-----------------------------------------------------

mat eval_dens_comp_nested(int d,
                          vec grid, 
                          vec nested_c,
                          vec active_location,
                          vec active_scale, 
                          vec eta_location, 
                          vec eta_scale,
                          vec eta_weights,
                          mat mij,
                          mat Tlj, 
                          vec sigmal){
  
  mat res(d, grid.n_elem);
  res.fill(0.0);
  double accu_mass;
  
  for(uword r = 0; r < d; r++){
    accu_mass = 0.0;
    
    for(uword j = 0; j < sigmal.n_elem; j++){
      accu_mass += sigmal(j) * Tlj(nested_c(r), j);
      res.row(r) += sigmal(j) * Tlj(nested_c(r), j) * normpdf(grid, active_location(j), sqrt(active_scale(j))).t();
    }
    
    for(uword j = 0; j < eta_weights.n_elem; j++){
      accu_mass += eta_weights(j) * mij(nested_c(r), j);
      res.row(r) += eta_weights(j) * mij(nested_c(r), j) * normpdf(grid, eta_location(j), sqrt(eta_scale(j))).t();
    }
    
    res.row(r) /= accu_mass;
  }
  return res;
}

//-----------------------------------------------------
// update sigma
//-----------------------------------------------------

double integrand_update_sigma(double y, 
                              int val_idx, 
                              mat freqs, 
                              vec u, 
                              double phi,
                              double sigma){
  int d = freqs.n_rows;
  double temp_val = (accu(freqs.col(val_idx)) - sigma - 1.0) * log(atan(y) * 2.0 / M_PI);
  
  for(uword i = 0; i < d; i++){
    temp_val -= (freqs(i, val_idx) + phi) * log((u(i) * 2.0 * atan(y) / M_PI) + 1.0);
  }
  temp_val += log(2.0 / (M_PI * (1.0 + pow(y, 2.0))));
  return(temp_val);
}

double logsumexp(vec v){
  double m_val = max(v);
  double sum = m_val + log(accu(exp(v - m_val))); 
  return sum;
}

double product_integral_update_sigma(mat freqs, 
                                     vec u, 
                                     double phi, 
                                     double sigma, 
                                     int M = 100){
  int k = freqs.n_cols; 
  double temp_int = 0.0;
  double output = 0.0;
  vec temp_val(3 * M / 2);
  
  for(uword l = 0; l < k; l++){
    vec yval = linspace(pow(0.1, 5), 1.0, M + 1);
    for(uword i = 1; i < M/2 + 1; i++){
      temp_val(i - 1) = integrand_update_sigma(yval(2*i - 2), l, freqs, u, phi, sigma);
      temp_val(M / 2 + i - 1) = log(4.0) + integrand_update_sigma(yval(2*i - 1), l, freqs, u, phi, sigma);
      temp_val(M + i - 1) = integrand_update_sigma(yval(2*i), l, freqs, u, phi, sigma);
    }
    temp_int += log(1.0 / (3.0 * M)) + logsumexp(temp_val);
  }
  return temp_int;
}

void update_sigma_MH(double &sigma, 
                     double phi,
                     double MH_var,
                     double a_sigma, 
                     double b_sigma,
                     mat freqs,
                     vec active_location,
                     vec u,  
                     vec sigmal, 
                     int M = 100){
  
  int k = active_location.n_elem;
  double y_old = tan(sigma * M_PI - M_PI_2);
  double y_new = randn() * sqrt(MH_var) + y_old;
  double sigma_new = atan(y_new) / M_PI + 0.5;
  
  double acc_rate = ((pow(sigma_new, a_sigma - 1) * pow(1 - sigma_new, b_sigma - 1)) / 
                     (pow(sigma, a_sigma - 1) * pow(1 - sigma, b_sigma - 1))) * 
                     ((1.0 + pow(y_old, 2.0)) / (1 + pow(y_new, 2.0))) * 
                     exp(product_integral_update_sigma(freqs, u, phi, sigma_new, M)- 
                     product_integral_update_sigma(freqs, u, phi, sigma, M)) * 
                     exp(k * (log(sigma_new) + lgamma(phi + sigma) + lgamma(1 - sigma) -
                     log(sigma) - lgamma(phi + sigma_new) - lgamma(1 - sigma_new)));
  
  if(randu() <= acc_rate){
    sigma = sigma_new;
  }
}

//-----------------------------------------------------
//   MAIN FUNCTION   ----------------------------------
//-----------------------------------------------------

// [[Rcpp::export]]
Rcpp::List main_nested(vec Y, 
                       vec group, 
                       int d,
                       int niter, 
                       int nburn,
                       int thin,
                       double m0,
                       double k0,
                       double a0,
                       double b0,
                       double m1,
                       double s21,
                       double tau1,
                       double tau2,
                       double a1,
                       double b1, 
                       double sigma, 
                       double phi,
                       double a_sigma, 
                       double b_sigma,
                       double beta,
                       vec grid,
                       int M = 100, 
                       double tol = 0.00001, 
                       int max_iter = 1000, 
                       double epsilon = 0.001, 
                       bool eval_density = false,
                       int IS_nval = 100, 
                       int nupd = 0, 
                       double MH_var = 0.1){
  int nitem = (niter - nburn) / thin;
  mat res_clust(nitem, Y.n_elem);
  mat res_nested(nitem, d);
  vec res_sigma(nitem);
  vec res_n_eta(nitem);
  
  // density output
  cube dens_out(1,1,1);
  if(eval_density){
    dens_out.resize(d, grid.n_rows, nitem);
    dens_out.fill(0.0);
  }
  
  // initialize the clusters
  vec clust(Y.n_elem);
  clust.fill(0);
  
  // initialize the vector of parameters
  vec active_location(1);
  vec active_scale(1);
  vec nested_weigths(d);
  vec nested_c = regspace(0, d - 1);
  
  active_location.fill(1.0);
  active_scale.fill(1.0);
  nested_weigths.fill(1.0 / d);
  
  // initialize the quantities for the CRM
  vec eta_weights;
  vec eta_location;
  vec eta_scale;
  
  // initialize freqs
  mat freqs, freqs_nested;
  update_freqs(d, group, clust, freqs);
  update_freqs_nested(d, group, nested_c, clust, freqs_nested);
  // initialize the augmenting quantity U
  vec u = randg(d);
  
  // initialize the compound weights
  mat mij;
  mat Tlj(d,1);
  Tlj.fill(1.0);
  
  vec sigmal(1);
  sigmal.fill(1.0);
  
  // initialize the times
  int start_s = clock();
  int current_s;
  if(nupd == 0){
    nupd = niter / 10;
  }
  
  // main loop
  int res_index = 0;
  Rcpp::Rcout << "\nStarting the algorithm...\n" ;
  for(uword iter = 0; iter < niter; iter ++){
    
    update_sigmal_nested(freqs_nested, freqs, nested_c, sigma, phi, u, Tlj, sigmal);
    update_Tlj_nested(freqs_nested, d, phi, nested_c, u, sigmal, Tlj);
    update_eta(u, sigma, phi, m0, k0, a0, b0,
               eta_weights, eta_location, eta_scale,
               M, tol, max_iter, epsilon);
    update_mji_nested(eta_weights, u, phi, nested_c, d, mij);
    update_nested_distribution(nested_weigths, beta, nested_c);
    update_nested_allocation(nested_weigths, nested_c, Y, group, eta_weights, eta_location,
                             eta_scale, mij, Tlj, sigmal, active_location, active_scale, d);
    update_cluster_nested(Y, clust, group, nested_c, eta_weights, eta_location, eta_scale, mij, Tlj, 
                   sigmal, active_location, active_scale);
    
    // save the density values
    if((iter >= nburn) & ((iter + 1) % thin == 0) & eval_density){
      dens_out.slice(res_index) = eval_dens_comp_nested(d, grid, nested_c, active_location, active_scale,  
                     eta_location, eta_scale, eta_weights, mij, Tlj, sigmal);
    }
    
    para_clean(active_location, active_scale, clust);
    update_freqs(d, group, clust, freqs);
    update_freqs_nested(d, group, nested_c, clust, freqs_nested);
    accelerate_location_scale(Y, clust, active_location, active_scale, 
                              m0, k0, a0, b0);
    hyper_accelerate(active_location, active_scale, m0, k0, a0, b0,
                     m1, s21, tau1, tau2, a1, b1);
    // update_sigma_MH(sigma, phi, MH_var, a_sigma, b_sigma, freqs_nested, active_location, u, sigmal, M);
    update_u(u, d, eta_weights, mij, Tlj, sigmal, freqs_nested, nested_c);
    
    // save the results
    if((iter >= nburn) & ((iter + 1) % thin == 0)){
      res_clust.row(res_index) = clust.t();
      res_nested.row(res_index) = nested_c.t();
      res_sigma(res_index) = sigma;
      res_n_eta(res_index) = eta_weights.n_elem;
      res_index += 1;
    }
    
    // print the time and stop check  
    if((iter + 1) % nupd == 0){
      current_s = clock();
      Rcpp::Rcout << "Completed:\t" << (iter + 1) << "/" << niter << " - in " <<
        double(current_s-start_s)/CLOCKS_PER_SEC << " sec\n";
    }
    Rcpp::checkUserInterrupt();
  }
  
  Rcpp::List res;
  res["time"]   = double(clock() - start_s)/CLOCKS_PER_SEC;
  res["clust"]  = res_clust;
  res["nested"] = res_nested;
  res["dens"]   = dens_out;
  res["sigma"]  = res_sigma;
  res["n_eta"]  = res_n_eta;
  return res;
}

