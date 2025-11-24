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
  vec probs;
  double mw = max(lweights);
  
  probs = exp(lweights - mw);
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
// freqs: matrix, each row a clust, each column an element
// Tlj: matrix, each row a clust, each column an element

//-----------------------------------------------------
// Update freqs
// compute the frequencies for each group and each global active cluster

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

//-----------------------------------------------------
// Update Tlj's
// sample the group-and-cluster specific weights for the existent atoms
// the full conditional is a Gamma(n_lj + phi, sigma_l u_j + 1)

void update_Tlj(mat freqs,      // dim freq(d,k)
                int d,
                double phi, 
                vec u, 
                vec sigmal,     // vector of size k
                mat &Tlj){      // dim Tlj(d,k)
  
  int k = freqs.n_cols;
  Tlj.resize(d, k);
  Tlj.fill(0.0);
  
  // Rcpp::Rcout << "\n\n-----------\n\n" << Tlj << "\n\n" << freqs << "\n\n" << sigmal << "\n\n" << u << "\n\n";
  for(uword l = 0; l < k; l++){
    for(uword j = 0; j < d; j++){
      Tlj(j,l) = randg(distr_param(freqs(j,l) + phi, 1.0 / (sigmal(l) * u(j) + 1.0)));
    }
  }
}

//-----------------------------------------------------
// Update sigmal
// sample the atom specific weight, where 
// the full conditional is a Gamma(n_l. - sig, sum(u_j * T_lj))

void update_sigmal(mat freqs,     // dim freq(d,k)
                   double sig,    // overall discount
                   double phi,
                   vec u, 
                   vec &sigmal,   // vector of size k 
                   int IS_nval = 100){      
  
  int k = freqs.n_cols;
  int d = freqs.n_rows; 
  sigmal.resize(k);
  sigmal.fill(0.0);
  
  
  vec temp_location, temp_weight;
  int temp_val; 
  
  // Rcpp::Rcout << "\n\n--------------------\n\n" <<
  //   freqs << "\n\n" <<
  //     u << "\n\n" <<
  //       sigmal << "\n\n" << k << "\n\n";
  for(uword l = 0; l < k; l++){
    temp_location = randg(IS_nval, distr_param(1.0, 0.1));
    temp_weight.resize(IS_nval);
    
    for(uword i = 0; i < IS_nval; i++){
      temp_weight(i) = (accu(freqs.col(l)) - sig - 1.0) * log(temp_location(i)) - 
        dot((freqs.col(l) + phi), log(temp_location(i) * u + 1.0)) - 0.1 * temp_location(i) + log(0.1);  
    }
    // Rcpp::Rcout << "here" << "\n\n";
    temp_val = rint_not_unif_log(temp_weight);
    // Rcpp::Rcout << "here" << "\n\n" << temp_val << "\n\n";
    sigmal(l) = temp_location(temp_val);
  }
}

//-----------------------------------------------------
// Update mji's
// sample the weights 

void update_mji(vec eta_weights, 
                vec u, 
                double phi,
                int d,
                mat &mji){
  
  int max_idx_eps = eta_weights.n_elem;
  mji.resize(d,max_idx_eps);
  // Rcpp::Rcout << "\n\n" << d << "\t" << max_idx_eps;
  // Rcpp::Rcout << "\n\n" << u.n_elem << "\n\n";
  
  for(uword j = 0; j < d; j++){
    for(uword i = 0; i < max_idx_eps; i++){
      mji(j,i) = randg(distr_param(phi, 1 / (u(j) * eta_weights(i) + 1)));
    }
  }
}

//-----------------------------------------------------
// function integrand 
// transformed integrand which map in (0,1)

// // [[Rcpp::export]]
// double integrand_levy_inversion(double S,
//                                 vec u,
//                                 double J,
//                                 double sigma,
//                                 double phi){
//   int d = u.n_elem;
//   double logint = (phi * d + sigma - 1) * log(S);
//   for(uword i = 0; i < d; i++){
//     logint -= phi * log(S + exp(log(u(i)) + log(J)));
//   }
//   return exp(logint);
// }

double integrand_levy_inversion(double y,
                                vec u,
                                double J,
                                double sigma,
                                double phi){
  int d = u.n_elem;
  double logint = -(1.0 + sigma) * log(atan(y) * 2.0 / M_PI) + log(2.0 / (M_PI * (1 + pow(y, 2.0))));
  for(uword i = 0; i < d; i++){
    logint -= phi * log((u(i) * atan(y) * 2.0 / M_PI) + 1.0);
  }
  return exp(logint);
}

//-----------------------------------------------------
// function Simpson's rule integration 
// approximate the integral for the tranformed integrand

// // [[Rcpp::export]]
// double integrate_levy_inversion(double temp_S,
//                                 vec u,
//                                 double J,
//                                 double sigma,
//                                 double phi,
//                                 int M = 100){
// 
//   vec xval = linspace(pow(0.1, 25), 1.0, M + 1);
//   double output = 0.0;
// 
//   for(uword i = 1; i < M/2 + 1; i++){
//     output += integrand_levy_inversion(xval(2*i - 2), u, J, sigma, phi) +
//       4.0 * integrand_levy_inversion(xval(2*i - 1), u, J, sigma, phi) +
//       integrand_levy_inversion(xval(2*i), u, J, sigma, phi);
//   }
//   output *= 1.0 / (3.0 * M);
// 
//   // correct the integral
//   output *= sigma * exp(lgamma(phi) - lgamma(sigma + phi) - lgamma(1 - sigma));
//   output -= temp_S * pow(J, sigma);
//   // return the value
//   return output;
// }

double integrate_levy_inversion(double temp_S,
                                vec u,
                                double J,
                                double sigma,
                                double phi,
                                int M = 100){

  vec xval = linspace(atan(J) * 2.0 / M_PI, 1.0, M + 1);
  double output = 0.0;

  for(uword i = 1; i < M/2 + 1; i++){
    output += integrand_levy_inversion(xval(2*i - 2), u, J, sigma, phi) +
      4.0 * integrand_levy_inversion(xval(2*i - 1), u, J, sigma, phi) +
      integrand_levy_inversion(xval(2*i), u, J, sigma, phi);
  }
  output *= (1.0 - atan(J) * 2.0 / M_PI) / (3.0 * M);

  // correct the integral
  output *= sigma * exp(lgamma(phi) - lgamma(sigma + phi) - lgamma(1 - sigma));
  output -= temp_S;
  // return the value
  return output;

  // Rcpp::Rcout << J << "\t" << output << "\t" << sigma << "\t" << temp_S << "\n\n";
}

//-----------------------------------------------------
// bisection to find the weight

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
  
  // Rcpp::Rcout << "\n----\n----\n----\n" << sigma << "\t" << phi << "\n" << u.t() << "\n"
  //   << "\t" << S << "\n----\n----\n----\n";
  
  // Rcpp::Rcout << "\nQUI_eta"; 
  while(cond){
    // Rcpp::Rcout << "\t CONT_TRUE"; 
    S += randg(distr_param(1, 1));
    eta_weights.resize(index + 1);
    eta_weights(index) = bisection_rule(S, u, sigma, phi, M, tol, max_iter);
    
    // Rcpp::Rcout << "\n----\n" << eta_weights(index) << "\t" << S << "\n----\n";
    if(eta_weights(index) < epsilon){
      cond = false;
    }
    index += 1; 
  }
  // Rcpp::Rcout << "\n\n" << eta_weights.n_elem;
  eta_location.resize(eta_weights.n_elem);
  eta_scale.resize(eta_weights.n_elem);
  
  for(uword j = 0; j < eta_weights.n_elem; j++){
    eta_scale(j) = 1 / randg(distr_param(a0, 1 / b0)); 
    eta_location(j) = randn() * (eta_scale(j) / k0) + m0;
  }
}

//-----------------------------------------------------
// update U

void update_u(vec &u, 
              vec eta_weights, 
              mat mij,
              mat Tlj, 
              vec sigmal, 
              mat freqs){
  int d = u.n_elem;
  double mass; 
  for(uword j = 0; j < d; j++){
    mass = dot(Tlj.row(j), sigmal) + dot(mij.row(j), eta_weights);
    // mass = dot(mij.row(j), eta_weights);
    u(j) = randg(distr_param(accu(freqs.row(j)), 1 / mass));
  }
}

//-----------------------------------------------------
// update cluster

void update_cluster(vec Y, 
                    vec &clust, 
                    vec group, 
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
      prob(j) = sigmal(j) * Tlj(group(i), j) * normpdf(Y(i), active_location(j), sqrt(active_scale(j)));
    }
    
    for(uword j = 0; j < eta_weights.n_elem; j++){
      prob(sigmal.n_elem + j) = eta_weights(j) * mij(group(i), j) * normpdf(Y(i), eta_location(j), sqrt(eta_scale(j)));
    }
    
    // Rcpp::Rcout << prob.t() << "\n\n";
    clust(i) = rint_not_unif(prob);
  }
}

//-----------------------------------------------------
// clean parameters

void para_clean_ICS(vec &mu,
                    vec &s2,
                    vec &clust) {
  int k = mu.n_elem;
  double tmu, ts2;
  
  // for all the used parameters
  for(uword i = 0; i < k; i++){
    
    // if a cluster is empty
    if((int) sum(clust == i) == 0){
      
      // find the last full cluster, then swap
      for(uword j = k; j > i; j--){
        if((int) sum(clust == j) != 0){
          
          // swap the correpsonding elements
          clust( find(clust == j) ).fill(i);
          
          tmu = mu[i];
          mu[i] = mu[j];
          mu[j] = tmu;
          
          ts2 = s2[i];
          s2[i] = s2[j];
          s2[j] = ts2;
          
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
}

//-----------------------------------------------------
// acceleration step

void accelerate_location_scale(vec Y, 
                               vec clust, 
                               vec group, 
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

mat eval_dens_comp(int d,
                   vec grid, 
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
      accu_mass += sigmal(j) * Tlj(r, j);
      res.row(r) += sigmal(j) * Tlj(r, j) * normpdf(grid, active_location(j), sqrt(active_scale(j))).t();
    }
    
    for(uword j = 0; j < eta_weights.n_elem; j++){
      accu_mass += eta_weights(j) * mij(r, j);
      res.row(r) += eta_weights(j) * mij(r, j) * normpdf(grid, eta_location(j), sqrt(eta_scale(j))).t();
    }
    
    res.row(r) /= accu_mass;
  }
  return res;
}

//-----------------------------------------------------
// evaluate the density ancestor

vec eval_ancestor(int d,
                  vec grid, 
                  vec active_location,
                  vec active_scale, 
                  vec eta_location, 
                  vec eta_scale,
                  vec eta_weights,
                  vec sigmal){
  
  vec res(grid.n_elem);
  res.fill(0.0);
  double accu_mass;
  
  accu_mass = 0.0;
  
  for(uword j = 0; j < sigmal.n_elem; j++){
    accu_mass += sigmal(j);
    res += sigmal(j) * normpdf(grid, active_location(j), sqrt(active_scale(j)));
  }
  
  for(uword j = 0; j < eta_weights.n_elem; j++){
    accu_mass += eta_weights(j);
    res += eta_weights(j) * normpdf(grid, eta_location(j), sqrt(eta_scale(j)));
  }
  
  res /= accu_mass;
  
  return res;
}
//-----------------------------------------------------
// update sigma

double integrand_update_sigma(double y, 
                              int val_idx, 
                              mat freqs, 
                              vec u, 
                              double phi,
                              double sigma){
  int d = freqs.n_rows;
  double temp_val = (accu(freqs.col(val_idx)) - sigma - 1.0) * log(atan(y) * 2.0 / M_PI);
  
  // Rcpp::Rcout << "\n-----------\n" << temp_val << "\t";
  for(uword i = 0; i < d; i++){
    temp_val -= (freqs(i, val_idx) + phi) * log((u(i) * 2.0 * atan(y) / M_PI) + 1.0);
    // Rcpp::Rcout << temp_val << "\t";
  }
  temp_val += log(2.0 / (M_PI * (1.0 + pow(y, 2.0))));
  // Rcpp::Rcout << temp_val << "\t";
  // return exp(temp_val);
  return(temp_val);
}

double logsumexp(vec v){
  double m_val = max(v);
  double sum = m_val + log(accu(exp(v - m_val))); 
  return sum;
}

// double product_integral_update_sigma(mat freqs, 
//                                      vec u, 
//                                      double phi, 
//                                      double sigma, 
//                                      int M = 100){
//   int k = freqs.n_cols; 
//   double temp_int = 1.0;
//   double output = 0.0;
//   
//   for(uword l = 0; l < k; l++){
//     vec yval = linspace(pow(0.1, 5), 1.0, M + 1);
//     output = 0.0;
//     for(uword i = 1; i < M/2 + 1; i++){
//       output += integrand_update_sigma(yval(2*i - 2), l, freqs, u, phi, sigma) + 
//         4.0 * integrand_update_sigma(yval(2*i - 1), l, freqs, u, phi, sigma) + 
//         integrand_update_sigma(yval(2*i), l, freqs, u, phi, sigma);
//     }
//     temp_int *= 1.0 / (3.0 * M) * output;
//   }
//   return temp_int;
// }

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
  
  // int k = active_location.n_elem;
  // double y_old = tan(sigma * M_PI - M_PI_2);
  // double y_new = randn() * sqrt(MH_var) + y_old;
  // double sigma_new = atan(y_new) / M_PI + 0.5;
  
  // // double E = mean(pow(sigmal, sigma - sigma_p));
  // double acc_rate = (pow(sigma_new, a_sigma - 1) * pow(1 - sigma_new, b_sigma - 1)) / 
  //                    (pow(sigma, a_sigma - 1) * pow(1 - sigma, b_sigma - 1));
  //                    
  // Rcpp::Rcout << "\n--------\n" << acc_rate << "\n\n";
  // acc_rate *= 
  //   ((1.0 + pow(y_old, 2.0)) / (1 + pow(y_new, 2.0))); 
  // 
  // 
  // Rcpp::Rcout << "\n\n" << acc_rate << "\n"
  //             << product_integral_update_sigma(freqs, u, phi, sigma_new, M) << "\n"
  //             << product_integral_update_sigma(freqs, u, phi, sigma, M) << "\n";
  // acc_rate *= 
  //   product_integral_update_sigma(freqs, u, phi, sigma_new, M) / 
  //   product_integral_update_sigma(freqs, u, phi, sigma, M); 
  // 
  // 
  // Rcpp::Rcout << "\n\n" << acc_rate << "\n\n";
  // acc_rate *= 
  //   exp(k * (log(sigma_new) + lgamma(phi + sigma) + lgamma(1 - sigma) -
  //       log(sigma) - lgamma(phi + sigma_new) - lgamma(1 - sigma_new)));
  // 
  // if(randu() <= acc_rate){
  //   sigma = sigma_new;
  // }
  
  int k = active_location.n_elem;
  double y_old = tan(sigma * M_PI - M_PI_2);
  double y_new = randn() * sqrt(MH_var) + y_old;
  double sigma_new = atan(y_new) / M_PI + 0.5;
  
  // double E = mean(pow(sigmal, sigma - sigma_p));
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
  
  // Rcpp::Rcout << "\n\n" << acc_rate << "\n\n";
}

// void update_phi_MH(double sigma, 
//                    double &phi,
//                    double MH_var,
//                    mat freqs,
//                    vec active_location,
//                    vec u,  
//                    vec sigmal, 
//                    int M = 100){
//   
//   int k = active_location.n_elem;
//   double y_old = log(phi);
//   double y_new = randn() * sqrt(MH_var) + y_old;
//   double phi_new = exp(y_new);
//   
//   // double E = mean(pow(sigmal, sigma - sigma_p));
//   double acc_rate = exp(y_new - y_old) * 
//     product_integral_update_sigma(freqs, u, phi_new, sigma, M) / 
//     product_integral_update_sigma(freqs, u, phi, sigma, M) * 
//     exp(k * (lgamma(phi_new) + lgamma(phi + sigma) -
//     lgamma(phi) - lgamma(phi_new + sigma)));
//   
//   if(randu() <= acc_rate){
//     phi = phi_new;
//   }
//   // Rcpp::Rcout << "\n\n" << phi << "\n\n";
// }

void update_phi_MH(double sigma, 
                   double &phi,
                   double MH_var,
                   double a_phi, 
                   double b_phi,
                   mat freqs,
                   vec active_location,
                   vec u,  
                   vec sigmal, 
                   int M = 100){
  
  int k = active_location.n_elem;
  double y_old = log(phi);
  double y_new = randn() * sqrt(MH_var) + y_old;
  double phi_new = exp(y_new);
  
  // double E = mean(pow(sigmal, sigma - sigma_p));
  // double acc_rate = exp(y_new - y_old) * 
  //   product_integral_update_sigma(freqs, u, phi_new, sigma, M) / 
  //   product_integral_update_sigma(freqs, u, phi, sigma, M) * 
  //   exp(k * (lgamma(phi_new) + lgamma(phi + sigma) -
  //   lgamma(phi) - lgamma(phi_new + sigma)));
  
  double acc_rate = pow(phi_new / phi, a_phi) * exp(- b_phi * (phi_new - phi)) *  
    exp(product_integral_update_sigma(freqs, u, phi_new, sigma, M) -
    product_integral_update_sigma(freqs, u, phi, sigma, M)) *
    exp(k * (lgamma(phi_new) + lgamma(phi + sigma) -
    lgamma(phi) - lgamma(phi_new + sigma)));
  
  for(uword i = 0; i < freqs.n_cols; i++){
    for(uword j = 0; j < freqs.n_rows; j++){
      if(freqs(j,i) > 0){
        acc_rate *= exp(lgamma(phi_new + freqs(j,i)) - lgamma(phi_new) -
          lgamma(phi + freqs(j,i)) + lgamma(phi));
      }
    }
  }
  
  if(randu() <= acc_rate){
    phi = phi_new;
  }
  // Rcpp::Rcout << "\n" << phi << "\n";
}

//-----------------------------------------------------
//-----------------------------------------------------
//-----------------------------------------------------

// [[Rcpp::export]]
Rcpp::List main_compound(vec Y, 
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
                         double a_phi, 
                         double b_phi,
                         vec grid,
                         int M = 100, 
                         double tol = 0.00001, 
                         int max_iter = 1000, 
                         double epsilon = 0.001, 
                         bool eval_density = false,
                         int IS_nval = 100, 
                         int nupd = 0, 
                         double MH_var = 0.1, 
                         double MH_var2 = 0.1){
  int nitem = (niter - nburn) / thin;
  mat res_clust(nitem, Y.n_elem);
  vec res_sigma(nitem);
  vec res_phi(nitem);
  vec res_n_eta(nitem);
  
  // density output
  cube dens_out(1,1,1);
  mat ancestor_out(1,1);
  if(eval_density){
    dens_out.resize(d, grid.n_rows, nitem);
    dens_out.fill(0.0);
    
    ancestor_out.resize(nitem, grid.n_rows);
    ancestor_out.fill(0.0);
  }
  
  // initialize the clusters
  vec clust(Y.n_elem);
  clust.fill(0);
  
  // initialize the vector of parameters
  vec active_location(1);
  vec active_scale(1);
  active_location.fill(1.0);
  active_scale.fill(1.0);
  
  // initialize the quantities for the CRM
  vec eta_weights;
  vec eta_location;
  vec eta_scale;
  
  // initialize freqs
  mat freqs;
  update_freqs(d, group, clust, freqs);
  
  // initialize the augmenting quantity U
  vec u = randg(d);
  
  // initialize the compound weights
  mat mij;
  mat Tlj;
  vec sigmal(1);
  sigmal.fill(1.0);
  
  // initialize the times
  int start_s = clock();
  int current_s;
  if(nupd == 0){
    nupd = niter / 10;
  }
  
  // Rcpp::Rcout << "\nqui0";
  // main loop
  int res_index = 0;
  Rcpp::Rcout << "\nStarting the algorithm...\n" ;
  for(uword iter = 0; iter < niter; iter ++){
    // Rcpp::Rcout << "\nqui1";
    update_sigmal(freqs, sigma, phi, u, sigmal, IS_nval);
    // Rcpp::Rcout << "\nqui2";
    // Rcpp::Rcout << "\n\n " << freqs << "\n\n " <<  u << "\n\n " << sigmal << "\n\n" << Tlj << "\n\n";
    update_Tlj(freqs, d, phi, u, sigmal, Tlj);
    // Rcpp::Rcout << "\nqui3\t\t";
    update_eta(u, sigma, phi, m0, k0, a0, b0,
               eta_weights, eta_location,eta_scale,
               M, tol, max_iter, epsilon);
    // Rcpp::Rcout << "\nqui4";
    update_mji(eta_weights, u, phi, d, mij);
    // Rcpp::Rcout << "\nqui5";
    update_cluster(Y, clust, group, eta_weights, eta_location, eta_scale, mij, Tlj, 
                   sigmal, active_location, active_scale);
    
    if((iter >= nburn) & ((iter + 1) % thin == 0) & eval_density){
      dens_out.slice(res_index) = eval_dens_comp(d, grid, active_location, active_scale,  
                     eta_location, eta_scale, eta_weights, mij, Tlj, sigmal);
      
      ancestor_out.row(res_index) = eval_ancestor(d, grid, active_location,
                   active_scale, eta_location, eta_scale, eta_weights, sigmal).t();
    }
    
    active_location.resize(max(clust) + 1);
    active_scale.resize(max(clust) + 1);
    active_location.fill(0.0);
    active_scale.fill(0.0);
    // Rcpp::Rcout << "\nqui6";
    para_clean_ICS(active_location, active_scale, clust);
    // Rcpp::Rcout << "\nqui7";
    update_freqs(d, group, clust, freqs);
    // Rcpp::Rcout << "\nqui8";
    accelerate_location_scale(Y, clust, group, active_location, active_scale, 
                              m0, k0, a0, b0);
    // Rcpp::Rcout << "\nqui9";
    hyper_accelerate(active_location, active_scale, m0, k0, a0, b0,
                     m1, s21, tau1, tau2, a1, b1);
    // Rcpp::Rcout << "\nqui10";
    update_sigma_MH(sigma, phi, MH_var, a_sigma, b_sigma, freqs, active_location, u, sigmal, M);
    // Rcpp::Rcout << "\nqui11";
    update_u(u, eta_weights, mij, Tlj, sigmal, freqs);
    // 
    update_phi_MH(sigma, phi, MH_var2, a_phi, b_phi, freqs, active_location, u, sigmal, M);
    // Rcpp::Rcout << "\n" << sigma;
    // save the results
    if((iter >= nburn) & ((iter + 1) % thin == 0)){
      res_clust.row(res_index) = clust.t();
      res_sigma(res_index) = sigma;
      res_phi(res_index) = phi;
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
  res["dens"]   = dens_out;
  res["sigma"]  = res_sigma;
  res["phi"]    = res_phi;
  res["n_eta"]  = res_n_eta;
  res["ancestor"] = ancestor_out;
  return res;
}

