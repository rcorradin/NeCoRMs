// This code is mainly based on the implementation
// of Dr. Francesco Denti, which can be found at
//        https://github.com/Fradenti/CommonAtomModel/blob/master/Functions/CAM

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

// ----------------------------
// Update Uj
void Update_Uj(vec &Uj, 
               vec &xi_z,
               int &NN_z,
               vec zj, 
               double kappa){
  vec tvec(2);
  for(uword j = 0; j < Uj.n_elem; j++){
    Uj(j) = randu() * (1 - kappa) * pow(kappa, (zj(j) + 1));
  }
  tvec(0) = max(1 + floor((log(Uj) - log(1 - kappa)) / log(kappa)));
  tvec(1) = max(zj);
  NN_z = max(tvec);
  xi_z.resize(NN_z);
  for(uword j = 0; j < NN_z; j++){
    xi_z(j) = (1 - kappa) * pow(kappa, (j + 1));
  }
}

// ----------------------------
// Update Uij

void Update_Uij(vec &Uij, 
                vec &xi_c,
                int &NN_c,
                vec cij, 
                double kappa){
  vec tvec(2);
  for(uword j = 0; j < Uij.n_elem; j++){
    Uij(j) = randu() *  (1 - kappa) * pow(kappa, (cij(j) + 1));
  }
  tvec(0) = max(1 + floor((log(Uij) - log(1 - kappa)) / log(kappa)));
  tvec(1) = max(cij);
  NN_c = max(tvec);
  xi_c.resize(NN_c);
  for(uword j = 0; j < NN_c; j++){
    xi_c(j) = (1 - kappa) * pow(kappa, (j + 1));
  }
}

// ----------------------------
// update distributional sticks

void Update_Distributional_Sticks(vec &pi, 
                                  vec &pi_v, 
                                  vec zj, 
                                  int NN_z, 
                                  double alpha){
  pi.resize(NN_z);
  pi_v.resize(NN_z);
  double tempX, tempY;
  for(uword j = 0; j < NN_z; j++){
    tempX = randg(distr_param(1 + accu(find(zj == j)), 1.0));
    tempY = randg(distr_param(alpha + accu(find(zj > j)), 1.0));
    pi_v(j) = tempX / (tempX + tempY);
    if(j == 0){
      pi(j) = pi_v(j);
    } else {
      pi(j) = pi_v(j) * ((1 - pi_v(j-1)) * pi(j - 1)) / pi_v(j-1);
    }
  }
}

// ----------------------------
// update obs sticks

void Update_Observational_Sticks(mat &omega, 
                                 mat &v_omega, 
                                 vec cij, 
                                 vec zj_pg_rep,
                                 int NN_c, 
                                 int NN_z, 
                                 double beta){
  
  
  omega.resize(NN_c, NN_z);
  v_omega.resize(NN_c, NN_z);
  omega.fill(0);
  v_omega.fill(0);
  double tempX, tempY;
  
  for(uword j = 0; j < NN_z; j++){
    for(uword i = 0; i < NN_c; i++){
      tempX = randg(distr_param(1 + accu(cij.elem(find(zj_pg_rep == j)) == i), 1.0));
      tempY = randg(distr_param(beta + accu(cij.elem(find(zj_pg_rep == j)) >  i), 1.0));
      v_omega(i,j) = tempX / (tempX + tempY);
      if(i == 0){
        omega(i,j) = v_omega(i,j);
      } else {
        omega(i,j) = v_omega(i,j) * ((1 - v_omega(i-1,j)) * omega(i-1,j)) / v_omega(i-1,j);
      }
    }
  }
}


// ----------------------------
// update Zj

void Update_Zj(vec &zj_pg, 
               vec Uj,    
               vec xi_z, 
               vec xi_c,
               vec pi_z,
               vec cij,
               mat omega,
               vec y_group,
               int NN_z, 
               int J){
  
  vec p_zj_k(NN_z), subCij; 
  mat subOmega; 
  uvec ind, indCij;
  
  for(uword q = 0; q < J; q++){
    ind = find(y_group == q);
    subCij = cij.elem(ind);
    indCij = arma::conv_to<arma::uvec>::from(subCij);
    subOmega = omega.rows(indCij);
    for(uword k = 0; k < NN_z; k++){
      p_zj_k(k) = log(Uj(q) < xi_z(k)) - log(xi_z(k)) + 
        log(pi_z(k)) + accu(log(subOmega.col(k)));
    }
    if(is_finite(max(p_zj_k))){
      zj_pg(q) = rint_not_unif_log(p_zj_k);  
    }else{
      p_zj_k.fill(0.0);
      zj_pg(q) = rint_not_unif_log(p_zj_k);  
    }     
  }
}

// ----------------------------
// update zj_pg_rep

void Update_zj_pg_rep(vec &zj_pg_rep,
                      vec zj_pg, 
                      vec group){
  for(uword i = 0; i < zj_pg_rep.n_elem; i++){
    zj_pg_rep(i) = zj_pg(group(i));
  }
}

// ----------------------------
// update Cij

void Update_Cij(vec &cij,
                vec y_obser,
                vec Uij,
                vec xi_c,
                mat omega,
                vec zj_pg_rep,
                mat theta,
                int N, 
                int NN_c){
  vec p(NN_c);
  
  for(uword i = 0; i < N; i++){
    for(uword k = 0; k < NN_c; k++){
      p(k) = log(xi_c(k) > Uij(i)) + log(omega(k, zj_pg_rep(i))) - 
        log(xi_c(k)) + log_normpdf(y_obser(i), theta(k,0), sqrt(theta(k,1)));    
    }
    
    if(is_finite(max(p))){
      cij(i) = rint_not_unif_log(p); 
    }else{
      p.fill(0.0);
      cij(i) = rint_not_unif_log(p); 
    }
  }
}

// ----------------------------
// update theta

void Update_theta(mat &theta,
                  vec Y,
                  vec cij,
                  double a0, 
                  double b0, 
                  double k0, 
                  double m0,
                  int NN_c, 
                  int J){
  
  theta.resize(NN_c,2);
  double astar, bstar, mustar, kstar, data_m;
  vec tY;
  
  for(uword i = 0; i < NN_c; i++){
    int nj = accu(cij == i);
    tY = Y.elem(find(cij == i));
    
    if(nj > 0){
      data_m = accu(tY) / nj;
      kstar = (k0 + nj);
      mustar = ((m0 * k0) + nj * data_m) / kstar;
      astar = a0 + (nj / 2.0);
      bstar = b0 + (accu(pow(tY - data_m, 2)) + (nj * k0 * pow(data_m - m0, 2)) / (kstar)) / 2.0;
    } else {
      kstar = k0;
      mustar = m0;
      astar = a0;
      bstar = b0;
    }
    
    theta(i,1) = 1.0 / randg(distr_param(astar, 1.0 / bstar));
    theta(i,0) = randn() * sqrt(theta(i,1) / kstar) + mustar;  
  }
}

// ----------------------------
// evaluate the density on a grid

mat eval_dens(vec grid, 
              mat theta,
              mat omega,
              vec Zj,
              int NN_c, 
              int J){
  mat dens_out(J, grid.n_rows);
  dens_out.fill(0.0);
  double taccu;
  
  for(uword j = 0; j < J; j++){
    taccu = accu(omega.col(Zj(j)));
    for(uword i = 0; i < NN_c; i++){
      dens_out.row(j) += omega(i,Zj(j)) / taccu * normpdf(grid, theta(i,0), sqrt(theta(i,1))).t();
    }
  }
  return(dens_out);
}

// --------------------
// update hyperparam alpha

void update_alpha(double &alpha, 
                  vec zj, 
                  double a_alpha, 
                  double b_alpha){
  vec uzj = unique(zj);
  int k = uzj.n_elem; 
  
  double tempX = randg(distr_param(alpha + 1, 1.0));
  double tempY = randg(distr_param(zj.n_elem, 1.0));
  double logeta = log(tempX / (tempX + tempY));
  
  double Q = (a_alpha + k - 1) / (zj.n_elem * (b_alpha - logeta));
  double pi_eta = Q / (1 + Q);
  double u = randu();
  if (u < pi_eta) {
    alpha = randg(distr_param(a_alpha + k, 1.0 / (b_alpha - logeta)));
  } else{
    alpha = randg(distr_param(a_alpha + k - 1, 1.0 / (b_alpha - logeta)));
  }
}

// --------------------
// update hyperparam beta

void update_beta(double &beta,
                 vec zj,
                 vec zj_pg_rep, 
                 vec cij, 
                 mat omega_v,
                 double gamma, 
                 double delta){
  double Mbar; 
  double astar = gamma; 
  double bstar = delta; 
  double Sbar = max(zj);
  
  for(uword j = 0; j < zj.n_elem; j++){
    if(accu(zj_pg_rep == j) > 0){
      Mbar = max(cij.elem(find(zj_pg_rep == j)));
      astar += Mbar; 
      for(uword i = 0; i < Mbar; i++){
        bstar -= log(1 - omega_v(i,j));
      }
    } 
  }
  beta = randg(distr_param(astar, 1.0 / bstar));
}

// --------------------
// MAIN FUNCTION

// [[Rcpp::export]]
Rcpp::List main_CAM(vec Y, 
                    vec group, 
                    int J,
                    int niter, 
                    int nburn,
                    int thin,
                    double m0,
                    double k0,
                    double a0,
                    double b0,
                    double p_alpha_1,
                    double p_alpha_2,
                    double p_beta_1, 
                    double p_beta_2,
                    double kappa,
                    vec grid,
                    bool eval_density = false,
                    int nupd = 0){
  
  int N = Y.n_elem;
  int nitem = (niter - nburn) / thin;
  mat res_clust(nitem, Y.n_elem);
  mat res_nested(nitem, J);
  
  double alpha = randg(distr_param(p_alpha_1, 1.0 / p_alpha_1));
  double beta = randg(distr_param(p_beta_1, 1.0 / p_beta_2));
  
  // density output
  cube dens_out(1,1,1);
  if(eval_density){
    dens_out.resize(J, grid.n_rows, nitem);
    dens_out.fill(0.0);
  }
  
  // initialize CAM quantities 
  vec pi, pi_v, cij(N), zj_pg, zj_pg_rep(N), Uj(J), Uij(N), xi_z, xi_c;
  int NN_z, NN_c;
  mat v_omega, omega, theta; 
  
  // starting values
  zj_pg = regspace(0, J - 1);
  cij.fill(0);
  Update_zj_pg_rep(zj_pg_rep, zj_pg, group);
  
  
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
    
    Update_Uj(Uj, xi_z, NN_z, zj_pg, kappa);
    Update_Uij(Uij, xi_c, NN_c, cij, kappa);
    Update_Distributional_Sticks(pi, pi_v, zj_pg, NN_z, alpha);
    Update_Observational_Sticks(omega, v_omega, cij, zj_pg_rep, NN_c, NN_z, beta);
    Update_theta(theta, Y, cij, a0, b0, k0, m0, NN_c, J);
    Update_Zj(zj_pg, Uj, xi_z, xi_c, pi, cij, omega, group, NN_z, J);
    Update_zj_pg_rep(zj_pg_rep, zj_pg, group);
    Update_Cij(cij, Y, Uij, xi_c, omega, zj_pg_rep, theta, N, NN_c);
    update_alpha(alpha, zj_pg, p_alpha_1, p_alpha_2);
    update_beta(beta, zj_pg, zj_pg_rep, cij, v_omega, p_beta_1, p_beta_2);
    
    // save the results
    if((iter >= nburn) & ((iter + 1) % thin == 0)){
      res_clust.row(res_index) = cij.t();
      res_nested.row(res_index) = zj_pg.t();
      if(eval_density){
        dens_out.slice(res_index) = eval_dens(grid, theta, omega, zj_pg, NN_c, J);
      }
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
  return res;
}

