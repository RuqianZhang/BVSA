#include <RcppArmadillo.h>
#include <iostream>
#include <Rcpp/Benchmark/Timer.h>
#include <stdio.h>

using namespace Rcpp;
using namespace std;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

vec vec_invgamma(int n, double a, double b);
vec vec_binom(vec prob);
mat mvrnormC(int n, vec mu, mat sigma);
vec rPolyaGamma(int num, double h, double z);
uvec find_active_except(vec indices, int j);

// [[Rcpp::export]]
List BVSA_C(mat X, vec Y, vec T, mat Z,
	int iter=2000, int burn=1000, double tau_gamma0=0.05, double tau_gamma1=1, double q_gamma=0.2,
	double tau_beta0=0.05, double tau_beta1=1, double q_beta=0.1, int max_size=30,
	double a0=2, double b0=1, double sigmasq_alpha=1){
	int N = X.n_rows, q_X = X.n_cols-1, q_Z = Z.n_cols-1;

	mat trace_gamma(q_X+1, iter);
	mat trace_I_gamma(q_X+1, iter);
	mat trace_beta(q_Z+1, iter);
	mat trace_I_beta(q_Z+1, iter);
	mat trace_alpha(2, iter);
	mat trace_delta(N, iter);
  mat trace_sigma(1, iter);
  
	// Initiate by priors
	// sigmasq
	double sigmasq = vec_invgamma(1, a0, b0)(0);
	// I_beta
	vec p_beta = join_cols(linspace(1,1,1), linspace(q_beta, q_beta, q_Z));
	vec I_beta = vec_binom(p_beta);
	// beta
	vec mu_beta = linspace(0, 0, q_Z+1);
	vec cor_beta_diag(q_Z+1);
	for (int i = 0; i < q_Z+1; ++i){
	  cor_beta_diag(i) = I_beta(i)*sigmasq*pow(tau_beta1,2)
    +(1-I_beta(i))*sigmasq*pow(tau_beta0,2);
	};
	mat cor_beta = diagmat(cor_beta_diag);
	mat beta = mvrnormC(1, mu_beta, cor_beta).t();
  // alpha
  vec alpha(2);
  alpha(0) = R::rnorm(0, sqrt(sigmasq*sigmasq_alpha));
  alpha(1) = R::rnorm(10, sqrt(sigmasq*sigmasq_alpha));
  // I_gamma
  vec p_gamma = join_cols(linspace(1,1,1), linspace(q_gamma, q_gamma, q_X));
  vec I_gamma = vec_binom(p_gamma);
  
  // gamma
  vec mu_gamma = linspace(0, 0, q_X+1);
  vec cor_gamma_diag(q_X+1);
  for (int i = 0; i < q_X+1; ++i){
    cor_gamma_diag(i) = I_gamma(i)*pow(tau_gamma1,2)
    +(1-I_gamma(i))*pow(tau_gamma0,2);
  };
  mat cor_gamma = diagmat(cor_gamma_diag);
  mat gamma = mvrnormC(1, mu_gamma, cor_gamma).t();
  // delta
  vec p_delta(N);
  vec logit_lin = exp(-X*diagmat(I_gamma)*gamma);
  for (uword i = 0; i < N; ++i){
    p_delta(i) = 1/(1+logit_lin(i));
  };
  vec delta = vec_binom(p_delta);
  // polya-gamma
  vec w = rPolyaGamma(N, 1, 0);
  
  
  // Gibbs update
  vec kappa, m_betaA, m_gammaA, tmp_prob;
  uvec index_gamma, index_beta, index_ex;
  mat W, D, ID, Z_A, X_A, V_betaA, V_gammaA, beta_A, gamma_A;
  int k, l;
  double a, b, R, m_1, m_2, V_1, V_2, prob_threshold;

  Timer timer;
  
  for (int i = 0; i < iter; ++i){
    
    timer.step("start");
    
    kappa = (delta-0.5)/w;
    W = diagmat(w);
    k = sum(I_gamma);
    l = sum(I_beta);
    
    // update gamma
    index_gamma = find(I_gamma == 1);
    X_A = X.cols(index_gamma);
    V_gammaA = inv_sympd(X_A.t()*W*X_A
                               + diagmat(linspace(pow(tau_gamma1,-2), pow(tau_gamma1,-2), k)));
    m_gammaA = V_gammaA*X_A.t()*W*kappa;
    gamma_A = mvrnormC(1, m_gammaA, V_gammaA).t();
    gamma = mvrnormC(1, linspace(0, 0, q_X+1),
                    diagmat(linspace(pow(tau_gamma0, 2), pow(tau_gamma0, 2), q_X+1))).t();
    gamma.elem(index_gamma) = gamma_A;
    
    trace_gamma.col(i) = gamma;
    
    // update p_gamma & I_gamma
    for (int j = 1; j < q_X+1; ++j){
      index_ex = find_active_except(I_gamma, j);
      R = q_gamma*tau_gamma0/((1-q_gamma)*tau_gamma1)
        *sum(exp((kappa.t()-gamma(index_ex).t()*X.cols(index_ex).t())*W*X.col(j)*gamma(j)
            -0.5*(X.col(j).t()*W*X.col(j)+pow(tau_gamma1, -2)-pow(tau_gamma0, -2))*pow(gamma(j),2)));
      if (is_finite(R)) {p_gamma(j) = R/(1+R);} else {p_gamma(j) = 1;}
      I_gamma(j) = R::rbinom(1, p_gamma(j));
    };
    // if size > max_size, choose the covariates among active ones by ordering of the prob's.
    if (sum(I_gamma) > max_size+1) {
      I_gamma = linspace(0, 0, q_X+1);
      // get p(max_size), if p > p(max_size), I_gamma = 1
      tmp_prob = p_gamma;
      std::sort(tmp_prob.begin(), tmp_prob.end(), std::greater<double>());
      prob_threshold = tmp_prob(max_size);

      for (int j = 0; j < q_X+1; ++j) {
        if (p_gamma(j) >= prob_threshold) {
          I_gamma(j) = 1;
        }else {I_gamma(j) = 0;}
        
        if (sum(I_gamma) > max_size+1) {break;}
      }
    }
    
    trace_I_gamma.col(i) = I_gamma;
    
    timer.step("update gamma");
    
    
    // update W
    index_gamma = find(I_gamma == 1);
    X_A = X.cols(index_gamma);
    for (int j = 0; j < N; ++j){
      w(j) = rPolyaGamma(1, 1, sum(X_A.row(j)*gamma(index_gamma)))(0);
    }
    W = diagmat(w);
    
    timer.step("update delta");
    
    
    // update delta
    index_beta = find(I_beta == 1);
    Z_A = Z.cols(index_beta);
    for (int j = 0; j < N; ++j){
      logit_lin(j) = -1/(2*sigmasq)*(alpha(1)-2*Y(j)+2*sum(Z_A.row(j)*beta(index_beta))+alpha(0))
                   *T(j)*(alpha(0)-alpha(1))+sum(X_A.row(j)*gamma(index_gamma));
      p_delta(j) = 1/(1+exp(-logit_lin(j)));
    }
    delta = vec_binom(p_delta);
    trace_delta.col(i) = delta;
    
    
    // update beta
    D = diagmat(delta);
    ID = diagmat(linspace(1, 1, N))-D;
    V_betaA = inv_sympd(diagmat(linspace(1/sigmasq, 1/sigmasq, l))
                *(Z_A.t()*Z_A+diagmat(linspace(pow(tau_beta1,-2), pow(tau_beta1,-2), l))));
    m_betaA = V_betaA*diagmat(linspace(1/sigmasq, 1/sigmasq, l))
            *Z_A.t()*(Y-D*T*alpha(0)-ID*T*alpha(1));
    beta_A = mvrnormC(1, m_betaA, V_betaA).t();
    beta = mvrnormC(1, linspace(0, 0, q_Z+1),
                   diagmat(linspace(sigmasq*pow(tau_beta0,2), sigmasq*pow(tau_beta0,2), q_Z+1))).t();
    beta.elem(index_beta) = beta_A;
    trace_beta.col(i) = beta;
    
    // update alpha
    V_1 = 1/sigmasq*(sum(T.t()*D*T)+1/sigmasq_alpha);
    m_1 = 1/(V_1*sigmasq)*sum(T.t()*D*(Y-Z_A*beta(index_beta)-ID*T*alpha(1)));
    alpha(0) = R::rnorm(m_1, sqrt(1/V_1));
    V_2 = 1/sigmasq*(sum(T.t()*ID*T)+1/sigmasq_alpha);
    m_2 = 1/(V_2*sigmasq)*sum(T.t()*ID*(Y-Z_A*beta(index_beta)-D*T*alpha(0)));
    alpha(1) = R::rnorm(m_2, sqrt(1/V_2));
    trace_alpha.col(i) = alpha;
    
    // update p_beta & I_beta
    for (int j = 1; j < q_Z+1; ++j){
      index_ex = find_active_except(I_beta, j);
      R = q_beta*tau_beta0/((1-q_beta)*tau_beta1)
        *exp(sum(1/sigmasq*beta(j)*Z.col(j).t()
               *(Y-D*T*alpha(0)-ID*T*alpha(1)-Z.cols(index_ex)*beta(index_ex))
               -1/(2*sigmasq)*(Z.col(j).t()*Z.col(j)+pow(tau_beta1,-2)
                -pow(tau_beta0,-2))*pow(beta(j),2)));
      if (is_finite(R)) {p_beta(j) = R/(1+R);} else {p_beta(j) = 1;}
      I_beta(j) = R::rbinom(1, p_beta(j));
    };
    if (sum(I_beta) > max_size+1) {
      I_beta = linspace(0, 0, q_Z+1);
      // get p(max_size), if p > p(max_size), I_gamma = 1
      tmp_prob = p_beta;
      std::sort(tmp_prob.begin(), tmp_prob.end(), std::greater<double>());
      prob_threshold = tmp_prob(max_size);

      for (int j = 0; j < q_X+1; ++j) {
        if (p_beta(j) >= prob_threshold) {
          I_beta(j) = 1;
        }else {I_beta(j) = 0;}
        
        if (sum(I_beta) > max_size+1) {break;}
      }
    }
    
    trace_I_beta.col(i) = I_beta;
    
    timer.step("update beta");
    
    // update sigmasq
    index_beta = find(I_beta == 1);
    Z_A = Z.cols(index_beta);
    a = a0+N/2+(q_Z+1)/2+1;
    b = b0+0.5*sum(((Y-Z_A*beta(index_beta)-D*T*alpha(0)-ID*T*alpha(1)).t()
                  *(Y-Z_A*beta(index_beta)-D*T*alpha(0)-ID*T*alpha(1))
                  +1/sigmasq_alpha*(pow(alpha(0),2)+pow(alpha(1),2))
                  +beta.t()*diagmat(pow(tau_beta1,-2)*I_beta
                  +pow(tau_beta0,-2)*(linspace(1, 1, q_Z+1)-I_beta))*beta));
    sigmasq = vec_invgamma(1, a, b)(0);
    trace_sigma.col(i) = sigmasq;

    timer.step("end");
    
  }; // end Gibbs
  
  
  // burn-in and average
  mat burnin_I_gamma, burnin_I_beta, burnin_alpha, burnin_gamma, burnin_beta, burnin_delta, burnin_sigma;
  vec sigmasq_pred;
  
  burnin_I_gamma = trace_I_gamma.cols(burn, iter-1);
  p_gamma = burnin_I_gamma*linspace(1,1,iter-burn)/(iter-burn);
  I_gamma = linspace(0, 0, q_X+1);
  I_gamma(find(p_gamma >= 0.5)).ones();

  burnin_I_beta = trace_I_beta.cols(burn, iter-1);
  p_beta = burnin_I_beta*linspace(1,1,iter-burn)/(iter-burn);
  I_beta = linspace(0, 0, q_Z+1);
  I_beta(find(p_beta >= 0.5)).ones();

  burnin_alpha= trace_alpha.cols(burn, iter-1);
  alpha = burnin_alpha*linspace(1,1,iter-burn)/(iter-burn);
  burnin_gamma= trace_gamma.cols(burn, iter-1);  
  gamma = burnin_gamma*linspace(1,1,iter-burn)/(iter-burn);
  burnin_beta= trace_beta.cols(burn, iter-1);
  beta = burnin_beta*linspace(1,1,iter-burn)/(iter-burn);
  burnin_delta= trace_delta.cols(burn, iter-1);
  delta = burnin_delta*linspace(1,1,iter-burn)/(iter-burn);
  burnin_sigma= trace_sigma.cols(burn, iter-1);
  sigmasq_pred = burnin_sigma*linspace(1,1,iter-burn)/(iter-burn);

  if (alpha(0) <= alpha(1)){
    double tmp = alpha(0);
    alpha(0) = alpha(1);
    alpha(1) = tmp;
    delta = linspace(1, 1, N)-delta;
    gamma = -gamma;
  }

	return List::create(
		_["I_gamma"] = I_gamma,
		_["p_gamma"] = p_gamma,
    _["I_beta"] = I_beta,
    _["p_beta"] = p_beta,
    _["alpha"] = alpha,
    _["gamma"] = gamma,
    _["beta"] = beta,
    _["delta"] = delta,
	  _["timer"] = timer,
    _["trace_I_gamma"] = trace_I_gamma,
    _["trace_I_beta"] = trace_I_beta,
    _["trace_alpha"] = trace_alpha,
    _["trace_gamma"] = trace_gamma,
    _["trace_beta"] = trace_beta,
    _["trace_delta"] = trace_delta,
    _["sigmasq"] = sigmasq_pred,
    _["trace_sigmasq"] = trace_sigma
	);
	
}



// find the indices of active variables and remove certain one
uvec find_active_except(vec indices, int j){
  indices(j) = 0;
  uvec index = find(indices == 1);
  return index;
}


// random Polya-Gamma
vec rPolyaGamma(int num = 1, double h = 1, double z = 0){
  Function pg_sample("rpg", Environment::namespace_env("BayesLogit"));
  
  NumericVector tmp = pg_sample(num, h, z);
  vec res = tmp;
  return res;
}


// invgamma
vec vec_invgamma(int n, double a, double b){
  vec invgamma(n);
  for (int i = 0; i < n; ++i){
  	invgamma(i) = 1/R::rgamma(a, 1/b);
  };
  
  return invgamma;
}


// rbinom
vec vec_binom(vec prob){
  int dim = prob.size();
  vec out(dim);
  
  for (int i = 0; i < dim; ++i){
    out(i) = R::rbinom(1, prob(i));
  };

  return out;
}


// mvrnorm
mat mvrnormC(int n, vec mu, mat sigma) {
  int ncols = sigma.n_cols;
  mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}


