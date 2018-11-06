#define ARMA_64BIT_WORD 1

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double ll_pen_nb(arma::vec beta, arma::sp_mat X, arma::vec y, arma::vec offset,
		double theta, double lambda, arma::sp_mat S, double ll_factor,
		double lambda_factor, int n) {
  // theta must be a matrix to use lgamma on it
  arma::mat theta_mat(1,1);
  theta_mat.fill(theta);

  // mu = exp(eta) = exp(offset + X * beta)
  arma::vec lin = vectorise(X * beta);
  arma::vec eta = offset + lin;
  arma::vec mu = exp(eta);

  // some other parts
  arma::vec aux1 = theta + y;
  arma::vec aux2 = theta + mu;
  arma::vec lga = lgamma(aux1);
  arma::mat lgt = lgamma(theta_mat);
  double lgts = as_scalar(lgt);
  arma::vec lf = lgamma(y + 1);

  // first part of equation with the gamma and factorial fucntion
  double ls = sum(lga - (lf + lgts));

  // second part of the equation
  arma::vec rs = (trans(y) * eta) + (n * theta * log(theta)) - (trans(aux1) * log(aux2));
  double rss = as_scalar(rs);

  // penalization term
  arma::mat pen = (trans(beta) * S) * beta;
  double pens = as_scalar(pen);

  // all together normalized by a factor
  double res = (-1/ll_factor) * (ls + rss) + ((1/lambda_factor) * (lambda * pens));
  return res;
}

// [[Rcpp::export]]
arma::vec gr_ll_pen_nb(arma::vec beta, arma::sp_mat X, arma::sp_mat XT, arma::vec y,
			   arma::vec offset, double theta,
			   double lambda, arma::sp_mat S) {
  
  // mu = exp(eta) = exp(offset + X * beta)
  arma::vec lin = vectorise(X * beta);
  arma::vec eta = offset + lin;
  arma::vec mu = exp(eta);
  arma::vec z = (y - mu)/(1 + (mu/theta));
  arma::vec pen = S * beta;
  arma::vec gr = XT *z;
  arma::vec res = ((-1) * gr) + (2 * lambda * pen);
  return res;
}

// [[Rcpp::export]]
arma::vec negbin_hessian(arma::vec y, arma::vec mu, double theta) {
  arma::vec ls = mu % (1 + y/theta);
  arma::vec rs = square(1 + mu/theta);
  arma::vec dd = ls/rs;
  return dd;
}

// [[Rcpp::export]]
double ll_pen_qbd(arma::vec beta, arma::sp_mat X, arma::vec y, arma::vec offset,
		     double theta, double lambda, arma::sp_mat S, double ll_factor,
		     double lambda_factor, int n) {
  
  arma::vec lin = vectorise(X * beta);
  arma::vec eta = offset + lin;
  arma::vec mu = exp(eta)/(exp(eta) + 1);

  arma::vec a = log(mu);
  arma::vec b = -log(lgamma(y)) - log(lgamma(1 - y));
  arma::vec aux = mu + y * theta;
  arma::vec c = log(aux) * (y - 1);
  arma::vec d = (1 - y) * (log(1 - aux));
  double res = accu(a + b + c + d);
  return res;
}

// [[Rcpp::export]]
arma::sp_mat compute_pen_hessian(arma::vec beta, arma::sp_mat X, arma::sp_mat XT, arma::vec offset,
				 arma::vec y, arma::sp_mat S, double lambda, double theta, int hessid) {
  arma::vec lin = vectorise(X * beta);
  arma::vec eta = offset + lin;
  arma::vec mu = exp(eta);

  // 1 = negative binomial
  // 2 = quasi-binomial
  arma::vec dd = ones<vec>(size(y));

  if(hessid == 1) {
    dd = negbin_hessian(y, mu, theta);
  }

  arma::sp_mat D = diagmat(sp_mat(dd));

  arma::sp_mat H = XT * D * X;
  arma::sp_mat pen = 2 * lambda * S;
  arma::sp_mat res = H + pen;
  return res;
}

// [[Rcpp::export]]
arma::sp_mat compute_stdError(arma::sp_mat X, arma::sp_mat H) {
  arma::sp_mat V = (X * H) % X;
  arma::sp_mat Diag = sum(V, 1);
  arma::sp_mat res = sqrt(Diag);
  return res;
}
