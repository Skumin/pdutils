// [[Rcpp::depends(RcppArmadillo)]]
#include <vector>
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

NumericMatrix linsolve(const arma::mat& A, const arma::mat& B) {
  arma::colvec coef = arma::solve(A, B);
  return Rcpp::wrap(coef);
}

NumericMatrix mmult(NumericMatrix A, NumericMatrix B) {
  arma::mat Am = Rcpp::as< arma::mat >(A);
  arma::mat Bm = Rcpp::as< arma::mat >(B);
  return Rcpp::wrap( Am * Bm );
}

IntegerVector sequenceVec(int n) {
  std::vector<int> ivec (n, 1);
  std::iota(ivec.begin(), ivec.end(), 0);
  return wrap(ivec);
}

IntegerVector eraseInd(IntegerVector ivec, int eras) {
  IntegerVector newvec = clone(ivec);
  std::vector<int> newvec1 = as<std::vector<int> >(newvec);
  newvec1.erase (newvec1.begin()+eras);
  return wrap(newvec1);
}

// [[Rcpp::export]]
NumericMatrix mutate_deleq(NumericMatrix mat, NumericMatrix boxbounds, double fParam) {
  NumericMatrix newmat(mat.nrow(), mat.ncol());
  int matsize = newmat.nrow();
  int parsize = newmat.ncol();
  IntegerVector rngg1 = sequenceVec(matsize);
  for (int i = 0; i < matsize; i++) {
    IntegerVector rngg = eraseInd(rngg1, i);
    int boundsok = 0;
    while (boundsok == 0) {
      IntegerVector rs = sample(rngg, 3, false);
      NumericVector num = mat(rs[0], _ ) + fParam * (mat(rs[1], _ ) - mat(rs[2], _ ));
      IntegerVector tst(parsize);
      for (int j = 0; j < parsize; j++) {
        if (num[j] >= boxbounds(j, 0) && num[j] <= boxbounds(j, 1)) {
          tst[j] = 1;
        }
      }
      if (sum(tst) == parsize) {
        boundsok = 1;
        newmat(i, _ ) = num;
      }
    }
  }
  return newmat;
}

// [[Rcpp::export]]
NumericMatrix project_population(NumericMatrix mat, NumericMatrix Emat, NumericMatrix constr) {
  NumericMatrix Mmat = mmult(Emat, transpose(Emat));
  NumericMatrix projmat(mat.nrow(), mat.ncol());
  for (int i = 0; i < mat.nrow(); i++) {
    NumericMatrix num(mat.ncol(), 1);
    num( _ , 0) = mat(i, _ );
    NumericMatrix tmp = mmult(Emat, num);
    NumericVector ztmp = tmp( _ , 0) - constr( _ , 0);
    NumericMatrix z(ztmp.size(), 1);
    z( _ , 0) = ztmp;
    NumericMatrix u = linsolve(Rcpp::as< arma::mat >(Mmat), Rcpp::as< arma::mat >(z));
    NumericMatrix v = mmult(transpose(Emat), u);
    projmat(i, _ ) = num( _ , 0) - v( _ , 0);
  }
  return projmat;
}

// [[Rcpp::export]]
NumericMatrix gen_init_pop(int NP, NumericMatrix boxbounds, NumericMatrix Emat, NumericMatrix constr) {
  int dm = boxbounds.nrow();
  NumericMatrix xi(NP, dm);
  NumericMatrix Mmat = mmult(Emat, transpose(Emat));
  NumericMatrix y = linsolve(Rcpp::as< arma::mat >(Mmat), Rcpp::as< arma::mat >(constr));
  NumericMatrix x0 = mmult(transpose(Emat), y);
  for (int i = 0; i < NP; i++) {
    int boundsok = 0;
    while (boundsok == 0) {
      NumericMatrix d(dm, 1);
      for (int j = 0; j < dm; j++) {
        d(j, 0) = R::runif(0, 1) * (boxbounds(j, 1) - boxbounds(j, 0)) + boxbounds(j, 0);
      }
      d( _ , 0) = d( _ , 0)/sum(d( _ , 0));
      NumericMatrix z = mmult(Emat, d);
      NumericMatrix u = linsolve(Rcpp::as< arma::mat >(Mmat), Rcpp::as< arma::mat >(z));
      NumericMatrix v = mmult(transpose(Emat), u);
      NumericVector num = x0( _ , 0) + d( _ , 0) - v( _ , 0);
      IntegerVector ids(dm);
      for (int j = 0; j < dm; j++) {
        if (num[j] >= boxbounds(j, 0) && num[j] <= boxbounds(j, 1)) {
          ids[j] = 1;
        }
      }
      if (sum(ids) == ids.size()) {
        xi(i, _ ) = num;
        boundsok = 1;
      }
    }
  }
  return xi;
}

// [[Rcpp::export]]
NumericMatrix gen_init_pop_x0(int NP, NumericMatrix boxbounds, NumericMatrix Emat, NumericMatrix constr, NumericMatrix x0) {
  int dm = boxbounds.nrow();
  NumericMatrix xi(NP, dm);
  NumericMatrix Mmat = mmult(Emat, transpose(Emat));
  for (int i = 0; i < NP; i++) {
    int boundsok = 0;
    while (boundsok == 0) {
      NumericMatrix d(dm, 1);
      for (int j = 0; j < dm; j++) {
        d(j, 0) = R::runif(0, 1) * (boxbounds(j, 1) - boxbounds(j, 0)) + boxbounds(j, 0);
      }
      d( _ , 0) = d( _ , 0)/sum(d( _ , 0));
      NumericMatrix z = mmult(Emat, d);
      NumericMatrix u = linsolve(Rcpp::as< arma::mat >(Mmat), Rcpp::as< arma::mat >(z));
      NumericMatrix v = mmult(transpose(Emat), u);
      NumericVector num = x0( _ , 0) + d( _ , 0) - v( _ , 0);
      IntegerVector ids(dm);
      for (int j = 0; j < dm; j++) {
        if (num[j] >= boxbounds(j, 0) && num[j] <= boxbounds(j, 1)) {
          ids[j] = 1;
        }
      }
      if (sum(ids) == ids.size()) {
        xi(i, _ ) = num;
        boundsok = 1;
      }
    }
  }
  return xi;
}
