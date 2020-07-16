// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppEigen)]]

#include <Rcpp.h>
#include <RcppEigen.h>

#include <Eigen/Dense>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <cmath>
#include <numeric>
#include <chrono>

using std::begin;
using std::end;
using boost::math::tools::bisect;
using boost::math::binomial_coefficient;

namespace Eigen
  {
  namespace internal
  {
    template<typename Scalar>
    struct scalar_normal_dist_op
    {
      static boost::mt19937 rng;
      mutable boost::normal_distribution<Scalar> norm;

      EIGEN_EMPTY_STRUCT_CTOR(scalar_normal_dist_op)

        template<typename Index>
        inline const Scalar operator() (Index, Index = 0) const { return norm(rng); }
    };

    template<typename Scalar> boost::mt19937 scalar_normal_dist_op<Scalar>::rng;

    template<typename Scalar>
    struct functor_traits<scalar_normal_dist_op<Scalar> >
    {
      enum {Cost = 50 * NumTraits<Scalar>::MulCost, PacketAccess = false, IsRepeatable = false};
    };
  }
}

double pt_prob_less_k(int nums, int defs, double pd, double rho, Eigen::MatrixXd reals)
{
  Eigen::MatrixXd piSt_one_per(reals.rows(), reals.cols());
  boost::math::normal n_dist;

  double q_pd = quantile(n_dist, pd);
  double sqrt_rho = sqrt(rho);
  double sqrt_om_rho = sqrt(1. - rho);

  for (int j = 0; j < reals.cols(); j++)
  {
    for (int i = 0; i < reals.rows(); i++)
    {
      piSt_one_per(i, j) = 1. - cdf(n_dist, (q_pd - sqrt_rho * reals(i, j)) / sqrt_om_rho);
    }
  }

  Eigen::VectorXd piSt = Eigen::VectorXd::Ones(piSt_one_per.rows());
  piSt -= piSt_one_per.rowwise().prod();

  Eigen::MatrixXd s1 = Eigen::MatrixXd::Zero(reals.rows(), defs + 1);

  double binom_coeff;

  for (int j = 0; j < defs + 1; j++)
  {
    binom_coeff = boost::math::binomial_coefficient<double>(nums, j);
    for (int i = 0; i < piSt.rows(); i++)
    {
      s1(i, j) = binom_coeff * pow(piSt(i), j) * pow(1. - piSt(i), nums - j);
    }
  }

  Eigen::VectorXd s2 = s1.rowwise().sum();
  return s2.mean();
}

Eigen::MatrixXd mv_normal(Eigen::MatrixXd corM, int simulations)
{
  int size = corM.rows();
  Eigen::internal::scalar_normal_dist_op<double> randN;
  Eigen::internal::scalar_normal_dist_op<double>::rng.seed((unsigned int)std::chrono::high_resolution_clock::now().time_since_epoch().count());

  Eigen::MatrixXd normTransform(size, size);
  Eigen::LLT<Eigen::MatrixXd> cholSolver(corM);

  // We can only use the cholesky decomposition if
  // the covariance matrix is symmetric, pos-definite.
  // But a covariance matrix might be pos-semi-definite.
  // In that case, we'll go to an EigenSolver
  if (cholSolver.info() == Eigen::Success)
  {
    // Use cholesky solver
    normTransform = cholSolver.matrixL();
  }
  else
  {
    // Use eigen solver
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(corM);
    normTransform = eigenSolver.eigenvectors() * eigenSolver.eigenvalues().cwiseSqrt().asDiagonal();
  }

  Eigen::MatrixXd samples = (normTransform * Eigen::MatrixXd::NullaryExpr(size, simulations, randN));

  return samples.transpose();
}

// [[Rcpp::export]]
double pt_multi_pd(int nums, int defs, double rho, double tau, int pers, double ci, int simulations)
{
  Eigen::MatrixXd corM(pers, pers);

  for (int j = 0; j < pers; j++)
  {
    for (int i = 0; i < pers; i++)
    {
      corM(i, j) = pow(tau, abs(i - j));
    }
  }

  Eigen::MatrixXd rSt = mv_normal(corM, simulations);

  double x0 = 0.0000001;
  double x1 = 0.9999999;
  double x, f, e;
  double f0 = -1. * pt_prob_less_k(nums, defs, x0, rho, rSt) - ci + 1.;

  e = 0.0001;

  do
  {
    x = (x0 + x1) / 2.;
    f = -1. * pt_prob_less_k(nums, defs, x, rho, rSt) - ci + 1.;

    if(f0 * f < 0.)
    {
      x1 = x;
    }
    else
    {
      x0 = x;
    }
  } while (fabs(f) > e);

  return x;
}

// [[Rcpp::export]]
Eigen::VectorXd pt_multi_pd_full(Eigen::VectorXd nums, Eigen::VectorXd defs, double rho, double tau, int pers, double ci, int simulations)
{
  Eigen::MatrixXd corM(pers, pers);

  for (int j = 0; j < pers; j++)
  {
    for (int i = 0; i < pers; i++)
    {
      corM(i, j) = pow(tau, abs(i - j));
    }
  }

  Eigen::MatrixXd rSt = mv_normal(corM, simulations);

  Eigen::VectorXd nums_all(nums.size());
  Eigen::VectorXd defs_all(nums.size());
  Eigen::VectorXd roots(nums.size());

  nums_all(0) = nums(0);
  defs_all(0) = defs(0);

  for (int i = 1; i < nums.size(); i++)
  {
    nums_all(i) = nums_all(i - 1) + nums(i);
    defs_all(i) = defs_all(i - 1) + defs(i);
  }

  Eigen::VectorXd nums_all_rev = nums_all.reverse();
  Eigen::VectorXd defs_all_rev = defs_all.reverse();

  double x0, x1;
  double x, f, f0;
  double e = 0.00001;

  int nms, def;

  for (int i = 0; i < nums.size(); i++)
  {
    nms = nums_all_rev(i);
    def = defs_all_rev(i);

    x0 = 0.0000001;
    x1 = 0.9999999;
    f0 = -1. * pt_prob_less_k(nms, def, x0, rho, rSt) - ci + 1.;

    do
    {
      x = (x0 + x1) / 2.;
      f = -1. * pt_prob_less_k(nms, def, x, rho, rSt) - ci + 1.;

      if(f0 * f < 0.)
      {
        x1 = x;
      }
      else
      {
        x0 = x;
      }
    } while (fabs(f) > e);

    roots(i) = x;
  }

  return roots.reverse();
}
