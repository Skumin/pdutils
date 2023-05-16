// Code based on the "DeLong placements" C++ code from the pROC R package
#include <Rcpp.h>
using namespace Rcpp;

bool _cmp(std::pair<int, double> l, std::pair<int, double> r)
{
  return l.second < r.second;
}

// [[Rcpp::export]]
List delong_placements(NumericVector pds_def, NumericVector pds_nondef)
{

  int i, j, k, m, n, mdupl, ndupl, L;

  std::vector<double> cases = as<std::vector<double>>(pds_def);
  std::vector<double> controls = as<std::vector<double>>(pds_nondef);

  m = cases.size();
  n = controls.size();
  L = m + n;

  // concatenate cases and controls into a vector of L pairs of the form
  // (index, value), also save class labels (1 for cases, 0 for controls)
  std::vector<std::pair<int, double>> Z;
  std::vector<bool> labels;
  for (i = 0; i < m; i++)
  {
    Z.push_back(std::pair<int, double>(i, cases.at(i)));
    labels.push_back(true);
  }
  Rcpp::checkUserInterrupt();
  for (j = 0; j < n; j++)
  {
    Z.push_back(std::pair<int, double>(m + j, controls.at(j)));
    labels.push_back(false);
  }
  Rcpp::checkUserInterrupt();

  // sort Z from smallest to largest value, so Z holds the order indices and
  // order statistics of all classifiers
  std::sort(Z.begin(), Z.end(), _cmp);
  Rcpp::checkUserInterrupt();

  // the following calculates the "Delong-placements" X and Y in a single pass
  // over the vector Z, instead of having to double loop over all pairs of
  // (X_i, Y_j)
  std::vector<double> XY(L, 0.0);  // vector to hold the unnormalised X and Y values
  std::vector<int> X_inds, Y_inds; // temporary vectors to save indices of duplicates
  m = n = i = 0;                   // initialisation
  while (i < L)
  {
    X_inds.clear();
    Y_inds.clear();
    mdupl = ndupl = 0;
    if (i % 10000 == 0)
      Rcpp::checkUserInterrupt();
    while (1)
    {
      j = Z.at(i).first;
      if (labels.at(j))
      {
        mdupl++;
        X_inds.push_back(j);
      }
      else
      {
        ndupl++;
        Y_inds.push_back(j);
      }
      if (i == L - 1)
      {
        break;
      }
      if (Z.at(i).second != Z.at(i + 1).second)
      {
        break;
      }
      i++;
    }
    for (k = 0; k < mdupl; k++)
    {
      XY.at(X_inds.at(k)) = n + ndupl / 2.0;
    }
    for (k = 0; k < ndupl; k++)
    {
      XY.at(Y_inds.at(k)) = m + mdupl / 2.0;
    }
    n += ndupl;
    m += mdupl;
    i++;
  }

  double sum = 0.0;
  std::vector<double> X, Y;
  Rcpp::checkUserInterrupt();

  for (i = 0; i < L; i++)
  {
    if (labels.at(i))
    {
      sum += XY.at(i);
      X.push_back(XY.at(i) / n);
    }
    else
    {
      Y.push_back(1.0 - XY.at(i) / m);
    }
  }

  List ret;
  ret["theta"] = sum / m / n;
  ret["X"] = X;
  ret["Y"] = Y;
  return (ret);
}
