// This is a modified copy of fastglm/R/fit_glm.R

#include "ref_fastglm_fit_glm.h"
#include "ref_fastglm_fit_glm_dense.h"

using namespace Rcpp;

extern "C"
{
  SEXP logit_link(SEXP mu);
}

NumericVector linkfun_gaussian(const NumericVector &mu)
{
  return mu;
}

NumericVector linkfun_binomial(const NumericVector &mu)
{
  SEXP res = logit_link(wrap(mu));
  return as<NumericVector>(res);
}

NumericVector linkfun_poisson(const NumericVector &mu)
{
  return Rcpp::log(mu);
}

List fastglm(NumericMatrix x,
             SEXP y,
             std::string family,
             Nullable<NumericVector> start,
             Nullable<NumericVector> weights,
             Nullable<NumericVector> offset,
             Nullable<NumericVector> etastart,
             Nullable<NumericVector> mustart,
             int method,
             double tol,
             int maxit)
{
  // Determine number of observations and whether y is a matrix.
  int nobs;
  bool yIsMatrix = false;
  NumericMatrix yMat = as<NumericMatrix>(y);
  NumericVector yVec;
  if (yMat.ncol() > 1)
  {
    nobs = yMat.nrow();
    yIsMatrix = true;
  }
  else
  {
    yVec = yMat(_, 0);
    nobs = yVec.size();
  }

  // Set default weights (rep(1, nobs)) and offset (rep(0, nobs)) if not provided.
  NumericVector wt = weights.isNotNull() ? as<NumericVector>(weights)
                                         : NumericVector(nobs, 1.0);
  NumericVector off = offset.isNotNull() ? as<NumericVector>(offset)
                                         : NumericVector(nobs, 0.0);

  // Save a copy of mustart if provided.
  Nullable<NumericVector> mukeep = mustart;

  // Variables to hold computed "mustart" and an auxiliary vector n.
  NumericVector must;
  NumericVector n;

  if (family == "binomial")
  {
    if (!yIsMatrix)
    {
      // Single-column y: assume numeric (ignoring factor conversion).
      n = NumericVector(nobs, 1.0);
      // For observations with zero weight, set y to 0.
      for (int i = 0; i < nobs; i++)
      {
        if (wt[i] == 0)
          yVec[i] = 0;
      }
      // Check that y is in [0,1].
      for (int i = 0; i < nobs; i++)
      {
        if (yVec[i] < 0 || yVec[i] > 1)
          stop("y values must be 0 <= y <= 1");
      }
      must = NumericVector(nobs);
      for (int i = 0; i < nobs; i++)
      {
        must[i] = (wt[i] * yVec[i] + 0.5) / (wt[i] + 1);
      }
      // Warn if weighted successes are not integer.
      for (int i = 0; i < nobs; i++)
      {
        double m_val = wt[i] * yVec[i];
        if (std::abs(m_val - std::round(m_val)) > 0.001)
          warning("non-integer #successes in a binomial glm!");
      }
    }
    else
    {
      // y is a two-column matrix.
      if (yMat.ncol() != 2)
        stop("For binomial family, y must have 1 or 2 columns");
      // Warn if counts are non-integer.
      for (int i = 0; i < yMat.nrow(); i++)
      {
        for (int j = 0; j < yMat.ncol(); j++)
        {
          if (std::abs(yMat(i, j) - std::round(yMat(i, j))) > 0.001)
            warning("non-integer counts in a binomial glm!");
        }
      }
      NumericVector y1 = yMat(_, 0); // first column
      NumericVector y2 = yMat(_, 1); // second column
      n = NumericVector(nobs);
      for (int i = 0; i < nobs; i++)
      {
        n[i] = y1[i] + y2[i];
      }
      yVec = NumericVector(nobs);
      for (int i = 0; i < nobs; i++)
      {
        if (n[i] == 0)
          yVec[i] = 0;
        else
          yVec[i] = y1[i] / n[i];
      }
      // Adjust weights.
      for (int i = 0; i < nobs; i++)
      {
        wt[i] = wt[i] * n[i];
      }
      must = NumericVector(nobs);
      for (int i = 0; i < nobs; i++)
      {
        must[i] = (n[i] * yVec[i] + 0.5) / (n[i] + 1);
      }
    }
  }
  else if (family == "poisson")
  {
    if (yIsMatrix)
      stop("Poisson family expects y to be a vector");
    for (int i = 0; i < nobs; i++)
    {
      if (yVec[i] < 0)
        stop("negative values not allowed for the 'Poisson' family");
    }
    n = NumericVector(nobs, 1.0);
    must = NumericVector(nobs);
    for (int i = 0; i < nobs; i++)
    {
      must[i] = yVec[i] + 0.1;
    }
  }
  else if (family == "gaussian")
  {
    if (yIsMatrix)
      stop("Gaussian family expects y to be a vector");
    n = NumericVector(nobs, 1.0);
    must = NumericVector(nobs);
    for (int i = 0; i < nobs; i++)
    {
      must[i] = yVec[i];
    }
  }
  else
  {
    stop("Unsupported family");
  }

  // If mustart was provided, override the computed value.
  if (mukeep.isNotNull())
    must = as<NumericVector>(mukeep);

  // Define the link function as a lambda.
  // For gaussian, the identity; for binomial, call logit_link_cpp; for poisson, log.
  NumericVector (*linkfun)(const NumericVector &);
  if (family == "gaussian")
  {
    linkfun = &linkfun_gaussian;
  }
  else if (family == "binomial")
  {
    linkfun = &linkfun_binomial;
  }
  else if (family == "poisson")
  {
    linkfun = &linkfun_poisson;
  }
  else
  {
    stop("Unsupported family");
  }

  // Ensure y is numeric: if y was a matrix, we now use the processed yVec.
  // Otherwise, yVec is already set.

  // Compute eta: use etastart if provided; else, if start is provided compute offset + x %*% start;
  // otherwise, use linkfun(must).
  NumericVector eta;
  if (etastart.isNotNull())
  {
    eta = as<NumericVector>(etastart);
  }
  else if (start.isNotNull())
  {
    NumericVector startVec = as<NumericVector>(start);
    if (startVec.size() != x.ncol())
      stop("Length of start must equal number of columns of x");
    int nrow = x.nrow();
    NumericVector xstart(nrow, 0.0);
    for (int i = 0; i < nrow; i++)
    {
      double sum = 0.0;
      for (int j = 0; j < x.ncol(); j++)
      {
        sum += x(i, j) * startVec[j];
      }
      xstart[i] = sum;
    }
    eta = off; // start with offset
    for (int i = 0; i < nrow; i++)
    {
      eta[i] += xstart[i];
    }
  }
  else
  {
    eta = linkfun(must);
  }

  // If start is not provided, default to a zero vector of appropriate length.
  NumericVector startVec;
  if (start.isNotNull())
  {
    startVec = as<NumericVector>(start);
  }
  else
  {
    startVec = NumericVector(x.ncol(), 0.0);
  }

  // Finally, call fastglm with the processed arguments.
  // Note: we use yVec as the final response.
  List res = fit_glm(x, yVec, wt, off, startVec, eta, method, tol, maxit, family);
  return res;
}
