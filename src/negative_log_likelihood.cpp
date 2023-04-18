#include "fastcpd.h"

//' Solve logistic/poisson regression using Gradient Descent Extension to the
//' multivariate case
//' This function is not meant to be called directly by the user.
//'
//' @param data A data frame containing the data to be segmented.
//' @param theta Estimate of the parameters. If null, the function will estimate
//'   the parameters.
//' @param family Family of the model.
//' @param lambda Lambda for L1 regularization. Only used for lasso.
//' @param cv Whether to perform cross-validation to find the best lambda.
//' @keywords internal
//'
//' @return Negative log likelihood of the corresponding data with the given
//'   family.
// [[Rcpp::export]]
Rcpp::List negative_log_likelihood(
    arma::mat data,
    Rcpp::Nullable<arma::colvec> theta,
    std::string family,
    double lambda,
    bool cv = false
) {
    if (theta.isNull() && family == "lasso" && cv) {
        Rcpp::Environment glmnet = Rcpp::Environment::namespace_env("glmnet");
        Rcpp::Function cv_glmnet = glmnet["cv.glmnet"], predict_glmnet = glmnet["predict.glmnet"];
        Rcpp::List out = cv_glmnet(
            data.cols(1, data.n_cols - 1),
            data.col(0),
            Rcpp::Named("family") = "gaussian"
        );
        Rcpp::S4 out_coef = predict_glmnet(
            out["glmnet.fit"],
            Rcpp::Named("s") = out["lambda.1se"],
            Rcpp::Named("type") = "coefficients",
            Rcpp::Named("exact") = false
        );
        arma::vec glmnet_i = Rcpp::as<arma::vec>(out_coef.slot("i"));
        arma::vec glmnet_x = Rcpp::as<arma::vec>(out_coef.slot("x"));
        arma::vec par = arma::zeros(data.n_cols - 1);
        for (unsigned int i = 1; i < glmnet_i.n_elem; i++) {
            par(glmnet_i(i) - 1) = glmnet_x(i);
        }
        return Rcpp::List::create(Rcpp::Named("par") = par,
                                  Rcpp::Named("value") = R_NilValue);
    } else if (theta.isNull() && family == "lasso" && !cv) {
        Rcpp::Environment stats = Rcpp::Environment::namespace_env("stats"), glmnet = Rcpp::Environment::namespace_env("glmnet");
        Rcpp::Function deviance = stats["deviance"], glmnet_ = glmnet["glmnet"], predict_glmnet = glmnet["predict.glmnet"];
        Rcpp::List out = glmnet_(
            data.cols(1, data.n_cols - 1),
            data.col(0),
            Rcpp::Named("family") = "gaussian",
            Rcpp::Named("lambda") = lambda
        );
        Rcpp::S4 out_par = out["beta"];
        arma::vec par_i = Rcpp::as<arma::vec>(out_par.slot("i"));
        arma::vec par_x = Rcpp::as<arma::vec>(out_par.slot("x"));
        arma::vec par = arma::zeros(data.n_cols - 1);
        for (unsigned int i = 0; i < par_i.n_elem; i++) {
            par(par_i(i)) = par_x(i);
        }
        double value = Rcpp::as<double>(deviance(out));
        arma::vec fitted_values = Rcpp::as<arma::vec>(predict_glmnet(out, data.cols(1, data.n_cols - 1), Rcpp::Named("s") = lambda));
        arma::vec residuals = data.col(0) - fitted_values;
        return Rcpp::List::create(Rcpp::Named("par") = par,
                                  Rcpp::Named("value") = value / 2,
                                  Rcpp::Named("residuals") = residuals);
    } else if (theta.isNull()) {
        // Estimate theta in binomial/poisson/gaussian family
        arma::mat x = data.cols(1, data.n_cols - 1);
        arma::vec y = data.col(0);
        Rcpp::Environment fastglm = Rcpp::Environment::namespace_env("fastglm");
        Rcpp::Function fastglm_ = fastglm["fastglm"];
        Rcpp::List out = fastglm_(x, y, family);
        arma::vec par = Rcpp::as<arma::vec>(out["coefficients"]);
        arma::vec residuals = Rcpp::as<arma::vec>(out["residuals"]);
        double value = out["deviance"];
        return Rcpp::List::create(Rcpp::Named("par") = par,
                                  Rcpp::Named("value") = value / 2,
                                  Rcpp::Named("residuals") = residuals);
    } else if (family == "lasso" || family == "gaussian") {

        arma::colvec theta_nonnull = Rcpp::as<arma::colvec>(theta);
        // Calculate negative log likelihood in gaussian family
        arma::vec y = data.col(0);
        arma::mat x = data.cols(1, data.n_cols - 1);
        double penalty = lambda * arma::accu(arma::abs(theta_nonnull));
        return Rcpp::List::create(Rcpp::Named("value") = arma::accu(arma::pow(y - x * theta_nonnull, 2)) / 2 + penalty);

    } else if (family == "binomial") {

        // Calculate negative log likelihood in binomial family
        arma::colvec theta_nonnull = Rcpp::as<arma::colvec>(theta);
        arma::vec y = data.col(0);
        arma::mat x = data.cols(1, data.n_cols - 1);
        arma::colvec u = x * theta_nonnull;
        return Rcpp::List::create(Rcpp::Named("value") = arma::accu(-y % u + arma::log(1 + arma::exp(u))));
    } else {
        // Calculate negative log likelihood in poisson family
        arma::colvec theta_nonnull = Rcpp::as<arma::colvec>(theta);
        arma::vec y = data.col(0);
        arma::mat x = data.cols(1, data.n_cols - 1);
        arma::colvec u = x * theta_nonnull;
        Rcpp::NumericVector y_factorial = Rcpp::lfactorial(Rcpp::wrap(y));
        return Rcpp::List::create(Rcpp::Named("value") = arma::accu(-y % u + arma::exp(u) + arma::vec(y_factorial.begin(), y_factorial.size(), false)));
    }
}
