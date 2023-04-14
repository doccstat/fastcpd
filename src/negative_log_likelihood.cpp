#include "fastcpd.h"

//' Solve logistic/poisson regression using Gradient Descent Extension to the
//' multivariate case
//'
//' @param data A data frame containing the data to be segmented.
//' @param theta Estimate of the parameters. If null, the function will estimate
//'   the parameters.
//' @param family Family of the model.
//' @param lambda Lambda for L1 regularization. Only used for lasso.
//' @param cv Whether to perform cross-validation to find the best lambda.
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
    if (theta.isNull() && family == "lasso") {
        if (cv) {
            Rcpp::Environment glmnet("package:glmnet");
            Rcpp::Function cv_glmnet = glmnet["cv.glmnet"];
            Rcpp::Function predict_glmnet = glmnet["predict.glmnet"];
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
            for (int i = 1; i < glmnet_i.n_elem; i++) {
                par(glmnet_i(i) - 1) = glmnet_x(i);
            }
            return Rcpp::List::create(Rcpp::Named("par") = par,
                                      Rcpp::Named("value") = R_NilValue);
        } else {
            Rcpp::Environment glmnet("package:glmnet");
            Rcpp::Function glmnet_fit = glmnet["glmnet"];
            Rcpp::List out = glmnet_fit(
                data.cols(1, data.n_cols - 1),
                data.col(0),
                Rcpp::Named("family") = "gaussian",
                Rcpp::Named("lambda") = lambda
            );
            Rcpp::Function glmnet_predict = glmnet["predict.glmnet"];
            arma::vec fitted_values = Rcpp::as<arma::vec>(glmnet_predict(out, data.cols(1, data.n_cols - 1), Rcpp::Named("s") = lambda));
            double value = out["deviance"];
            return Rcpp::List::create(Rcpp::Named("par") = Rcpp::as<arma::vec>(out["beta"]).col(0),
                                      Rcpp::Named("value") = value / 2,
                                      Rcpp::Named("residuals") = data.col(0) - fitted_values);
        }
    } else if (theta.isNull()) {
        // Estimate theta in binomial/poisson/gaussian family
        arma::mat x = data.cols(1, data.n_cols - 1);
        arma::vec y = data.col(0);
        Rcpp::Environment fastglm("package:fastglm");
        Rcpp::Function fastglm_fit = fastglm["fastglm"];
        Rcpp::List out = fastglm_fit(x, y, family);
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
