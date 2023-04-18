#include "fastcpd.h"

//' Vanilla PELT implementation.
//' This function is not meant to be called directly by the user.
//'
//' @param data A data frame containing the data to be segmented.
//' @param beta Initial cost value.
//' @param segment_count Number of segments for initial guess.
//' @param trim Trimming for the boundary change points.
//' @param momentum_coef Momentum coefficient to be applied to each update.
//' @param k Function on number of epochs in SGD.
//' @param family Family of the model. Can be "binomial", "poisson", "lasso" or
//'   "gaussian". If not provided, the user must specify the cost function and
//'   its gradient (and Hessian).
//' @param epsilon Epsilon to avoid numerical issues. Only used for binomial and
//'   poisson.
//' @param min_prob Minimum probability to avoid numerical issues. Only used for
//'   poisson.
//' @param winsorise_minval Minimum value to be winsorised. Only used for
//'   poisson.
//' @param winsorise_maxval Maximum value to be winsorised. Only used for
//'   poisson.
//' @param p Number of parameters to be estimated.
//' @param cost Cost function to be used. If not specified, the default is
//'   the negative log-likelihood for the corresponding family.
//' @param cp_only Whether to return only the change points or with the cost
//'   values for each segment.
//'
//' @return Change points and corresponding cost values.
// [[Rcpp::export]]
Rcpp::List fastcpd_vanilla(
    arma::mat data,
    double beta,
    double segment_count,
    double trim,
    double momentum_coef,
    Rcpp::Function k,
    std::string family,
    double epsilon,
    double min_prob,
    double winsorise_minval,
    double winsorise_maxval,
    double p,
    Rcpp::Function cost,
    bool cp_only
) {

    int n = data.n_rows;
    std::vector<int> r_t_set = {0, 1};
    std::vector<std::vector<int>> cp_set(n + 1);
    std::fill(cp_set.begin(), cp_set.end(), std::vector<int>(1));
    cp_set[0] = {};
    std::vector<double> f_t(n + 1);
    f_t[0] = -beta;
    std::fill(f_t.begin() + 1, f_t.end(), 0);

    for (int t = 2; t <= n; t++) {
        int r_t_count = r_t_set.size();
        // number of cost values is the same as number of elemnts in R_t
        arma::vec cval(r_t_count);
        cval.fill(0);

        // for tau in R_t\{t-1}
        for (int i = 0; i < r_t_count - 1; i++) {
            int tau = r_t_set[i];

            if (t - tau >= 1) {
                cval[i] = Rcpp::as<double>(cost(data.rows(tau, t - 1)));
            }
        }

        cval[r_t_count - 1] = 0;
        arma::vec f_subset(r_t_set.size());
        for (unsigned int i = 0; i < r_t_set.size(); i++) {
            f_subset[i] = f_t[r_t_set[i]];
        }

        arma::vec obj = cval + f_subset + beta;
        double min_val = obj.min();
        int tau_star = r_t_set[obj.index_min()];

        std::vector<int> tau_star_cp_set = cp_set[tau_star + 1];
        cp_set[t] = tau_star_cp_set;
        cp_set[t].push_back(tau_star);
        // for (int i = 0; i < cp_set[t].size(); i++) {
        //     std::cout << cp_set[t][i] << " ";
        // }
        // std::cout << std::endl;
        // std::cout << "cval: " << cval << std::endl;
        // std::cout << "f_subset: " << f_subset << std::endl;
        // std::cout << "obj: " << obj << std::endl;
        // std::cout << "min_val: " << min_val << std::endl;
        // std::cout << "pruned_left: " << ((cval + f_subset) <= min_val) << std::endl;

        arma::uvec pruned_left = (cval + f_subset) <= min_val;
        arma::vec new_r_t_set = arma::vec(pruned_left.n_elem + 1);
        for (unsigned int i = 0; i < pruned_left.n_elem; i++) {
            new_r_t_set[i] = r_t_set[pruned_left[i]];
        }
        new_r_t_set[pruned_left.n_elem] = t;
        r_t_set = arma::conv_to<std::vector<int> >::from(new_r_t_set);

        f_t[t] = min_val;
    }

    std::vector<int> cp_set_final = cp_set[n];

    arma::uvec trim_cond1 = arma::find(arma::conv_to<arma::vec>::from(cp_set_final) >= trim * n);
    arma::uvec trim_cond2 = arma::find(arma::conv_to<arma::vec>::from(cp_set_final) <= (1 - trim) * n);

    std::vector<int> cp_set_trimmed = arma::conv_to<std::vector<int> >::from(arma::intersect(trim_cond1, trim_cond2));
    cp_set_trimmed.insert(cp_set_trimmed.begin(), 0);
    std::sort(cp_set_trimmed.begin(), cp_set_trimmed.end());

    std::vector<int> segment_indices;
    for (unsigned int i = 1; i < cp_set_trimmed.size(); i++) {
        if (cp_set_trimmed[i] - cp_set_trimmed[i - 1] < trim * n) {
            segment_indices.push_back(i - 1);
        }
    }

    std::vector<int> cp_set_trimmed_i = cp_set_trimmed, cp_set_trimmed_i1 = cp_set_trimmed;
    if (segment_indices.size() > 0) {
        for (unsigned int i = 0; i < segment_indices.size(); i++) {
            cp_set_trimmed_i.erase(cp_set_trimmed_i.begin() + segment_indices[i]) + 1;
            cp_set_trimmed_i1.erase(cp_set_trimmed_i.begin() + segment_indices[i]);
        }
        for (unsigned int i = 0; i < cp_set_trimmed_i.size(); i++) {
            cp_set_trimmed_i[i] = floor((cp_set_trimmed_i[i] + cp_set_trimmed_i1[i]) / 2);
        }
    }
    std::vector<int> cp_set_final_final = std::vector<int>(cp_set_trimmed_i.begin() + 1, cp_set_trimmed_i.end());

    if (cp_only) {
        return Rcpp::List::create(Rcpp::Named("cps") = cp_set);
    } else {
        cp_set_final_final.insert(cp_set_final_final.begin(), 0);
        cp_set_final_final.push_back(n);
        std::vector<double> cost_values;

        for (unsigned int i = 0; i < cp_set_final_final.size() - 1; i++) {
            double cost_value = Rcpp::as<double>(cost(data.rows(cp_set_final_final[i], cp_set_final_final[i + 1] - 1)));
            cost_values.push_back(cost_value);
        }
        return Rcpp::List::create(Rcpp::Named("cps") = cp_set_final_final,
                                  Rcpp::Named("cost") = cost_values);
    }
}
