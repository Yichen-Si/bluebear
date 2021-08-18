#include "cthmm.h"

void cthmm::forward() {
    alpha.col(0) = init_state_prob * Emtx.col(obs[0]);
    alpha.col(0) = alpha.col(0) / alpha.col(0).sum();
    for (int32_t j = 1; j < n_obs; ++j) {
        float d = distance[j] * 1.;
        ArrayXd p_stay = exp(-theta * d);
        for (int32_t i = 0; i < n_state; ++i) {
            alpha(i,j) = (p_stay(i) * alpha(i,j-1) + (1. - p_stay(i)) * (alpha.col(j-1) * Amtx.col(i)).sum()) * Emtx(i,obs[j]);
        }
        alpha.col(j) /= alpha.col(j).sum();
    }
};

void cthmm::backward() {
    beta.col(n_obs-1) = 1./n_state;
    for (int32_t j = n_obs-2; j >= 0; --j) {
        float d = distance[j+1] * 1.;
        ArrayXd p_stay = exp(-theta * d);
        for (int32_t i = 0; i < n_state; ++i) {
            beta(i,j) = (1. - p_stay(i)) * (beta.col(j+1) * Emtx.col(obs[j+1]) * Amtx.row(i).transpose()).sum();
            beta(i,j) += p_stay(i) * beta(i, j+1) * Emtx(i, obs[j+1]);
        }
        beta.col(j) /= beta.col(j).sum();
    }
};

void cthmm::EM(int32_t max_iter_EM, int32_t max_iter_NR, double tol) {
    int32_t n_iter = 0;
    ArrayXd update_size = ArrayXd::Zero(3) + 1.;
    double max_update = 1.;
    Eigen::IOFormat MtxFmt(3, Eigen::DontAlignCols, "\t", "\n");
    while ( n_iter < max_iter_EM && max_update > tol ) {
        // E-step
        printf("%d-th iteration, E step...\n", n_iter);
        forward();
        backward();
        marginal = alpha * beta;
        for (int32_t t = 0; t < n_obs; ++t) {
            marginal.col(t) /= marginal.col(t).sum();
        }
        // M-step
        printf("%d-th iteration, M step...\n", n_iter);
        // 1 update emisstion
        ArrayXXd Emtx_new = ArrayXXd::Zero(n_state, n_category);
        for (int32_t i = 0; i < n_state; ++i) {
            for (int32_t j = 0; j < n_obs; ++j) {
                Emtx_new(i,obs[j]) += marginal(i, j);
            }
            Emtx_new.row(i) /= Emtx_new.row(i).sum();
        }
std::cout << "Emtx\n" << Emtx_new.format(MtxFmt) << '\n';
        // 2.1 update transition - prob. cond. change
        init_state_prob = marginal.col(0);
        ArrayXXd Amtx_new = ArrayXXd::Zero(n_state, n_state);
        for (int32_t t = 1; t < n_obs; ++t) {
            Amtx_new += (alpha.col(t-1).matrix() * beta.col(t).matrix().transpose()).array()
                      * ((1.-exp(-theta*distance[t])).matrix() * Emtx.col(obs[t]).matrix().transpose()).array();
        }
        Amtx_new *= Amtx;
        for (int32_t i = 0; i < n_state; ++i) {
            Amtx_new.row(i) += Amtx_new.row(i).maxCoeff() * 1e-3;
            Amtx_new(i,i) = 0;
            Amtx_new.row(i) /= Amtx_new.row(i).sum();
        }
std::cout << "Amtx\n" << Amtx_new.format(MtxFmt) << '\n';

        // 2.2 update transition - jump rate
        ArrayXd theta_new = theta;
        ArrayXd theta_1st = ArrayXd::Zero(n_state);
        ArrayXd theta_2nd = ArrayXd::Zero(n_state);
        ArrayXd theta_delta = ArrayXd::Zero(n_state) + 1.;
        int32_t n_iter_NR = 0;
        while ( theta_delta.abs().maxCoeff() > tol && n_iter_NR < max_iter_NR ) {
            for (int32_t t = 1; t < n_obs; ++t) {
                ArrayXd p_stay_org = exp(-theta * distance[t]);
                ArrayXd p_stay_new = exp(-theta_new * distance[t]);
                for (int32_t i = 0; i < n_state; ++i) {
                    double tmp = 0;
                    theta_1st(i) -= distance[t] * alpha(i, t-1) * beta(i, t) * p_stay_org(i) * Emtx(i, obs[t]);
                    for (int32_t j = 0; j < n_state; ++j) {
                        if (j != i) {
                            tmp += Amtx(i,j) * beta(j, t) * Emtx(j,obs[t]);
                        }
                    }
                    theta_1st(i) += tmp * alpha(i, t-1) * (1.-p_stay_org(i)) * distance[t] * p_stay_new(i) / (1.-p_stay_new(i));
                    theta_2nd(i) += tmp * alpha(i, t-1) * (1.-p_stay_org(i)) * (-pow(distance[t], 2)) * p_stay_new(i) / pow(1.-p_stay_new(i), 2);
                }
            }
            theta_delta = - theta_1st / theta_2nd;
            theta_new += theta_delta;
            for (int32_t i = 0; i < n_state; ++i) {
                if (theta_new(i) < 1./(1e7*dist_scale)) {
                    theta_new(i) = (theta_new(i) - theta_delta(i)) * 1.5;
                }
            }
            n_iter_NR++;
        }
        if (n_iter_NR >= max_iter_NR && theta_delta.maxCoeff() > tol) {
            printf("Max iteration reached in N-R in the %d-th EM iteration, last rate change is %.3e.\n", n_iter, theta_delta.maxCoeff());
        } else {
            printf("%d iterations of N-R are used in the %d-th EM iteration, last rate change is %.3e.\n", n_iter_NR, n_iter, theta_delta.maxCoeff());
        }
std::cout << "Theta\n" << theta_new.format(MtxFmt) << '\n';
        update_size(0) = (theta_new - theta).abs().maxCoeff();
        update_size(1) = (Emtx_new - Emtx).abs().maxCoeff();
        update_size(2) = (Amtx_new - Amtx).abs().maxCoeff();
        max_update = update_size.maxCoeff();
        Emtx = Emtx_new;
        Amtx = Amtx_new;
        theta = theta_new;
        n_iter++;
    }
};

void cthmm::leave_one_out() {
    loo.resize(n_state, n_obs);
    for (int32_t t = 1; t < n_obs-1; ++t) {
        ArrayXd a(n_state);
        ArrayXd b(n_state);
        ArrayXd p_stay = exp(-theta * distance[t]);
        ArrayXd p_next = exp(-theta * distance[t+1]);
        for (int32_t i = 0; i < n_state; ++i) {
            a(i) = (alpha.col(t-1) * (1. - p_stay(i)) * Amtx.col(i)).sum()
                  + alpha(i, t-1) * p_stay(i);
            b(i) = (beta.col(t+1) * (1. - p_next(i)) * Amtx.row(i).transpose() * Emtx.col(obs[t+1])).sum()
                  + beta(i, t+1) * p_next(i) * Emtx(i, obs[t+1]);
        }
        loo.col(t) = (a * b) / (a * b).sum();
    }
    for (int32_t i = 0; i < n_state; ++i) {
        int32_t t = 0;
        ArrayXd p_next = exp(-theta * distance[t+1]);
        loo(i, 0) = (beta.col(t+1) * (1. - p_next(i)) * Amtx.row(i).transpose() * Emtx.col(obs[t+1])).sum()
                   + beta(i, t+1) * p_next(i) * Emtx(i, obs[t+1]);
        t = n_obs-1;
        ArrayXd p_stay  = exp(-theta * distance[t]);
        loo(i, t) = (alpha.col(t-1) * (1. - p_stay(i)) * Amtx.col(i)).sum()
                   + alpha(i, t-1) * p_stay(i);
    }
    loo.col(0) /= loo.col(0).sum();
    loo.col(n_obs-1) /= loo.col(n_obs-1).sum();
};

void cthmm::viterbi() {
    phi.resize(n_state, n_obs);
    max_ll = init_state_prob.log() + Emtx.col(obs[0]).log();
    ArrayXd delta = max_ll;
    for (int32_t t = 1; t < n_obs; ++t) {
        ArrayXd p_leave = (1.-exp(-theta * distance[t])).log();
        ArrayXd tmp = delta - theta * distance[t];
        for (int32_t i = 0; i < n_state; ++i) {
            phi(i, t) = i;
            for (int32_t j = 0; j < n_state; ++j) {
                if (j != i && delta(j) + p_leave(j) + log(Amtx(j,i)) > tmp(i) ) {
                    phi(i, t) = j;
                }
            }
            max_ll(i) = (phi(i,t) == i) ? tmp(i) : delta(phi(i, t)) + p_leave(phi(i, t)) + log(Amtx(phi(i, t),i) );
        }
        max_ll += Emtx.col(obs[t]).log();
        delta = max_ll;
    }
    // Backtrack
    viterbi_path.resize(n_obs);
    viterbi_path[n_obs-1] = 0;
    for (int32_t i = 1; i < n_state; ++i) {
        if (max_ll(i) > max_ll(viterbi_path[n_obs-1])) {
            viterbi_path[n_obs-1] = i;
        }
    }
    for (int32_t t = n_obs-1; t > 0; --t) {
        viterbi_path[t-1] = phi(viterbi_path[t], t);
    }
}
