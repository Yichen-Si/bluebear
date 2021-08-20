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

void cthmm::update_matrix() {
    marginal = alpha * beta;
    for (int32_t t = 0; t < n_obs; ++t) {
        marginal.col(t) /= marginal.col(t).sum();
    }
    // 1 update emisstion
    ArrayXXd Emtx_new = ArrayXXd::Zero(n_state, n_category);
    for (int32_t i = 0; i < n_state; ++i) {
        for (int32_t j = 0; j < n_obs; ++j) {
            Emtx_new(i,obs[j]) += marginal(i, j);
        }
        Emtx_new.row(i) += Emtx_new.row(i).maxCoeff() * 0.01;
        Emtx_new.row(i) /= Emtx_new.row(i).sum();
    }
    // 2 update transition - prob. cond. change
    init_state_prob = marginal.col(0);
    ArrayXXd Amtx_new = ArrayXXd::Zero(n_state, n_state);
    for (int32_t t = 1; t < n_obs; ++t) {
        Amtx_new += (alpha.col(t-1).matrix() * beta.col(t).matrix().transpose()).array() * ((1.-exp(-theta*distance[t])).matrix() * Emtx.col(obs[t]).matrix().transpose()).array();
    }
    Amtx_new *= Amtx;
    for (int32_t i = 0; i < n_state; ++i) {
        Amtx_new.row(i) += Amtx_new.row(i).maxCoeff() * 1e-2;
        Amtx_new(i,i) = 0;
        Amtx_new.row(i) /= Amtx_new.row(i).sum();
    }
    update_size(1) = (Emtx_new - Emtx).abs().maxCoeff();
    update_size(2) = (Amtx_new - Amtx).abs().maxCoeff();
    Emtx = Emtx_new;
    Amtx = Amtx_new;
};

void cthmm::mixed_optim(int32_t max_iter, int32_t max_iter_inner, double tol, double tol_inner, int32_t criterion) {
    update_size = ArrayXd::Zero(3) + 1.;
    int32_t n_iter = 0;
    double max_update = 1.;
    Eigen::IOFormat MtxFmt(3, Eigen::DontAlignCols, "\t", "\n");
    forward();
    backward();
    update_matrix();
    double ll_old = 0, ll_new = 0;
    if (criterion == 1) {
        leave_one_out();
        leave_one_out_composite_map();
        ll_old = loo_map;
    } else if (criterion == 2) {
        viterbi();
        ll_old = viterbi_ll;
    } else {
        leave_one_out();
        leave_one_out_composite_posterior();
        ll_old = loo_post;
    }
    while ( n_iter < max_iter && max_update > tol ) {
        printf("%d-th iteration.\n", n_iter);
        // transition rate = arg max LOO composite likelihood
        ArrayXd theta_org = theta;
        ArrayXd theta_new = theta;
        for (int32_t state_idx = 0; state_idx < n_state; ++state_idx) {
            double x1, x2;
            x1 = 1./(5e6 * dist_scale);
            x2 = 1./(1e1 * dist_scale);
            double v = 1., v_prev = 1.;
            BrentObj brent(x1,x2);
            double arg = brent.arg;
            int32_t n_iter_inner = 0;
            theta = theta_org; // Update each rate separately
            while (brent.status > 0 && n_iter_inner < max_iter_inner) {
                theta(state_idx) = arg;
                forward(); backward();
                if (criterion == 1) {
                    leave_one_out();
                    leave_one_out_composite_map();
                    v = -loo_map;
                } else if (criterion == 2) {
                    viterbi();
                    v = -viterbi_ll;
                } else {
                    leave_one_out();
                    leave_one_out_composite_posterior();
                    v = -loo_post;
                }
                arg = brent.local_min_rc(x1,x2,v);
                n_iter_inner++;
                if (n_iter_inner > 1 && abs(v-v_prev) < tol_inner) {
                    break;
                }
                v_prev = v;
            } // Finish one transition rate
            theta_new(state_idx) = theta(state_idx);
            printf("Theta_%d - Output after %d iterations: a=%.3f, b=%.3f, arg=%.3f, value=%.3f\n", state_idx, n_iter_inner, 1./x1, 1./x2, 1./arg, v);
        } // Finish separate updating each transition rate
        theta = theta_new;
        forward(); backward();
        viterbi();
        leave_one_out();
        leave_one_out_composite_map();
        leave_one_out_composite_posterior();
        if (criterion == 1) {
            ll_new = loo_map;
        } else if (criterion == 2) {
            ll_new = viterbi_ll;
        } else {
            ll_new = loo_post;
        }
        History* rec = new History();
        rec->iter = n_iter;
        rec->failed = (ll_new > ll_old) ? 0 : 1;
        rec->loo_map = loo_map;
        rec->loo_post = loo_post;
        rec->viterbi_ll = viterbi_ll;
        rec->theta = theta_new;
        track.push_back(rec);
        if (ll_new > ll_old) {
            update_matrix();
            printf("Ascent succeed in the %d-th iteration - ll_new=%.3f, ll_old=%.3f, difference=%.2f\n", n_iter, ll_new, ll_old, ll_new-ll_old);
            std::cout << "Emtx\n" << Emtx.format(MtxFmt) << '\n';
            std::cout << "Amtx\n" << Amtx.format(MtxFmt) << '\n';
            std::cout << "Scale(bp):\t" << (theta.inverse().transpose()/dist_scale).format(MtxFmt) << '\n';
            update_size(0) = (theta - theta_org).abs().maxCoeff();
            max_update = update_size.maxCoeff();
            ll_old = ll_new;
        } else {
            theta = theta_org; // failed, roll back to previous state
            forward();
            backward();
            update_matrix();
            max_update = 1.;
        }
        n_iter++;
    }
};


double cthmm::min_obj_theta_indivisual(int32_t state_idx, double x) {
    // Compute part of -Q corresponding to \theta_i
    double res = 0;
    int32_t i = state_idx;
    for (int32_t t = 1; t < n_obs; ++t) {
        ArrayXXd smtx = ArrayXXd::Zero(n_state, n_state);
        for (int32_t ii = 0; ii < n_state; ii++) {
            smtx(ii,ii) = alpha(ii,t-1) * exp(-theta(ii) * distance[t]) * beta(ii,t) * Emtx(ii, obs[t]);
            for (int32_t jj = 0; jj < n_state; jj++) {
                if (jj==ii) {continue;}
                smtx(ii,jj) = alpha(ii,t-1) * (1.-exp(-theta(ii) * distance[t])) * Amtx(ii,jj) * beta(jj,t) * Emtx(jj, obs[t]);
            }
        }
        smtx /= smtx.sum();
        res += smtx(i,i) * (-x * distance[t]);
        smtx(i,i) = 0;
        res += smtx.row(i).sum() * log(1.-exp(-x * distance[t]));
    }
    return -res;
};

double cthmm::optim_brent_theta_individual(int32_t state_idx, int32_t max_iter, double tol) {
    // Optimize a single transition rate
    double x1, x2;
    // Initial interval to search for local min
    x1 = 1./(2e6 * dist_scale);
    x2 = 1./(1e1 * dist_scale);
    double v = 1., v_prev = 1.;
    BrentObj brent(x1,x2);
    double arg = brent.arg;
    int32_t n_iter = 0;
    while (brent.status > 0 && n_iter < max_iter) {
        v = min_obj_theta_indivisual(state_idx, arg);
        arg = brent.local_min_rc(x1,x2,v);
        n_iter++;
        if (n_iter > 1 && abs(v-v_prev) < tol) {
            break;
        }
        v_prev = v;
    }
    printf("Theta_%d - Output after %d iterations: a=%.3f, b=%.3f, arg=%.3f, value=%.3f\n", state_idx, n_iter, 1./x1/dist_scale, 1./x2/dist_scale, arg, v);
    return arg;
};

void cthmm::min_obj_theta(ArrayXd x, ArrayXd& res) {
    // Compute part of -Q corresponding to \theta
    res = ArrayXd::Zero(n_state);
    for (int32_t t = 1; t < n_obs; ++t) {
        ArrayXXd smtx = ArrayXXd::Zero(n_state, n_state);
        for (int32_t i = 0; i < n_state; i++) {
            smtx(i,i) = alpha(i,t-1) * exp(-theta(i) * distance[t]) * beta(i,t) * Emtx(i, obs[t]);
            for (int32_t j = 0; j < n_state; j++) {
                if (j==i) {continue;}
                smtx(i,j) = alpha(i,t-1) * (1.-exp(-theta(i) * distance[t])) * Amtx(i,j) * beta(j,t) * Emtx(j, obs[t]);
            }
        }
        smtx /= smtx.sum();
        for (int32_t i = 0; i < n_state; ++i) {
            res(i) += smtx(i,i) * (-x(i) * distance[t]);
            smtx(i,i) = 0;
            res(i) += smtx.row(i).sum() * log(1.-exp(-x(i) * distance[t]));
        }
    }
    res *= (-1.);
};

void cthmm::optim_brent_theta(int32_t max_iter, double tol, ArrayXd& arg) {
    // Optimize a single transition rate
    ArrayXd x1 = ArrayXd::Zero(n_state) + 1./(5e6*dist_scale);
    ArrayXd x2 = ArrayXd::Zero(n_state) + 1./(1e2*dist_scale);
    ArrayXd v = ArrayXd::Zero(n_state) + 1.;
    ArrayXd v_prev = v;
    std::vector<BrentObj*> brent_v;
    Eigen::IOFormat VecFmt(5, Eigen::DontAlignCols, "\t", "\t");
    for (int32_t i = 0; i < n_state; ++i) {
        BrentObj* bt = new BrentObj(x1(i), x2(i));
        arg(i) = bt->arg;
        brent_v.push_back(bt);
    }
    int32_t n_iter = 0, status = 1;
    while (status > 0 && n_iter < max_iter) {
        min_obj_theta(arg, v);
        status = 0;
        for (int32_t i = 0; i < n_state; ++i) {
            arg(i) = brent_v[i]->local_min_rc(x1(i),x2(i),v(i));
            status += brent_v[i]->status;
        }
        n_iter++;
        if (n_iter > 1 && (v-v_prev).abs().maxCoeff() < tol) {
            break;
        }
        v_prev = v;
    }
    printf("Partial Q value for theta after %d iterations: ", n_iter);
    std::cout << (-v).format(VecFmt) << '\n';
    for (auto bt : brent_v) {
        delete bt;
    }
};

void cthmm::EM(int32_t max_iter_EM, int32_t max_iter_inner, double tol_EM, double tol_inner, int32_t optim_method) {
    update_size = ArrayXd::Zero(3) + 1.;
    int32_t n_iter = 0;
    double max_update = 1.;
    Eigen::IOFormat MtxFmt(3, Eigen::DontAlignCols, "\t", "\n");
    while ( n_iter < max_iter_EM && max_update > tol_EM ) {
        // E-step
        printf("%d-th iteration, E step...\n", n_iter);
        forward();
        backward();
        // M-step
        printf("%d-th iteration, M step...\n", n_iter);
        // 2.2 update transition rate
        ArrayXd theta_new = theta;
        if (optim_method==0) { // Brent
            optim_brent_theta(max_iter_inner,tol_inner,theta_new);
        } else { // Newton
            NewtonRaphson(max_iter_inner, tol_inner, theta_new);
        }
        update_size(0) = (theta_new - theta).abs().maxCoeff();
        std::cout << "Scale(kb):\t" << (theta_new.inverse().transpose()).format(MtxFmt) << '\n';
        // Update emisstion & transition conditional probabilities
        update_matrix();
        max_update = update_size.maxCoeff();
        theta = theta_new;
        n_iter++;
        if (n_iter % record_freq == 0) {
            std::cout << "Emtx\n" << Emtx.format(MtxFmt) << '\n';
            std::cout << "Amtx\n" << Amtx.format(MtxFmt) << '\n';
            viterbi();
            leave_one_out();
            leave_one_out_composite_map();
            leave_one_out_composite_posterior();
            printf("\nLOO composite likelihood after %d iterations: %.2f\n", n_iter, loo_map);
            History* rec = new History();
            rec->iter = n_iter;
            rec->failed = 0;
            rec->loo_map = loo_map;
            rec->loo_post = loo_post;
            rec->viterbi_ll = viterbi_ll;
            rec->theta = theta;
            track.push_back(rec);
        }
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
        loo.col(t) = (a * b);
        loo.col(t)/= loo.col(t).sum();
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

void cthmm::leave_one_out_composite_posterior() {
    if (loo.rows() != n_state || loo.cols() != n_obs) {
        leave_one_out();
    }
    loo_post = 0;
    for (int32_t t = 1; t < n_obs-1; ++t) {
        loo_post += log((loo.col(t) * Emtx.col(obs[t])).sum());
    }
};

void cthmm::leave_one_out_composite_map() {
    if (loo.rows() != n_state || loo.cols() != n_obs) {
        leave_one_out();
    }
    loo_map = 0;
    for (int32_t t = 1; t < n_obs-1; ++t) {
        int32_t max_indx = 0;
        for (int32_t i = 1; i < n_state; ++i) {
            if (loo(i, t) > loo(max_indx, t)) {
                max_indx = i;
            }
        }
        loo_map += log( loo(max_indx, t) * Emtx(max_indx, obs[t]) );
    }
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
    viterbi_ll = max_ll.maxCoeff();
};


void cthmm::NewtonRaphson(int32_t max_iter_inner, double tol_inner, ArrayXd& theta_new) {
    ArrayXd theta_tmp = theta;
    ArrayXd theta_1st = ArrayXd::Zero(n_state);
    ArrayXd theta_2nd = ArrayXd::Zero(n_state);
    ArrayXd theta_delta = ArrayXd::Zero(n_state) + 1.;
    int32_t n_iter_inner = 0;
    ArrayXd theta_inv_delta = ArrayXd::Zero(n_state) + 1.;
    while ( theta_inv_delta.maxCoeff() > tol_inner && n_iter_inner < max_iter_inner ) {
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
        theta_tmp = theta_new + theta_delta;
        for (int32_t i = 0; i < n_state; ++i) {
            if (theta_tmp(i) < 1./(1e7*dist_scale)) {
                theta_tmp(i) = abs(theta_new(i)) * 1.5;
            }
            if (theta_tmp(i) > 1./(10*dist_scale)) {
                theta_tmp(i) = abs(theta_new(i)) * 0.8;
            }
            if ( std::isnan(theta_tmp(i)) ) {
                theta_tmp(i) = theta(i) * 1.2;
            }
        }
        theta_inv_delta = (theta_new.inverse() - theta_tmp.inverse()).abs();
        theta_new = theta_tmp;
        n_iter_inner++;
    }
    printf("%d iterations of N-R are used, last scale change is %.3e.\n", n_iter_inner, theta_inv_delta.maxCoeff());
};



double obj_map_wrap(const VectorXd& vals_inp, VectorXd* grad_out, void* opt_data)
{
std::cout << "Evaluation...\t";
    cthmm* hmm = (cthmm*) opt_data;
    hmm->theta = vals_inp.array().inverse();
    hmm->forward();
    hmm->backward();
    hmm->leave_one_out();
    hmm->leave_one_out_composite_map();
    double obj_val = hmm->loo_map;
std::cout << obj_val << "\n";
    return -obj_val;
}

void nm_loo_map_optim(cthmm* hmm_obj, int32_t max_iter, double tol) {
    // Use Nelder-Mead to optimize transition rate parameters
    int32_t n_iter = 0, n_iter_improved = 0;
    double max_update = 1.;
    Eigen::IOFormat MtxFmt(3, Eigen::DontAlignCols, "\t", "\n");
    hmm_obj->forward();
    hmm_obj->backward();
    hmm_obj->update_matrix();
    hmm_obj->leave_one_out();
    hmm_obj->leave_one_out_composite_map();
    double ll_old = hmm_obj->loo_map, ll_new = 0;
    double ll_older = ll_old;
    while ( n_iter < max_iter && max_update > tol ) {
        // transition rate = arg max LOO composite likelihood
        ArrayXd theta_org = hmm_obj->theta;
        VectorXd theta = (hmm_obj->theta).inverse().matrix();
        bool nm_success = optim::nm(theta, obj_map_wrap, hmm_obj);
        printf("%d-th iteration. Nelder-Mead for theta exit %d.\n", n_iter, nm_success);
        hmm_obj->viterbi();
        hmm_obj->leave_one_out_composite_posterior();
        ll_new = hmm_obj->loo_map;
        History* rec = new History();
        rec->iter = n_iter;
        rec->failed = (ll_new > ll_old) ? 0 : 1;
        rec->loo_map = hmm_obj->loo_map;
        rec->loo_post = hmm_obj->loo_post;
        rec->viterbi_ll = hmm_obj->viterbi_ll;
        rec->theta = hmm_obj->theta;
        hmm_obj->track.push_back(rec);
        if (ll_new > ll_old) {
            hmm_obj->update_matrix();
            printf("Ascent succeed in the %d-th iteration - ll_new=%.3f, ll_old=%.3f, difference=%.2f\n", n_iter, ll_new, ll_old, ll_new-ll_old);
            n_iter_improved++;
            std::cout << "Emtx\n" << (hmm_obj->Emtx).format(MtxFmt) << '\n';
            std::cout << "Amtx\n" << (hmm_obj->Amtx).format(MtxFmt) << '\n';
            std::cout << "Scale(kb):\t" << ((hmm_obj->theta).inverse().transpose()).format(MtxFmt) << '\n';
            hmm_obj->update_size(0) = (theta.array().inverse() - theta_org).abs().maxCoeff();
            max_update = (hmm_obj->update_size).maxCoeff();
            ll_old = ll_new;
        } else {
            printf("Ascent failed in the %d-th iteration - ll_new=%.3f, ll_old=%.3f, difference=%.2f\n", n_iter, ll_new, ll_old, ll_new-ll_old);
            hmm_obj->theta = theta_org; // failed, roll back to previous state
            hmm_obj->forward();
            hmm_obj->backward();
            hmm_obj->update_matrix();
            max_update = 1.;
        }
        n_iter++;
    }
    printf("Finish %d iterations, %d increase LOO MAP - ll_new=%.3f, ll_original=%.3f, difference=%.2f\n", n_iter, n_iter_improved, ll_new, ll_older, ll_new-ll_older);
};

// Box constrain
VectorXd constr_fn(const VectorXd& vals_inp, MatrixXd* jacob_out, void* constr_data) {
    std::cout << "Box...\n";
    Box* bound = (Box*) constr_data;
    VectorXd box(vals_inp.rows());
    for (int32_t i = 0; i < (int32_t) box.rows(); ++i) {
        box(i) = (vals_inp(i)-bound->lower(i))*(bound->upper(i)-vals_inp(i));
    }
    return box;
};
