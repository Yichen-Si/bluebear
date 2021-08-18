#ifndef CTHMM_H
#define CTHMM_H

#include <string>
#include <vector>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <map>
#include <cmath>
#include <cfloat>
#include <iostream>

#include <Eigen/Dense>
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::ArrayXXd;
using Eigen::ArrayXd;

class cthmm {
public:
    int32_t n_state, n_obs, n_category;
    std::vector<int32_t>& obs;
    std::vector<float>& distance;
    double dist_scale;
    ArrayXXd Amtx, Emtx;
    ArrayXXd alpha, beta, loo, marginal;
    ArrayXd theta, init_state_prob, max_ll;
    Eigen::Array<int8_t,Eigen::Dynamic,Eigen::Dynamic> phi;
    std::vector<int8_t> viterbi_path;

    cthmm(std::vector<int32_t>& _obs, std::vector<float>& _dist, double _s,
          std::vector<double>& _scale, ArrayXXd _A, ArrayXXd _E, ArrayXd _pi ) :
          obs(_obs), distance(_dist), dist_scale(_s), Amtx(_A), Emtx(_E), init_state_prob(_pi)
          {
              n_state    = _scale.size();
              n_obs      = obs.size();
              n_category = Emtx.cols();
              alpha = MatrixXd::Zero(n_state, n_obs);
              beta  = MatrixXd::Zero(n_state, n_obs);
              theta.resize(n_state);
              for (int32_t i = 0; i < n_state; ++i) {
                  theta(i) = 1./_scale[i];
              }
          }
    void forward();
    void backward();
    void EM( int32_t max_iter_EM = 20, int32_t max_iter_NR = 20, double tol = 1e-8 );
    void leave_one_out();
    void viterbi();
};




#endif
