#include "cramore.h"
#include "tsv_reader.h"
#include "hts_utils.h"
#include "utils.h"
#include "cthmm.h"
#include <iomanip>
#include <fstream>
#include <iostream>
#include <sstream>

#include <Eigen/Dense>
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::ArrayXXd;

// Goal: HMM for methylation status of CpG sites
int32_t cpgHMM(int32_t argc, char** argv) {
    std::string inTsv, initEmis, initTrans, reg, out, chrom;
    int32_t max_iter_outer = 20, max_iter_inner = 20;
    double  tol_outer = 1e-6, tol_inner = 1e-6;
    int32_t start, end, st, ed, n_sample;
    int32_t n_state, n_obs;
    double  chunk_size_double = 1e6;
    int32_t min_obs = (int) 1e4;
    int32_t ac_column = 6, pos_column = 2;
    bool    output_full_likelihood = 1, output_viterbi = 1, output_leave_one_out = 1;
    double  dist_scale = 1e-3;
    bool    optim_EM = 1;

  paramList pl;
    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input", NULL)
    LONG_STRING_PARAM("in-tsv",&inTsv, "Input file containing CpG position and allele count")
    LONG_INT_PARAM("position_col",&pos_column,"Which column contains CpG position")
    LONG_INT_PARAM("ac_col",&ac_column,"Which column contains sample allele count")
    LONG_STRING_PARAM("region",&reg,"Region to process (index required)")
    LONG_INT_PARAM("sample-size",&n_sample,"Total sample size the input AC is based on")
    LONG_STRING_PARAM("init-emission",&initEmis, "Input file containing initial emission probabilities")
    LONG_STRING_PARAM("init-transition",&initTrans,"Input file containing initial transition scale and transition probabilities for each state")

    LONG_PARAM_GROUP("Additional Options", NULL)
    LONG_INT_PARAM("max-iter-outer",&max_iter_outer,"Maximun iterations of EM")
    LONG_INT_PARAM("max-iter-inner",&max_iter_inner,"Maximun iterations of Newton's method inside each EM iteration")
    LONG_DOUBLE_PARAM("tol-outer",&tol_outer,"Tolerance to declare convergence in EM")
    LONG_DOUBLE_PARAM("tol-inner",&tol_inner,"Tolerance to declare convergence in nuemerical optimization transition rate (evaluated in scale, unit kb)")
    LONG_DOUBLE_PARAM("min-obs",&min_obs,"Minimum observation in a block")
    LONG_DOUBLE_PARAM("chunk-size",&chunk_size_double,"Maximun window to process at one time")
    LONG_INT_PARAM("optim-EM",&optim_EM,"If use EM as parameter learning method")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output file prefix")
    LONG_INT_PARAM("output-likelihood",&output_full_likelihood, "Whether to output full conditional likelihood matrix")
    LONG_INT_PARAM("output-viterbi",&output_viterbi, "Whether to output viterbi path")
    LONG_INT_PARAM("output-loo",&output_leave_one_out, "Whether to output leave-one-out prediction")

  END_LONG_PARAMS();
  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

if ( inTsv.empty() || reg.empty() || initTrans.empty() || initEmis.empty() || out.empty() ) {
    error("[E:%s:%d %s] --in-tsv, --region, --init-emission, --init-transition, --out are required parameters",__FILE__,__LINE__,__FUNCTION__);
}
Eigen::IOFormat MtxTsvFmt(4, Eigen::DontAlignCols, "\t", "\n");

std::string line;
std::vector<std::string> v;
split(v, ":-", reg);
chrom = v[0];
if (v.size() > 2) {
    if (!str2int32(v[1], start) || !str2int32(v[2], end)) {
      error("Invalid region.");
    }
} else {
    start = 0;
    end = (int) (300e6);
}

// Read initial parameters
// Transition
n_state = 0;
std::vector<double> trans_scale;
std::vector<std::vector<double> > Amtx_v;
std::ifstream rf;
rf.open(initTrans, std::ifstream::in);
while (std::getline(rf, line)) {
    if (line.at(0) != '#') {
        split(v, "\t", line);
        try {
            trans_scale.push_back(std::stod(v[1]));
            Amtx_v.push_back( std::vector<double>(v.size()-2) );
            for (uint32_t i = 2; i < v.size(); ++i) {
                Amtx_v[n_state][i-2] = std::stod(v[i]);
            }
            n_state++;
        } catch(...) {
            error("Invalid initial transition parameter.");
        }
    }
}
rf.close();
ArrayXd init_prob(n_state);
ArrayXXd Amtx(n_state, n_state);
printf ("Read %d states from initial transition parameters.\n", n_state);
for (int32_t i = 0; i < n_state; ++i) {
    init_prob(i) = trans_scale[i];
    if ((int32_t) Amtx_v[i].size() != n_state) {
        error("Invalid initial transition parameter.");
    }
    printf("%.1f\t", trans_scale[i]);
    for (int32_t j = 0; j < n_state; ++j) {
        Amtx(i,j) = Amtx_v[i][j];
    }
    Amtx(i,i) = 0;
    Amtx.row(i) /= Amtx.row(i).sum();
    trans_scale[i] *= dist_scale;
}
init_prob /= init_prob.sum();
std::cout << '\n' << Amtx.format(MtxTsvFmt) << '\n';

// Emission
rf.open(initEmis, std::ifstream::in);
std::vector<std::vector<double> > Emtx_v;
std::vector<int32_t> ac_cut;
while (std::getline(rf, line)) {
    split(v, "\t", line);
    if (line.at(0) == '#') {
        for (uint32_t i = 1; i < v.size(); ++i) {
            ac_cut.push_back(std::stoi(v[i]));
        }
    } else {
        std::vector<double> ev;
        double dsum = 0.0;
        for (uint32_t i = 1; i < v.size(); ++i) {
            ev.push_back(std::stod(v[i]));
            dsum += ev[i-1];
        }
        for (uint32_t i = 1; i < v.size(); ++i) {
            ev[i-1] /= dsum;
        }
        Emtx_v.push_back(ev);
    }
}
rf.close();
if ((int32_t) Emtx_v.size() != n_state || ac_cut[0] != 0) {
    error("Incompatible initial transition and emission.");
}
int32_t n_category = (int32_t) ac_cut.size();
ArrayXXd Emtx(n_state, n_category);
for (int32_t i = 0; i < n_state; ++i) {
    for (int32_t j = 0; j < n_category; ++j) {
        Emtx(i,j) = Emtx_v[i][j];
    }
    Emtx.row(i) /= Emtx.row(i).sum();
}
printf ("Read %ld states and %ld observation categories from initial emission parameters.\n", Emtx.rows(), Emtx.cols() );
std::cout << Emtx.format(MtxTsvFmt) << '\n';

// Read input observaitons
int32_t  chunk_size = (int) chunk_size_double;
std::vector<int32_t> position;
std::vector<float> distance;
std::vector<int32_t> obs;
st = start;
ed = end;
if (end - start > chunk_size) {
    ed = st+chunk_size-1;
}
printf("%s:%d-%d\n", chrom.c_str(), st, ed);
tsv_reader tr = (inTsv.c_str());
tr.jump_to(chrom.c_str(), st, ed);
std::map<int32_t, int32_t > ac_cut_ct;
for (uint32_t i = 0; i < ac_cut.size(); ++i) {
    ac_cut_ct[i] = 0;
}
while(tr.read_line(line) > 0) {
	position.push_back(tr.int_field_at(pos_column) );
    int32_t ac  = tr.int_field_at(ac_column);
    if (ac > n_sample) {
        ac = 2*n_sample - ac;
    }
    int32_t indx = binarySearch(ac_cut, 0, n_category-1, ac);
    ac_cut_ct[indx]++;
    obs.push_back( indx );
}
n_obs = obs.size();
if (n_obs < min_obs) {
    error("Region does not contain enough CpG sites.");
}
printf( "Read %d observations, count by categories:\n", n_obs );
for (uint32_t i = 0; i < ac_cut.size(); ++i) {
    std::cout << ac_cut[i] << '\t' << ac_cut_ct[i] << '\n';
}
distance.resize(n_obs);
for (int32_t j = 1; j < n_obs; ++j) {
    distance[j] = (position[j] - position[j-1])*dist_scale;
}
distance[0] = 0;
if (*(std::min_element(distance.begin(), distance.end())) < 0) {
    error("Input should be sorted with non-decreasing positions.");
}

// Initialize HMM object
notice("Initializing HMM object");
cthmm hmm_obj(obs, distance, dist_scale, trans_scale, Amtx, Emtx, init_prob);
// EM
notice("Start parameter learning");
if (optim_EM) {
    hmm_obj.EM(max_iter_outer, max_iter_inner, tol_outer, tol_inner);
} else {
    hmm_obj.mixed_optim(max_iter_outer, max_iter_inner, tol_outer, tol_inner);
}
// Final
notice("Finish parameter learning, start final LOO");
hmm_obj.leave_one_out();
notice("Finish LOO, start Viterbi");
hmm_obj.viterbi();
notice("Finish Viterbi");

// Output
std::ofstream mf;
std::string outf;
// Output parameters
outf = out + ".transition.tsv";
mf.open(outf.c_str(), std::ofstream::out);
mf << "#State\tScale";
for (int32_t i = 0; i < n_state; ++i) {
    mf << '\t' << i;
}
mf << '\n';
mf << std::setprecision(4) << std::fixed;
for (int32_t i = 0; i < n_state; ++i) {
    mf << std::to_string(i) << '\t' << 1./hmm_obj.theta(i)/dist_scale;
    for (int32_t j = 0; j < n_state; ++j) {
        mf << '\t' << hmm_obj.Amtx(i, j);
    }
    mf << '\n';
}
mf.close();
outf = out + ".emission.tsv";
mf.open(outf.c_str(), std::ofstream::out);
mf << "#State";
for (auto & v : ac_cut) {
    mf << '\t' << v;
}
mf << '\n';
for (int32_t i = 0; i < n_state; ++i) {
    mf << i << '\t' << hmm_obj.Emtx.row(i).format(MtxTsvFmt) << '\n';
}
mf.close();

// Output Viterbi path. chr,st,ed,state,n_obs,collapsed SFS
if (output_viterbi) {
    FILE *wf;
    outf = out + ".viterbi";
    wf = fopen(outf.c_str(), "w");
    int32_t v_st = position[0] - 1;
    int32_t n_tot = 0;
    std::vector<int32_t> sfs_window(n_category, 0);
    int8_t  v_state = hmm_obj.viterbi_path[0];
    for (int32_t i = 1; i < n_obs; ++i) {
        if (hmm_obj.viterbi_path[i] != v_state) {
            std::stringstream ss;
            ss << sfs_window[0];
            for (int32_t j = 1; j < n_category; ++j) {
                ss << ',' << sfs_window[j];
            }
            fprintf(wf, "%s\t%d\t%d\t%d\t%d\t%s\n", chrom.c_str(), v_st, position[i-1], n_tot, v_state,ss.str().c_str());
            v_state = hmm_obj.viterbi_path[i];
            v_st = position[i] - 1;
            n_tot = 0;
            std::fill(sfs_window.begin(), sfs_window.end(), 0);
        }
        n_tot++;
        sfs_window[obs[i]]++;
    }
    fclose(wf);
}

// Output conditional probabilities
if (output_full_likelihood) {
    outf = out + ".likelihood";
    mf.open(outf.c_str(), std::ofstream::out);
    for (int32_t i = 0; i < n_obs; ++i) {
        mf << position[i] << '\t' << hmm_obj.marginal.col(i).transpose().format(MtxTsvFmt) << '\n';
    }
    mf.close();
}

// Output leave-one-out probabilities
if (output_leave_one_out) {
    outf = out + ".loo";
    mf.open(outf.c_str(), std::ofstream::out);
    for (int32_t i = 0; i < n_obs; ++i) {
        mf << position[i] << '\t' << hmm_obj.loo.col(i).transpose().format(MtxTsvFmt) << '\n';
    }
    mf.close();
}


return 0;
}
