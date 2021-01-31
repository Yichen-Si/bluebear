#include "cramore.h"
#include "bcf_filter_arg.h"
#include "bcf_ordered_reader.h"
// #include "bcf_ordered_writer.h"
#include "hts_utils.h"
#include <bits/stdc++.h>
// #include "ibs0.h"
#include <sstream>


int32_t DoubletonBlockInfo(int32_t argc, char** argv) {

  std::string inVcf, inQuery, reg, out;
  int32_t verbose = 10000;

  bcf_vfilter_arg vfilt;
  bcf_gfilter_arg gfilt;

  paramList pl;
  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input", NULL)
    LONG_STRING_PARAM("in-vcf",&inVcf, "Input VCF/BCF file")
    // LONG_STRING_PARAM("region",&reg,"Genomic region to focus on")
    LONG_STRING_PARAM("in-query",&inQuery, "Position and individual pairs to query")

    LONG_PARAM_GROUP("Variant Filtering Options", NULL)
    LONG_MULTI_STRING_PARAM("apply-filter",&vfilt.required_filters, "Require at least one of the listed FILTER strings")
    LONG_STRING_PARAM("include-expr",&vfilt.include_expr, "Include sites for which expression is true")
    LONG_STRING_PARAM("exclude-expr",&vfilt.exclude_expr, "Exclude sites for which expression is true")

    LONG_PARAM_GROUP("Additional Options", NULL)
    // LONG_INT_PARAM("max-ac",&max_ac, "Maximum sample allele count to consider")
    // LONG_INT_PARAM("min-ac",&min_ac, "Minimum sample allele count to consider")
    LONG_INT_PARAM("verbose",&verbose, "Periodic report when scanning VCF")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output file")

  END_LONG_PARAMS();
  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  if (inVcf.empty() || inQuery.empty()) {
    error("[E:%s:%d %s] --in-vcf, --in-query are required parameters",__FILE__,__LINE__,__FUNCTION__);
  }

  // bcf reader
  std::vector<GenomeInterval> intervals;
  parse_intervals(intervals, "", "");
  BCFOrderedReader odr(inVcf, intervals);
  bcf1_t* iv = bcf_init();
  int32_t nsamples = bcf_hdr_nsamples(odr.hdr);
  char** id_ptr = bcf_hdr_get_samples(odr.hdr);
  std::vector<std::string> id_samples(id_ptr, id_ptr+nsamples);
  std::map<std::string, int32_t> id_index;
  for (int32_t i = 0; i < nsamples; ++i)
    id_index[id_samples[i]] = i;

  // handle filter string
  std::string filter_str;
  int32_t filter_logic = 0;
  if ( vfilt.include_expr.empty() ) {
    if ( vfilt.exclude_expr.empty() ) {
      // do nothing
    }
    else {
      filter_str = vfilt.exclude_expr;
      filter_logic |= FLT_EXCLUDE;
    }
  }
  else {
    if ( vfilt.exclude_expr.empty() ) {
      filter_str = vfilt.include_expr;
      filter_logic |= FLT_INCLUDE;
    }
    else {
      error("[E:%s:%d %s] Cannot use both --include-expr and --exclude-expr options",__FILE__,__LINE__,__FUNCTION__);
    }
  }

  filter_t* filt = NULL;
  if ( filter_logic != 0 )
    filter_init(odr.hdr, filter_str.c_str());

  // handle --apply-filtrs
  std::vector<int32_t> req_flt_ids;
  if ( !vfilt.required_filters.empty() ) {
    for(int32_t i=0; i < (int32_t)vfilt.required_filters.size(); ++i) {
      req_flt_ids.push_back(bcf_hdr_id2int(odr.hdr, BCF_DT_ID, vfilt.required_filters[i].c_str()));
    }
  }

  std::ifstream ifs;
  std::string line, header;
  std::vector<std::string> words;
  std::vector<int32_t> idvec(2,0);
  ifs.open(inQuery, std::ifstream::in);
  if (!ifs.is_open()) {
    error("Query file cannot be opened");
  }
  std::getline(ifs, header);
  header.erase(std::remove(header.begin(), header.end(), '\n'), header.end());
  header += "\tInformativeSites\th1\th1_guess\th2\th2_guess\tcenter_bp_1\tcenter_bp_2\tconfig_right_1\tconfig_left_1\tconfig_right_2\tconfig_left_2\n";

  // Output
  FILE * wf;
  wf = fopen(out.c_str(), "w");
  fputs(header.c_str(), wf);

  int32_t pos=0, nRec=0, st, ed;
  int32_t* p_gt = NULL;
  int32_t  n_gt = 0;
  int32_t* info_ac = NULL;
  int32_t* info_an = NULL;
  int32_t  n_ac = 0, n_an = 0;

  while(std::getline(ifs, line)) {
    words.clear();
    split(words, "\t", line);
    if (!str2int32(words[1], pos) ) {continue;}
    if (!str2int32(words[5], st) || !str2int32(words[6], ed)) {continue;}
    auto ptr1 = id_index.find(words[3]);
    auto ptr2 = id_index.find(words[4]);
    if (ptr1 == id_index.end() || ptr2 == id_index.end()) {continue;}
    std::vector<int32_t> indx{ptr1->second, ptr2->second};

    // std::string query_samples = words[3] + "," + words[4];
    reg = words[0] + ":" + words[5] + "-" + words[6];
    parse_intervals(intervals, "", reg.c_str());
    odr.jump_to_interval(intervals[0]);
    // bcf_hdr_set_samples(odr.hdr, query_samples.c_str(), 0);
    bcf_clear(iv);

    notice("Processing the %d-th region at %s.", nRec+1, reg.c_str());

    std::vector<int32_t> h,loc,mac;
    std::vector<int32_t> center_match(2,0);
    std::vector<int32_t> guess_hap(2,0);
    std::vector< std::vector<int32_t> > g_match(2, std::vector<int32_t>(0) );
    int32_t index_focal = -1, h1=0, h2=0;

    for(int32_t k=0; odr.read(iv); ++k) {
      // unpack FILTER column
      bcf_unpack(iv, BCF_UN_FLT);
      // check --apply-filters
      bool has_filter = req_flt_ids.empty() ? true : false;
      if ( ! has_filter ) {
        for(int32_t i=0; i < iv->d.n_flt; ++i) {
          for(int32_t j=0; j < (int32_t)req_flt_ids.size(); ++j) {
            if ( req_flt_ids[j] == iv->d.flt[i] )
              has_filter = true;
          }
        }
      }
      if ( ! has_filter ) { continue; }
      // check filter logic
      if ( filt != NULL ) {
        int32_t ret = filter_test(filt, iv, NULL);
        if ( filter_logic == FLT_INCLUDE ) { if ( !ret)  has_filter = false; }
        else if ( ret ) { has_filter = false; }
      }
      if ( ! has_filter ) { continue; }
      // extract genotype and apply genotype level filter
      if ( bcf_get_genotypes(odr.hdr, iv, &p_gt, &n_gt) < 0 ) {
        error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);
      }
      if (bcf_get_info_int32(odr.hdr, iv, "AC", &info_ac, &n_ac) < 0) {continue;}
      if (bcf_get_info_int32(odr.hdr, iv, "AN", &info_an, &n_an) < 0) {continue;}
      int32_t ac = 0, an = 0;
      ac = info_ac[0]; an = info_an[0];
      if (ac > an/2) {
        ac = an - ac;
      }
      if (ac==1) {continue;}

      bool flag = 0;
      std::vector<int32_t> vec(4,0);
      int32_t v_sum = 0;
      for (int32_t i = 0; i < 2; ++i) {
        for (int32_t j = 0; j < 2; ++j) {
          if (bcf_gt_is_missing( p_gt[indx[i]*2+j] )) {
            flag = 1;
          } else {
            vec[i*2+j] = bcf_gt_allele(p_gt[indx[i]*2+j]);
          }
          v_sum += vec[i*2+j];
        }
      }
      if (iv->pos+1 == pos) {
        index_focal = h.size();
        h1 = ( (vec[0] == 1 ) ? 1 : 2 );
        h2 = ( (vec[2] == 1 ) ? 1 : 2 );
        continue;
      }
      if (flag || (v_sum == 0) || (v_sum == 4)) {continue;}
      if ((vec[0]==vec[1]) && (vec[2]==vec[3]) && (vec[0] != vec[2])) {continue;}
      if ((vec[0]!=vec[1]) && (vec[2]!=vec[3])) {continue;}

      loc.push_back(iv->pos+1);
      mac.push_back(ac);
      if (vec[0] == vec[1]) {
        h.push_back(vec[0]);
        g_match[0].push_back(3);
        g_match[1].push_back( (vec[2] == vec[0]) ? 1 : 2 );
      } else {
        h.push_back(vec[2]);
        g_match[1].push_back(3);
        g_match[0].push_back( (vec[2] == vec[0]) ? 1 : 2 );
      }
    }

    notice("Processed the %d-th region at %s, found %d informative loci.", nRec+1, reg.c_str(), h.size());

    if (h.size() == 0) {continue;}
    if (index_focal < 0 || index_focal >= (int32_t) h.size()) {
      notice("Didn't find focal locus");
      continue;
    }
    int32_t b = 0;
    int32_t d = 0, c = 0, m = -1, p = -1;
    int32_t it = index_focal;
    int32_t pre_h = 0;
    int32_t last_informative_it = 0, current_start_it = 0;
    int32_t block0_start_it = -1, block0_end_it = -1;
    std::vector< std::vector<int32_t> > blocks;
    std::vector< std::string > out_config;
    std::stringstream ss;

    for (int32_t i = 0; i < 2; ++i) {
      blocks.clear();
      b=0;d=0;c=0;m=-1;p=-1;
      it = index_focal;
      // Block zero
      while (it >= 0) {
        if (g_match[i][it] != 3) {
          break;
        }
        it--;
      }
      if (g_match[i][it] != 3) {
        pre_h = g_match[i][it];
        last_informative_it = it;
        current_start_it = it;
        while(it >= 0) {
          if (g_match[i][it] != 3) {
            if (g_match[i][it] == pre_h) {
              current_start_it = it;
            } else {
              break;
            }
          }
          it--;
        }
        block0_start_it = current_start_it;

      } else { // No informative locus upstream

        it = index_focal;
        while (it < (int32_t) loc.size()) {
          if (g_match[i][it] != 3) {
            break;
          }
          it++;
        }
        if (g_match[i][it] != 3) {
          pre_h = g_match[i][it];
          last_informative_it = it;
          current_start_it = it;
          block0_start_it = current_start_it;
        } else { // No informative locus downstream either
          center_match[i] = 0;
          guess_hap[i] = 0;
          out_config.push_back(".");
          out_config.push_back(".");
          continue;
        }
      }

      it = last_informative_it;
      while (it < (int32_t) loc.size()) {
        if (g_match[i][it] != 3) {
          if (g_match[i][it] == pre_h) {
            last_informative_it = it;
          } else {
            break;
          }
        }
        it++;
      }
      block0_end_it = last_informative_it;
      center_match[i] = loc[block0_end_it] - loc[block0_start_it] + 1;
      guess_hap[i] = pre_h;

      // Downstream
      b = 0;
      it = last_informative_it;
      while (it < (int32_t) loc.size()) {
        if (g_match[i][it] != 3) {
          if (g_match[i][it] == pre_h) {
            last_informative_it = it;
          } else { // End one block, start a new one
            d = loc[last_informative_it] - loc[current_start_it] + 1;
            c = loc[it] - 1 - loc[last_informative_it];
            if (last_informative_it == current_start_it) {
              m = mac[current_start_it];
              p = loc[current_start_it];
            }
            std::vector<int32_t> rec = {b,pre_h,d,c,m,p};
            blocks.push_back(rec);
            d=0;c=0;m=-1;p=-1;
            b++;
            current_start_it = it;
            last_informative_it = it;
            pre_h = g_match[i][it];
          }
        }
        it++;
      }
      // Last block
      d = loc[last_informative_it] - loc[current_start_it] + 1;
      c = ed - 1 - loc[last_informative_it];
      if (last_informative_it == current_start_it) {
        m = mac[current_start_it];
        p = loc[current_start_it];
      }
      blocks.push_back( std::vector<int32_t> {b,pre_h,d,c,m,p} );
      ss.str(std::string());
      for (auto & v : blocks) {
        for (auto & u : v) {
          ss << u << ",";
        }
        ss.seekp(-1,ss.cur); ss<< ';';
      }
      out_config.push_back(ss.str());
      blocks.clear();

      // Upstream
      d=0;c=0;m=-1;p=-1;
      b=0;
      current_start_it = block0_end_it;
      last_informative_it = block0_start_it;
      it =  last_informative_it;
      pre_h = g_match[i][current_start_it];
      while (it >= 0) {
        if (g_match[i][it] != 3) {
          if (g_match[i][it] == pre_h) {
            last_informative_it = it;
          } else { // End one block, start a new one
            d = loc[current_start_it] - loc[last_informative_it] + 1;
            c = loc[last_informative_it] - loc[it] - 1;
            if (last_informative_it == current_start_it) {
              m = mac[current_start_it];
              p = loc[current_start_it];
            }
            blocks.push_back( std::vector<int32_t> {b,pre_h,d,c,m,p} );
            d=0;c=0;m=-1;p=-1;
            b--;
            current_start_it = it;
            last_informative_it = it;
            pre_h = g_match[i][it];
          }
        }
        it--;
      }
      // Last block
      d = loc[current_start_it] - loc[last_informative_it] + 1;
      c = loc[last_informative_it] - st - 1;
      if (last_informative_it == current_start_it) {
        m = mac[current_start_it];
        p = loc[current_start_it];
      }
      std::vector<int32_t> rec = {b,pre_h,d,c,m,p};
      blocks.push_back(rec);
      ss.str(std::string());
      for (auto & v : blocks) {
        for (auto & u : v) {
          ss << u << ",";
        }
        ss.seekp(-1,ss.cur); ss<< ';';
      }
      out_config.push_back(ss.str());
      blocks.clear();
  }

line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
fprintf(wf, "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\n", line.c_str(), (int32_t) loc.size(),h1,guess_hap[0],h2,guess_hap[1],center_match[0],center_match[1], out_config[0].c_str(),out_config[1].c_str(),out_config[2].c_str(),out_config[3].c_str());

nRec++;

// std::cout << line;

// std::cout << "\ng1 - downstream\n";
// for ( int32_t i = index_focal; i < loc.size(); ++i ) {
// std::cout << g_match[0][i];
// }
// std::cout << "\ng1 - upstream\n";
// for ( int32_t i = index_focal; i > 0; --i ) {
// std::cout << g_match[0][i];
// }

// std::cout<< "\nDownstream config - 1\n";
// std::cout << out_config[0];
// std::cout<< "\nUpstream config - 1\n";
// std::cout << out_config[1] << '\n';

// std::cout << "\ng2 - downstream\n";
// for ( int32_t i = index_focal; i < loc.size(); ++i ) {
// std::cout << g_match[1][i];
// }
// std::cout << "\ng2 - upstream\n";
// for ( int32_t i = index_focal; i > 0; --i ) {
// std::cout << g_match[1][i];
// }
// std::cout<< "\nDownstream config - 2\n";
// std::cout << out_config[2];
// std::cout<< "\nUpstream config - 2\n";
// std::cout << out_config[3] << '\n';

// std::cout << h1 << '\t' << h2 << '\n';

// if (nRec>50)
//   break;

  }
  ifs.close();
  fclose(wf);
  return 0;
}














