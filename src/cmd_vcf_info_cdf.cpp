#include "cramore.h"
#include "bcf_filter_arg.h"
#include "bcf_ordered_reader.h"
#include "utils.h"


// Goal: summarize empirical CDF of ibs0 length in vcf INFO
// TODO: context type specific; generic for numerical INFO

int32_t CDFInfoInt(int32_t argc, char** argv) {
  std::string inVcf, out, outf, reg;
  int32_t verbose = 10000;
  int32_t binsize = 2000;
  int32_t maxval = (int32_t) 5e6;
  int32_t minval = 0;
  int32_t maxac = 10;
  std::string tag = "AvgDist_bp";

  bcf_vfilter_arg vfilt;
  bcf_gfilter_arg gfilt;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Sites", NULL)
    LONG_STRING_PARAM("in-vcf",&inVcf, "Input VCF/BCF file")
    LONG_STRING_PARAM("region",&reg,"Genomic region to focus on")

    LONG_PARAM_GROUP("Variant Filtering Options", NULL)
    LONG_MULTI_STRING_PARAM("apply-filter",&vfilt.required_filters, "Require at least one of the listed FILTER strings")
    LONG_STRING_PARAM("include-expr",&vfilt.include_expr, "Include sites for which expression is true")
    LONG_STRING_PARAM("exclude-expr",&vfilt.exclude_expr, "Exclude sites for which expression is true")

    LONG_PARAM_GROUP("Additional Options", NULL)
    // LONG_INT_PARAM("kmer",&kmer, "The kmer to be considered. Cannot exceed what has been annotated in the VCF")
    LONG_INT_PARAM("bin-size",&binsize, "Bin size for empirical CDF")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output VCF file name")
    LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")

  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  // sanity check of input arguments
  if ( inVcf.empty() || out.empty() ) {
    error("[E:%s:%d %s] --in-vcf, --out are required parameters",__FILE__,__LINE__,__FUNCTION__);
  }

  int32_t nbins = (maxval - minval) / binsize;

  std::vector<int64_t> tmp (nbins,0);
  std::vector<int64_t> tvcdf[maxac+1];
  std::vector<int64_t> tscdf[maxac+1];
  std::vector<int64_t> cpgcdf[maxac+1];
  for (int32_t i = 0; i <= maxac; ++i) {
    tvcdf[i] = tmp;
    tscdf[i] = tmp;
    cpgcdf[i] = tmp;
  }

  // bcf reader
  std::vector<GenomeInterval> intervals;
  if ( !reg.empty() ) {
    parse_intervals(intervals, "", reg);
  }
  BCFOrderedReader odr(inVcf, intervals);
  bcf1_t* iv = bcf_init();

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

  int32_t *info_ac = NULL;
  int32_t n_ac = 0, ac = 0, ndist = 0;
  char *cpg = NULL;
  int32_t *dist = NULL;
  int32_t ncpg = 0;
  int32_t nVariant = 0;

  std::map<char,char> bpair = { {'A','T'}, {'T','A'}, {'C','G'}, {'G','C'}, {'-','-'} };
  std::set<std::string> tsset {"AG", "CT","GA","TC"};

  notice("Started reading site information from VCF file");

  for(int32_t k=0; odr.read(iv); ++k) {  // read marker
    // periodic message to user
    if ( nVariant % verbose == 0 ) {
      notice("Processing %d markers at %s:%d. Recorded %d variants", k, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1, nVariant);
    }

    // unpack FILTER column
    bcf_unpack(iv, BCF_UN_FLT);

    // check --apply-filters
    bool has_filter = req_flt_ids.empty() ? true : false;
    if ( ! has_filter ) {
      //notice("%d %d", iv->d.n_flt, (int32_t)req_flt_ids.size());
      for(int32_t i=0; i < iv->d.n_flt; ++i) {
        for(int32_t j=0; j < (int32_t)req_flt_ids.size(); ++j) {
          if ( req_flt_ids[j] == iv->d.flt[i] )
            has_filter = true;
        }
      }
    }
    if ( ! has_filter ) {continue;}

    // check filter logic
    if ( filt != NULL ) {
      int32_t ret = filter_test(filt, iv, NULL);
      if ( filter_logic == FLT_INCLUDE ) { if ( !ret)  has_filter = false; }
      else if ( ret ) { has_filter = false; }
    }
    if ( ! has_filter ) {continue;}

    if (!bcf_is_snp(iv)) {continue;}
    if (bcf_get_info_int32(odr.hdr, iv, "AC", &info_ac, &n_ac) < 0) {continue;}
    ac = info_ac[0];
    if (ac <= 1 || ac > maxac) {continue;}
    if (bcf_get_info_int32(odr.hdr, iv, tag.c_str(), &dist, &ndist) < 0) {continue;}
    if (bcf_get_info_string(odr.hdr, iv, "CpG", &cpg, &ncpg) < 0) {continue;}
    std::string mut = "";
    mut += (iv->d).allele[0][0];
    mut += (iv->d).allele[1][0];
    // std::string ref = (iv->d).allele[0][0];
    // std::string alt = (iv->d).allele[1][0];
    // std::string mut = ref+alt;
    int32_t it = dist[0]/binsize;
    if (it >= nbins) {it = nbins - 1;}
    if (it < 0) {continue;}

    if (tsset.find(mut) == tsset.end()) { // TV
      tvcdf[ac][it]++;
    } else { // TS
      tscdf[ac][it]++;
    }
    if (cpg[0] == 'T') { // CpG
      cpgcdf[ac][it]++;
    }
    nVariant++;
  }

  notice("Finished recording %d variants", nVariant);

  outf = out+".ibs0.avgbtw.cdf.txt";
  std::ofstream wf(outf.c_str(), std::ios::out);

  for (int32_t ac = 2; ac <= maxac; ++ac) {
    for (int32_t i = 0; i < nbins; ++i) {
      wf << ac << "\tTV\t" << (i+1)*binsize << '\t' << tvcdf[ac][i] << '\n';
      wf << ac << "\tTS\t" << (i+1)*binsize << '\t' << tscdf[ac][i] << '\n';
      wf << ac << "\tCpG\t" << (i+1)*binsize << '\t' << cpgcdf[ac][i] << '\n';
    }
  }
  wf.close();

  return 0;
}





// Goal: summarize kmer specific empirical CDF of ibs0 length in vcf INFO

int32_t CDFInfoFloatByKmer(int32_t argc, char** argv) {
  std::string inVcf, out, outf, reg;
  int32_t verbose = 10000;
  double binsize = 0.002;
  double maxval = 5.0;
  double minval = 0;
  int32_t maxac = 10;
  int32_t kmer = 3;
  int32_t offset = 0;
  std::string tag = "AvgDist_cM";

  bcf_vfilter_arg vfilt;
  bcf_gfilter_arg gfilt;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Sites", NULL)
    LONG_STRING_PARAM("in-vcf",&inVcf, "Input VCF/BCF file")
    LONG_STRING_PARAM("region",&reg,"Genomic region to focus on")

    LONG_PARAM_GROUP("Variant Filtering Options", NULL)
    LONG_MULTI_STRING_PARAM("apply-filter",&vfilt.required_filters, "Require at least one of the listed FILTER strings")
    LONG_STRING_PARAM("include-expr",&vfilt.include_expr, "Include sites for which expression is true")
    LONG_STRING_PARAM("exclude-expr",&vfilt.exclude_expr, "Exclude sites for which expression is true")

    LONG_PARAM_GROUP("Additional Options", NULL)
    LONG_INT_PARAM("kmer",&kmer, "The kmer to be considered. Cannot exceed what has been annotated in the VCF")
    LONG_DOUBLE_PARAM("bin-size",&binsize, "Bin size for empirical CDF")
    LONG_DOUBLE_PARAM("max-val",&maxval, "Upper bound for empirical CDF")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output VCF file name")
    LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")

  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  // sanity check of input arguments
  if ( inVcf.empty() || out.empty() ) {
    error("[E:%s:%d %s] --in-vcf, --out are required parameters",__FILE__,__LINE__,__FUNCTION__);
  }

  int32_t nbins = (maxval - minval) / binsize;

  std::vector<int64_t> tvcdf[maxac+1];
  for (int32_t i = 2; i <= maxac; ++i) {
    std::vector<int64_t> tmp (nbins,0);
    tvcdf[i] = tmp;
  }
  std::map<std::string, std::vector<std::vector<int64_t>* > > cdf;

  // mutation types
  std::vector<char> alphabet {'A','T','C','G'};
  std::set<char> alphaset(alphabet.begin(),alphabet.end());
  std::vector<std::string> foldmut{"A-T","A-C","A-G","C-A","C-G","C-T"};
  std::vector<std::string> tsmut{"A-G","C-T"};
  std::map<char,char> bpair = { {'A','T'}, {'T','A'}, {'C','G'}, {'G','C'}, {'-','-'} };

  // All possible k-mer configurations
  std::vector<std::string > motifs;
  int32_t ret = AllConfig(alphabet,kmer-1,motifs);
  for (int32_t i = 0; i < ret; ++i) {
    for (uint32_t j = 0; j < foldmut.size(); ++j) {
      std::string mtype = motifs[i].substr(0,(kmer-1)/2) + foldmut[j] + motifs[i].substr((kmer-1)/2);
      for (int32_t i = 0; i <= 1; ++i) {
        std::vector<int64_t> *tmp = new std::vector<int64_t>(0,0);
        cdf[mtype].push_back(tmp);
      }
      for (int32_t i = 2; i <= maxac; ++i) {
        std::vector<int64_t> *tmp = new std::vector<int64_t>(nbins,0);
        cdf[mtype].push_back(tmp);
      }
    }
  }

  // bcf reader
  std::vector<GenomeInterval> intervals;
  if ( !reg.empty() ) {
    parse_intervals(intervals, "", reg);
  }
  BCFOrderedReader odr(inVcf, intervals);
  bcf1_t* iv = bcf_init();

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

  int32_t *info_ac = NULL, *info_an = NULL;
  int32_t n_ac = 0, n_an = 0, ac = 0, an = 0, nctx = 0, nval = 0;
  float *val = NULL;
  char  *ctx = NULL;

  int32_t nVariant = 0;

  notice("Started reading site information from VCF file");

  for(int32_t k=0; odr.read(iv); ++k) {  // read marker
    // periodic message to user
    if ( k % verbose == 0 ) {
      notice("Processing %d markers at %s:%d. Recorded %d variants", k, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1, nVariant);
    }

    // unpack FILTER column
    bcf_unpack(iv, BCF_UN_FLT);

    // check --apply-filters
    bool has_filter = req_flt_ids.empty() ? true : false;
    if ( ! has_filter ) {
      //notice("%d %d", iv->d.n_flt, (int32_t)req_flt_ids.size());
      for(int32_t i=0; i < iv->d.n_flt; ++i) {
        for(int32_t j=0; j < (int32_t)req_flt_ids.size(); ++j) {
          if ( req_flt_ids[j] == iv->d.flt[i] )
            has_filter = true;
        }
      }
    }
    if ( ! has_filter ) {continue;}

    // check filter logic
    if ( filt != NULL ) {
      int32_t ret = filter_test(filt, iv, NULL);
      if ( filter_logic == FLT_INCLUDE ) { if ( !ret)  has_filter = false; }
      else if ( ret ) { has_filter = false; }
    }
    if ( ! has_filter ) {continue;}

    if (!bcf_is_snp(iv)) {continue;}
    if (bcf_get_info_int32(odr.hdr, iv, "AC", &info_ac, &n_ac) < 0) {continue;}
    if (bcf_get_info_int32(odr.hdr, iv, "AN", &info_an, &n_an) < 0) {continue;}
    ac = info_ac[0]; an = info_an[0];
    bool flip = 0;
    if (ac > an/2) {
      flip = 1;
      ac = an-ac;
    }
    if (ac <= 1 || ac > maxac) {continue;}
    if (bcf_get_info_float(odr.hdr, iv, tag.c_str(), &val, &nval) < 0) {continue;}
    int32_t it = val[0]/binsize;
    if (it >= nbins) {it = nbins - 1;}
    if (it < 0) {continue;}

    if (bcf_get_info_string(odr.hdr, iv, "CONTEXT", &ctx, &nctx) < 0) {continue;}
    std::string kctx(ctx);
    if (kmer < ((int32_t) kctx.size())) {
      offset = (kctx.size()-kmer)/2;
      kctx = kctx.substr(offset, kmer);
    }
    bool flag = 0;
    for (int32_t k = 0; k < kmer; ++k) {
      if ( alphaset.find(ctx[offset+k]) == alphaset.end() ) {
        flag = 1;
        break;
      }
    }
    if (flag) {continue;}
    char ref = iv->d.allele[0][0];
    char alt = iv->d.allele[1][0];
    if ( std::string(1, ref) != kctx.substr((kmer-1)/2, 1) ) {
      notice("Ref. allele != kmer center (%c/%c v.s. %s)", ref, alt, kctx.c_str());
    }
    char maj = ref, mir = alt;
    if (flip) {
      maj = alt;
      mir = ref;
    }
    std::string mtype(1, maj);
    mtype += "-";
    mtype += mir;

    if ( maj == 'T' || maj == 'G' ) {
      std::string fctx = "";
      std::string fmtype="";
      for (uint32_t k = kctx.size(); k > 0; --k) {
        fctx += bpair[kctx[k-1]];
      }
      for (uint32_t k = 0; k < mtype.size(); ++k) {
        fmtype += bpair[mtype[k]];
      }
      kctx = fctx;
      mtype = fmtype;
    }

    kctx = kctx.substr(0,(kmer-1)/2)+mtype+kctx.substr((kmer-1)/2+1);
    auto ptr = cdf.find(kctx);
    if (ptr != cdf.end()) {
      (*ptr->second[ac])[it]++;
    }
    if ( mtype != "A-G" && mtype != "C-T" ) { // TV
      tvcdf[ac][it]++;
    }
    nVariant++;
  }

  notice("Finished recording %d variants", nVariant);

  outf = out+".ibs0.cm.avgbtw.cdf.txt";
  std::ofstream wf(outf.c_str(), std::ios::out);

  for (int32_t ac = 2; ac <= maxac; ++ac) {
    for (int32_t i = 0; i < nbins; ++i) {
      wf << ac << "\tTV\t" << (i+1)*binsize << '\t' << tvcdf[ac][i] << '\n';
    }
  }
  for (auto & v : cdf) {
    for (int32_t ac = 2; ac <= maxac; ++ac) {
      for (int32_t i = 0; i < nbins; ++i) {
        wf << ac << '\t' << v.first << '\t' << (i+1)*binsize << '\t' << (*v.second[ac])[i] << '\n';
      }
    }
  }

  wf.close();

  return 0;
}
























// Goal: summarize empirical CDF of ibs0 length in vcf INFO by functional annotations

int32_t CDFInfoFloatByAnn(int32_t argc, char** argv) {
  std::string inVcf, out, outf, reg;
  int32_t verbose = 10000;
  double binsize = 0.002;
  double maxval = 3.0;
  double minval = 0;
  int32_t target_ac = 2;
  std::vector<std::string> keywords;

  // bcf_vfilter_arg vfilt;
  // bcf_gfilter_arg gfilt;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Sites", NULL)
    LONG_STRING_PARAM("in-vcf",&inVcf, "Input VCF/BCF file")
    LONG_STRING_PARAM("region",&reg,"Genomic region to focus on")

    // LONG_PARAM_GROUP("Variant Filtering Options", NULL)
    // LONG_MULTI_STRING_PARAM("apply-filter",&vfilt.required_filters, "Require at least one of the listed FILTER strings")
    // LONG_STRING_PARAM("include-expr",&vfilt.include_expr, "Include sites for which expression is true")
    // LONG_STRING_PARAM("exclude-expr",&vfilt.exclude_expr, "Exclude sites for which expression is true")

    LONG_PARAM_GROUP("Additional Options", NULL)
    LONG_INT_PARAM("ac",&target_ac,"Allele count to annotate")
    LONG_MULTI_STRING_PARAM("key-words",&keywords, "Functional annotation to summarize")
    LONG_DOUBLE_PARAM("bin-size",&binsize, "Bin size for empirical CDF")
    LONG_DOUBLE_PARAM("max-val",&maxval, "Upper bound for empirical CDF")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output VCF file name")
    LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")

  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  // sanity check of input arguments
  if ( inVcf.empty() || out.empty() ) {
    error("[E:%s:%d %s] --in-vcf, --out are required parameters",__FILE__,__LINE__,__FUNCTION__);
  }

  int32_t nbins = (maxval - minval) / binsize;
  // Type -> keywords -> count
  std::map<std::string, std::map<std::string, std::vector<int64_t> > > fn_ct;
  std::vector<std::string> type_key{"CpG","TS","TV"};
  std::vector<int64_t> cvec(nbins+1, 0);
  for (auto & v : type_key) {
    std::map<std::string, std::vector<int64_t> > tmp;
    for (auto & u : keywords) {
      tmp[u] = cvec;
    }
    tmp["Tot"] = cvec;
    fn_ct[v] = tmp;
  }

  // mutation types
  std::vector<char> alphabet {'A','T','C','G'};
  std::set<char> alphaset(alphabet.begin(),alphabet.end());
  std::set<std::string> tsmut{"A-G","C-T","G-A","T-C"};
  std::map<char,char> bpair = { {'A','T'}, {'T','A'}, {'C','G'}, {'G','C'}, {'-','-'} };

  // bcf reader
  std::vector<GenomeInterval> intervals;
  if ( !reg.empty() ) {
    parse_intervals(intervals, "", reg);
  }
  BCFOrderedReader odr(inVcf, intervals);
  bcf1_t* iv = bcf_init();

  float *info_cm = NULL;
  char  *ann = NULL;
  char  *cpg = NULL;
  float cm = 0.0;
  int32_t *info_ac = NULL;
  int32_t n_ac = 0, n_cm = 0, n_ann = 0, n_cpg = 0;
  int32_t nVariant = 0;

  notice("Started reading site information from VCF file");

  for(int32_t k=0; odr.read(iv); ++k) {  // read marker
    // periodic message to user
    bcf_unpack(iv, BCF_UN_INFO);
    if ( k % verbose == 0 ) {
      notice("Processing %d markers at %s:%d. Recorded %d variants", k, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1, nVariant);
    }
    if (!bcf_is_snp(iv)) {continue;}


    if (bcf_get_info_string(odr.hdr, iv, "ANN", &ann, &n_ann) < 0) {continue;}
    if (bcf_get_info_string(odr.hdr, iv, "CpG", &cpg, &n_cpg) < 0) {continue;}
    if (bcf_get_info_float(odr.hdr, iv, "AvgDist_cM", &info_cm, &n_cm) < 0) {continue;}
    if (bcf_get_info_int32(odr.hdr, iv, "AC", &info_ac, &n_ac) < 0) {continue;}
    if (info_ac[0] != target_ac) {continue;}
    cm = info_cm[0];
    if (cm < 0) {cm = 0.0;}

    std::string fun(ann);
    std::string type = "TV";
    char ref = iv->d.allele[0][0];
    char alt = iv->d.allele[1][0];
    std::string mtype(1, ref);
    mtype += "-";
    mtype += alt;
    if ( tsmut.find(mtype) != tsmut.end() ) {
      type = "TS";
    }
    if (cpg[0] == 'T') {
      type = "CpG";
    }

    int32_t b = (int32_t) std::round(cm / binsize);
    if (cm > maxval) {b = nbins;}
    if (cm < minval) {b = 0;}
    fn_ct[type]["Tot"][b] += 1;

    for (auto & v : keywords) {
      if ( fun.find(v) != std::string::npos ) {
        fn_ct[type][v][b] += 1;
      }
    }
    nVariant++;
  }

  notice("Finished recording %d variants", nVariant);

  outf = out;
  std::string suff = ".ibs0.cm.avgbtw.fn.ct";
  if (!reg.empty()) {
    outf += "_"+reg+suff;
  } else {
    outf += suff;
  }

  FILE *wf;
  wf = fopen(outf.c_str(),"w");
  for (auto & v : fn_ct) { // Type
    for (auto & u : v.second) { // Fn
      for (int32_t k = 0; k <= nbins; ++k) {
        fprintf(wf, "%d\t%s\t%s\t%.3f\t%ld\n", target_ac, v.first.c_str(), u.first.c_str(), (k*binsize),(long)u.second[k]);
      }
    }
  }
  fclose(wf);

  return 0;
}
























// Goal: summarize empirical CDF of ibs0 length from simulated data

int32_t CDFInfoFloatSimu(int32_t argc, char** argv) {
  std::string inVcf, out, outf, reg;
  int32_t verbose = 10000;
  double binsize = 0.002;
  double maxval = 5.0;
  double minval = 0;
  int32_t maxac = 10;
  std::string tag = "AvgDist_cM";

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Sites", NULL)
    LONG_STRING_PARAM("in-vcf",&inVcf, "Input VCF/BCF file")
    LONG_STRING_PARAM("region",&reg,"Genomic region to focus on")

    LONG_PARAM_GROUP("Additional Options", NULL)
    LONG_DOUBLE_PARAM("bin-size",&binsize, "Bin size for empirical CDF")
    LONG_DOUBLE_PARAM("max-val",&maxval, "Upper bound for empirical CDF")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output VCF file name")
    LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")

  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  // sanity check of input arguments
  if ( inVcf.empty() || out.empty() ) {
    error("[E:%s:%d %s] --in-vcf, --out are required parameters",__FILE__,__LINE__,__FUNCTION__);
  }

  int32_t nbins = (maxval - minval) / binsize;

  std::vector<int64_t> pmcdf[maxac+1];
  std::vector<int64_t> tvcdf[maxac+1];

  for (int32_t i = 2; i <= maxac; ++i) {
    std::vector<int64_t> tmp (nbins,0);
    tvcdf[i] = tmp;
    pmcdf[i] = tmp;
  }

  // bcf reader
  std::vector<GenomeInterval> intervals;
  if ( !reg.empty() ) {
    parse_intervals(intervals, "", reg);
  }
  BCFOrderedReader odr(inVcf, intervals);
  bcf1_t* iv = bcf_init();

  int32_t *info_ac = NULL, *pm_config = NULL;
  int32_t n_ac = 0, ac = 0, nval = 0, npm = 0;
  double *val = NULL;

  int32_t nVariant = 0;

  notice("Started reading site information from VCF file");

  for(int32_t k=0; odr.read(iv); ++k) {  // read marker
    // periodic message to user
    if ( k % verbose == 0 ) {
      notice("Processing %d markers at %s:%d. Recorded %d variants", k, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1, nVariant);
    }

    if (bcf_get_info_int32(odr.hdr, iv, "AC", &info_ac, &n_ac) < 0) {continue;}
    ac = info_ac[0];
    if (ac <= 1 || ac > maxac) {continue;}
    if (bcf_get_info_float(odr.hdr, iv, tag.c_str(), &val, &nval) < 0) {continue;}
    int32_t it = val[0]/binsize;
    if (it >= nbins) {it = nbins - 1;}
    if (it < 0) {continue;}

    if (bcf_get_info_int32(odr.hdr, iv, "PMconfig", &pm_config, &npm) < 0) {
      tvcdf[ac][it]++;
    } else {
      pmcdf[ac][it]++;
    }
    nVariant++;
  }

  notice("Finished recording %d variants", nVariant);

  outf = out+".ibs0.cm.avgbtw.cdf.txt";
  std::ofstream wf(outf.c_str(), std::ios::out);

  for (int32_t ac = 2; ac <= maxac; ++ac) {
    for (int32_t i = 0; i < nbins; ++i) {
      wf << ac << "\tNormal\t" << (i+1)*binsize << '\t' << tvcdf[ac][i] << '\n';
      wf << ac << "\tPM\t" << (i+1)*binsize << '\t' << pmcdf[ac][i] << '\n';
    }
  }
  wf.close();

  return 0;
}











