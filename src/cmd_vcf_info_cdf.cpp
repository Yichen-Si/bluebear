#include "cramore.h"
#include "bcf_filter_arg.h"
#include "bcf_ordered_reader.h"
#include "utils.h"


// Goal: summarize empirical CDF of ibs0 length in vcf INFO
// TODO: context type specific; generic for numerical INFO

int32_t CDFInfo(int32_t argc, char** argv) {
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
















