#include "cramore.h"
#include "bcf_filter_arg.h"
#include "bcf_ordered_reader.h"
#include "utils.h"

// Goal -- get Folded SFS

int32_t cmdVcfSFS(int32_t argc, char** argv) {
  std::string inVcf, out, outf, reg;
  int32_t verbose = 10000;
  int32_t use_info = 1;

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
    LONG_INT_PARAM("use-info",&use_info, "Whether to use AC/AN in info field")

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

  int32_t nVariant = 0;
  int32_t nsamples = bcf_hdr_nsamples(odr.hdr);
  int32_t sfs[nsamples+1] = {0};

  notice("Started Reading site information from VCF file, identifying %d samples", nsamples);

  if ( nsamples == 0 )
    error("FATAL ERROR: The VCF does not have any samples with genotypes");

  int32_t* p_gt = NULL;
  int32_t  n_gt = 0;
  int32_t nskip = 0, nmono = 0;
  int32_t *info_ac = NULL;
  int32_t *info_an = NULL;
  int32_t n_ac = 0, n_an = 0, ac = 0, an = 0;

  for(int32_t k=0; odr.read(iv); ++k) {  // read marker
    // periodic message to user
    if ( k % verbose == 0 )
      notice("Processing %d markers at %s:%d. Skipped %d filtered markers and %d uninformative markers, retaining %d variants", k, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1, nskip, nmono, nVariant);

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
    if ( ! has_filter ) { ++nskip; continue; }
    // check filter logic
    if ( filt != NULL ) {
      int32_t ret = filter_test(filt, iv, NULL);
      if ( filter_logic == FLT_INCLUDE ) { if ( !ret)  has_filter = false; }
      else if ( ret ) { has_filter = false; }
    }
    if ( ! has_filter ) { ++nskip; continue; }

    if (use_info) {
      if (bcf_get_info_int32(odr.hdr, iv, "AC", &info_ac, &n_ac) < 0) {continue;}
      if (bcf_get_info_int32(odr.hdr, iv, "AN", &info_an, &n_an) < 0) {continue;}
      ac = info_ac[0]; an = info_an[0];
      if (ac > an/2) {
        ac = an - ac;
      }
      if (ac == 0) {nmono++; continue;}
      sfs[(uint32_t) ac] ++;
      nVariant++;
    } else {
      // extract genotype
      if ( bcf_get_genotypes(odr.hdr, iv, &p_gt, &n_gt) < 0 ) {
        error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);
      }
      ac = 0; an=0;
      for(int32_t i=0; i < nsamples; ++i) {
        int32_t g1 = p_gt[2*i];
        int32_t g2 = p_gt[2*i+1];
        if ( !bcf_gt_is_missing(g1) && !bcf_gt_is_missing(g2) ) {
          ac += ((bcf_gt_allele(g1) > 0) ? 1 : 0) + ((bcf_gt_allele(g2) > 0) ? 1 : 0);
          an ++;
        }
      }
      if (ac > an/2) {
        ac = an - ac;
      }
      sfs[(uint32_t) ac] ++;
      nVariant++;
    }
  }

  notice("Finished Processing %d markers across %d samples, Skipping %d filtered markers and %d uninformative markers", nVariant, nsamples, nskip, nmono);

  outf = out+".folded.sfs";
  htsFile* wf = hts_open(outf.c_str(), "w");
  int32_t i = 0;
  for (i = 0; i < nsamples+1; i++) {
    hprintf(wf,"%d\t%d\n",i,sfs[i]);
  }
  hts_close(wf);

  return 0;
}












// Goal: context specific SFS (kmer contex defined by major allele)

int32_t KmerSFS(int32_t argc, char** argv) {
  std::string inVcf, out, outf, reg;
  int32_t verbose = 10000;
  int32_t use_info = 1;
  int32_t kmer = 5;
  std::map<std::string, std::vector<int64_t> > sfs;

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
    LONG_INT_PARAM("use-info",&use_info, "Whether to use AC/AN in info field")
    LONG_INT_PARAM("kmer",&kmer, "The kmer to be considered. Cannot exceed what has been annotated in the VCF")

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

  int32_t nskip = 0, nmono = 0;
  int32_t *info_ac = NULL;
  int32_t *info_an = NULL;
  int32_t n_ac = 0, n_an = 0, ac = 0, an = 0;
  int32_t nVariant = 0;
  int32_t nsamples = 0;

  // Need a way to set up sfs size when the VCF does not contain GT
  int32_t tmp = 0;
  BCFOrderedReader odrtmp(inVcf, intervals);
  while (odrtmp.read(iv)) {
    if (bcf_get_info_int32(odrtmp.hdr, iv, "AN", &info_an, &n_an)) {
      nsamples = std::max(nsamples,info_an[0]);
      tmp ++;
      if (tmp > 500)
        break;
    }
  }
  nsamples = (nsamples / 2 / 500) * 500 + 500;
  notice("Will initialize SFS to at most AC = %d.",nsamples);

  // Setup kmer sfs
  std::vector<char> alphabet {'A','T','C','G'};
  int32_t nletters = 4;
  std::vector<std::string> allmut;
  for (int32_t i = 0; i < nletters-1; ++i) {
    for (int32_t j = i+1; j < nletters; ++j) {
      std::stringstream mtype;
      mtype << alphabet[i] << "-" << alphabet[j];
      allmut.push_back(mtype.str());
      std::stringstream ctype;
      ctype << alphabet[j] << "-" << alphabet[i];
      allmut.push_back(ctype.str());
    }
  }

// std::cout << allmut.size() << '\t';
// for (auto & v : allmut) {
//   std::cout << v << '\t';
// }
// std::cout << '\n';

  std::vector<std::string > motifs;
  int32_t ret = AllConfig(alphabet,kmer-1,motifs);
  for (int32_t i = 0; i < ret; ++i) {
    for (uint32_t j = 0; j < allmut.size(); ++j) {
      std::string mtype = motifs[i].substr(0,(kmer-1)/2) + allmut[j] + motifs[i].substr((kmer-1)/2);
      std::vector<int64_t> tmp(nsamples+1, 0);
      sfs[mtype] = tmp;
    }
  }

// std::cout << motifs.size() << '\t';
// for (auto & v : motifs) {
//   std::cout << v << '\t';
// }
// std::cout << '\n';

  notice("Started reading site information from VCF file");

  for(int32_t k=0; odr.read(iv); ++k) {  // read marker
    // periodic message to user
    if ( (k+1) % verbose == 0 ) {
      notice("Processing %d markers at %s:%d. Skipped %d filtered markers and %d uninformative markers, retaining %d variants", k, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1, nskip, nmono, nVariant);
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
    if ( ! has_filter ) { ++nskip; continue; }

    // check filter logic
    if ( filt != NULL ) {
      int32_t ret = filter_test(filt, iv, NULL);
      if ( filter_logic == FLT_INCLUDE ) { if ( !ret)  has_filter = false; }
      else if ( ret ) { has_filter = false; }
    }
    if ( ! has_filter ) { ++nskip; continue; }

    if (bcf_get_info_int32(odr.hdr, iv, "AC", &info_ac, &n_ac) < 0) {continue;}
    if (bcf_get_info_int32(odr.hdr, iv, "AN", &info_an, &n_an) < 0) {continue;}
    ac = info_ac[0]; an = info_an[0];
    if (ac == 0) {nmono++; continue;}
    if (!bcf_is_snp(iv)) {nmono++; continue;}

    char *ctx = NULL;
    int32_t nctx = 0;
    if (bcf_get_info_string(odr.hdr, iv, "CONTEXT", &ctx, &nctx) < 0) {continue;}
    std::string kctx(ctx);
    if (kmer < ((int32_t) kctx.size())) {
      int32_t offset = (kctx.size()-kmer)/2;
      kctx = kctx.substr(offset, kmer);
    }
    char ref = iv->d.allele[0][0];
    char alt = iv->d.allele[1][0];
    if (ac > an/2) {
      ac = an - ac;
      kctx = kctx.substr(0,(kmer-1)/2)+alt+"-"+ref+kctx.substr((kmer-1)/2+1);
    } else {
      kctx = kctx.substr(0,(kmer-1)/2)+ref+"-"+alt+kctx.substr((kmer-1)/2+1);
    }
    auto ptr = sfs.find(kctx);
    if (ptr != sfs.end()) {
      (ptr->second)[ac] ++;
      nVariant++;
    } else {
      nmono++;
    }
  }

  notice("Finished Processing %d markers across %d samples, Skipping %d filtered markers and %d uninformative markers", nVariant, nsamples, nskip, nmono);

  outf = out+".ctx.kmer."+std::to_string(kmer)+".sfs";
  std::ofstream wf(outf.c_str(), std::ios::out);

  // Merge equivalent mutation types
  std::map<char,char> bpair = { {'A','T'}, {'T','A'}, {'C','G'}, {'G','C'}, {'-','-'} };
  std::vector<std::string> mset {"A-T","A-C","A-G","C-A","C-T","C-G"};
  std::vector<std::string> cset {"T-A","T-G","T-C","G-T","G-A","G-C"};
  for (int32_t i = 0; i < ret; ++i) {
    for (uint32_t j = 0; j < mset.size(); ++j) {
      std::string mtype = motifs[i].substr(0,(kmer-1)/2) + mset[j] + motifs[i].substr((kmer-1)/2);
      std::string ctype = "";
      for (uint32_t k = mtype.size(); k > 0; --k) {
        ctype += bpair[mtype[k-1]];
      }
      auto fs1 = sfs.find(mtype);
      auto fs2 = sfs.find(ctype);
      if (fs2 != sfs.end()) {
        for (uint32_t k = 0; k < (fs1->second).size(); ++k)
          fs1->second[k] += fs2->second[k];
      }
      wf << fs1->first;
      for (auto & w : fs1->second)
        wf << '\t' << w;
      wf << '\n';
    }
  }

  wf.close();

  return 0;
}
















