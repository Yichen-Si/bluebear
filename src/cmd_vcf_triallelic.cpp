#include "cramore.h"
#include "bcf_filter_arg.h"
#include "bcf_ordered_reader.h"
#include "utils.h"

// Goal: get full list of multi-allelic sites and their k-mer context

struct SNPbuff {
  int32_t pos, ac;
  std::string mutation;
  SNPbuff(int32_t _p, int32_t _x, std::string _m) : pos(_p), ac(_x), mutation(_m) {}
  SNPbuff() : pos(0), ac(0), mutation("") {}
  void update (int32_t _p, int32_t _x, std::string _m) {
    pos = _p; ac = _x; mutation = _m;
  }
};

void FlipNucleotide(std::string& org, std::string& img, std::map<char,char>& dic) {
  std::stringstream newstr;
  for (char &c : org)
    newstr << dic[c];
  img = newstr.str();
}

void FlipNucleotide(std::string& org, std::map<char,char>& dic) {
  std::stringstream newstr;
  for (char &c : org)
    newstr << dic[c];
  org = newstr.str();
}

int32_t MultiAllelicSites(int32_t argc, char** argv) {
  std::string inVcf, out, outf, reg;
  int32_t verbose = 10000;
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

  outf = out+".triallelic.ctx."+std::to_string(kmer)+".txt";
  htsFile* wf = hts_open(outf.c_str(), "w");

  int32_t nskip = 0, nmono = 0;
  int32_t *info_ac = NULL;
  int32_t *info_an = NULL;
  int32_t n_ac = 0, n_an = 0, ac = 0, an = 0;
  int32_t nVariant = 0;
  int32_t nsamples = 0;
  char *ctx = NULL;
  int32_t nctx = 0;

  std::set<char> alphabet {'A','T','C','G'};
  std::map<char,char> bpair = { {'A','T'}, {'T','A'}, {'C','G'}, {'G','C'}, {'-','-'} };
  // std::vector<SNPbuff> waitlist;
  SNPbuff buff;

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
    if (!bcf_is_snp(iv)) {nmono++; continue;}
    if (bcf_get_info_int32(odr.hdr, iv, "AC", &info_ac, &n_ac) < 0) {continue;}
    if (bcf_get_info_int32(odr.hdr, iv, "AN", &info_an, &n_an) < 0) {continue;}
    ac = info_ac[0]; an = info_an[0];
    if (ac == 0) {nmono++; continue;}
    if (bcf_get_info_string(odr.hdr, iv, "CONTEXT", &ctx, &nctx) < 0) {continue;}

    std::string ref(1, iv->d.allele[0][0]);
    std::string alt(1, iv->d.allele[1][0]);
    char maj = iv->d.allele[0][0];
    std::string m1, m2;
    m1 = ref+alt;
    if (ac > an/2) {
      ac = an-ac;
      m1 = alt+ref;
      maj = alt[0];
    }

    if (iv->pos != buff.pos) {
      buff.update(iv->pos,ac,m1);
      continue;
    }

    if (buff.mutation[0] != m1[0]) { // Error. Different major allele?
      notice("Error. Different major alleles? %s:%d\t%s\t%s", bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1, m1.c_str(), buff.mutation.c_str());
      buff.update(iv->pos,ac,m1);
      continue;
    }
    std::string kctx(ctx);
    if (kmer < ((int32_t) kctx.size())) {
      int32_t offset = (kctx.size()-kmer)/2;
      kctx = kctx.substr(offset, kmer);
    }
    bool miss = 0;
    for (auto & v : kctx) {
      if (alphabet.find(v) == alphabet.end()) {
        miss = 1;
        break;
      }
    }
    if (miss) {
      continue;
    }
    int32_t bac = buff.ac;
    m2 = buff.mutation;
    buff.update(iv->pos,ac,m1);

    if (maj != 'A' && maj != 'C') {
      maj = bpair[maj];
      FlipNucleotide(kctx, bpair);
      FlipNucleotide(m1, bpair);
      FlipNucleotide(m2, bpair);
    }
    kctx = kctx.substr(0,(kmer-1)/2)+maj+kctx.substr((kmer-1)/2+1);

    int32_t ac1=0,ac2=0,ac3=0;
    if (maj == 'A') {
      if (m1 == "AG") {
        ac1=ac;
        if (m2 == "AC") {
          ac2=bac;
        } else {
          ac3=bac;
        }
      } else if (m2 == "AG") {
        ac1 = bac;
        if (m1 == "AC") {
          ac2 = ac;
        } else {
          ac3 = ac;
        }
      } else {
        if (m1 == "AC") {
          ac2 = ac; ac3 = bac;
        } else {
          ac2 = bac; ac3 = ac;
        }
      }
    } else {
      if (m1 == "CT") {
        ac1=ac;
        if (m2 == "CA") {
          ac2=bac;
        } else {
          ac3=bac;
        }
      } else if (m2 == "CT") {
        ac1 = bac;
        if (m1 == "CA") {
          ac2 = ac;
        } else {
          ac3 = ac;
        }
      } else {
        if (m1 == "CA") {
          ac2 = ac; ac3 = bac;
        } else {
          ac2 = bac; ac3 = ac;
        }
      }
    }
    hprintf(wf, "%s\t%d\t%s\t%c\t%d\t%d\t%d\n", bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1, kctx.c_str(), maj, ac1, ac2, ac3);
    nVariant++;
  }

  notice("Detected %d tri-allelic markers across %d samples, Skipping %d filtered markers and %d uninformative markers", nVariant, nsamples, nskip, nmono);

  hts_close(wf);

  return 0;
}
















