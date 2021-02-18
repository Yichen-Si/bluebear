#include "cramore.h"
#include "bcf_filter_arg.h"
#include "bcf_ordered_reader.h"
#include "utils.h"
#include <algorithm>


template<typename T>
int32_t binarySearch(std::vector<T> & arr, int32_t l, int32_t r, T x)
{
    if (r > l) {
        int32_t mid = l + (r - l) / 2;
        if (arr[mid+1] > x && arr[mid] <= x)
            return mid;
        if (arr[mid] > x)
            return binarySearch(arr, l, mid - 1, x);
        return binarySearch(arr, mid + 1, r, x);
    }
    return r;
}

// Goal: mutation rate specific SFS (kmer contex defined by major allele)

int32_t mrSFS(int32_t argc, char** argv) {
  std::string inVcf, inBin, outf, reg;
  int32_t verbose = 10000;
  int32_t use_info = 1;
  double unit = 0.00005;
  int32_t ignore_multiallele = 0;

  bcf_vfilter_arg vfilt;
  bcf_gfilter_arg gfilt;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Sites", NULL)
    LONG_STRING_PARAM("in-vcf",&inVcf, "Input VCF/BCF file")
    LONG_STRING_PARAM("in-bin",&inBin, "Input file that contains mutation rate cutoffs to use")
    LONG_STRING_PARAM("region",&reg,"Genomic region to focus on")

    LONG_PARAM_GROUP("Variant Filtering Options", NULL)
    LONG_MULTI_STRING_PARAM("apply-filter",&vfilt.required_filters, "Require at least one of the listed FILTER strings")
    LONG_STRING_PARAM("include-expr",&vfilt.include_expr, "Include sites for which expression is true")
    LONG_STRING_PARAM("exclude-expr",&vfilt.exclude_expr, "Exclude sites for which expression is true")

    LONG_PARAM_GROUP("Additional Options", NULL)
    LONG_INT_PARAM("use-info",&use_info, "Whether to use AC/AN in info field")
    LONG_INT_PARAM("bi-allele",&ignore_multiallele, "If multi-allelic, only count the first one")
    LONG_DOUBLE_PARAM("unit",&unit, "The minimum unit / gcd of the cutoffs")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &outf, "Output VCF file name")
    LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")

  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  // sanity check of input arguments
  if ( inVcf.empty() || inBin.empty() || outf.empty() ) {
    error("[E:%s:%d %s] --in-vcf, --in-bin, --out are required parameters",__FILE__,__LINE__,__FUNCTION__);
  }

  int32_t max_rare = 50;
  std::vector<double> lower{0.0,0.001,0.005,0.01,0.05};
  std::vector<double> upper{0.001,0.005,0.01,0.05, 1.00};
  std::map<char,char> bpair = { {'A','T'}, {'T','A'}, {'C','G'}, {'G','C'}, {'-','-'} };

  std::map<char, std::vector<std::string> > mutation_tag;
  std::vector<std::string> AtoX{"A_C","A_G","A_T"};
  std::vector<std::string> CtoX{"C_A","C_G","C_T"};
  mutation_tag['A'] = AtoX;
  mutation_tag['C'] = CtoX;

  std::map<int32_t, std::vector<int64_t> > sfs;
  std::vector<int32_t> mr_keys;

  // Setup sfs
  std::ifstream ifs;
  std::string line;
  std::vector<std::string> words;
  double x;
  ifs.open(inBin, std::ifstream::in);
  if (!ifs.is_open()) {
    error("Info file cannot be opened");
  }
  while(std::getline(ifs, line)) {
    words.clear();
    split(words, "\t", line);
    try {
      x = std::stod(words[0]);
    }
    catch (const std::invalid_argument& ia) {
      continue;
    }
    mr_keys.push_back((int32_t) (x/unit));
  }
  ifs.close();

  std::sort(mr_keys.begin(),mr_keys.end());
  if (mr_keys[0] > 0) {
    mr_keys.insert(mr_keys.begin(), 0);
  }

  for (auto & v : mr_keys) {
    std::vector<int64_t> tmp(max_rare+1+lower.size(), 0);
    sfs[v] = tmp;
  }
  std::vector<int64_t> tmp(max_rare+1+lower.size(), 0);
  sfs[-1] = tmp; // For unknown mutation rate

  notice("Finish set up %d mutation rate bins", mr_keys.size() );



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

  notice("Started reading site information from VCF file");

  int32_t pre_pos = 0;
  int32_t multi_allele = 0;
  for(int32_t k=0; odr.read(iv); ++k) {  // read marker
    // periodic message to user
    if ( (k+1) % verbose == 0 ) {
      notice("Processing %d markers at %s:%d. Skipped %d filtered markers and %d uninformative markers, %d multi-allelic variants, retaining %d variants", k, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1, nskip, nmono, multi_allele, nVariant);
    }
    if (!bcf_is_snp(iv)) {nmono++; continue;}
    if (ignore_multiallele) {
      if (iv->pos == pre_pos) {multi_allele++; continue;}
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
      int32_t ftest = filter_test(filt, iv, NULL);
      if ( filter_logic == FLT_INCLUDE ) { if ( !ftest)  has_filter = false; }
      else if ( ftest ) { has_filter = false; }
    }
    if ( ! has_filter ) { ++nskip; continue; }

    if (bcf_get_info_int32(odr.hdr, iv, "AC", &info_ac, &n_ac) < 0) {continue;}
    if (bcf_get_info_int32(odr.hdr, iv, "AN", &info_an, &n_an) < 0) {continue;}
    ac = info_ac[0]; an = info_an[0];
    if (ac == 0) {nmono++; continue;}

    char ref = iv->d.allele[0][0];
    char alt = iv->d.allele[1][0];
    if (ac > an/2) {
      ac = an - ac;
    }
    if ( ref != 'A' && ref != 'C' ) {
      ref = bpair[ref];
      alt = bpair[ref];
    }
    std::string mtag = "MR_";
    mtag.append(1, ref);
    mtag.append(1, '_');
    mtag.append(1, alt);

    float *info_mr = NULL;
    int32_t n_mr = 0;
    int32_t key = -1;
    if (bcf_get_info_float(odr.hdr, iv, mtag.c_str(), &info_mr, &n_mr) < 0) {
      nmono += 1;
    } else {
      int32_t w = (int32_t) (info_mr[0] / unit);
      key = binarySearch<int32_t>(mr_keys,0,mr_keys.size()-1,w);
      if (key == -1) {key = 0;}
      key = mr_keys[key];
      nVariant++;
    }

    auto ptr = sfs.find( key );
    if (ac <= max_rare) {
      (ptr->second)[ac] ++;
    }
    double af = 1.0 * ac / an;
    for (uint32_t it = 0; it < lower.size(); ++it) {
      if (af > lower[it] && af <= upper[it]) {
        (ptr->second)[max_rare + it+1]++;
        break;
      }
    }
    pre_pos = iv->pos;
  }

  notice("Finished Processing %d markers. Skipping %d filtered markers and %d uninformative markers", nVariant, nskip, nmono);

  std::ofstream wf(outf.c_str(), std::ios::out);

  for (auto & v : mr_keys) {
    auto ptr = sfs.find(v);
    if (ptr->first == -1) {
      wf << "Unknown";
    } else {
      wf << ptr->first * unit;
    }
    for (auto & w : ptr->second)
      wf << '\t' << w;
    wf << '\n';
  }
  wf.close();
  return 0;
}







// Goal: mutation rate specific SFS (kmer contex defined by major allele)

int32_t CDFInfoFloatByMR(int32_t argc, char** argv) {
  std::string inVcf, inBin, outf, out, reg;
  int32_t verbose = 10000;
  double unit = 0.00005;
  double binsize = 0.002;
  double maxval = 2.0;
  double minval = 0;
  int32_t maxac = 10;
  int32_t ignore_multiallele = 0;
  std::string tag = "AvgDist_cM";

  bcf_vfilter_arg vfilt;
  bcf_gfilter_arg gfilt;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Sites", NULL)
    LONG_STRING_PARAM("in-vcf",&inVcf, "Input VCF/BCF file")
    LONG_STRING_PARAM("in-bin",&inBin, "Input file that contains mutation rate cutoffs to use")
    LONG_STRING_PARAM("region",&reg,"Genomic region to focus on")

    LONG_PARAM_GROUP("Variant Filtering Options", NULL)
    LONG_MULTI_STRING_PARAM("apply-filter",&vfilt.required_filters, "Require at least one of the listed FILTER strings")
    LONG_STRING_PARAM("include-expr",&vfilt.include_expr, "Include sites for which expression is true")
    LONG_STRING_PARAM("exclude-expr",&vfilt.exclude_expr, "Exclude sites for which expression is true")

    LONG_PARAM_GROUP("Additional Options", NULL)
    LONG_DOUBLE_PARAM("unit",&unit, "The minimum unit / gcd of the cutoffs")
    LONG_INT_PARAM("max-ac",&maxac,"Maximun AC")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output prefix")
    LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")

  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  // sanity check of input arguments
  if ( inVcf.empty() || inBin.empty() || out.empty() ) {
    error("[E:%s:%d %s] --in-vcf, --in-bin, --out are required parameters",__FILE__,__LINE__,__FUNCTION__);
  }

  int32_t nbins = (maxval - minval) / binsize;
  std::map<char,char> bpair = { {'A','T'}, {'T','A'}, {'C','G'}, {'G','C'}, {'-','-'} };
  std::vector<std::string> AtoX{"A_C","A_G","A_T"};
  std::vector<std::string> CtoX{"C_A","C_G","C_T"};
  std::map<char, std::vector<std::string> > mutation_tag;
  mutation_tag['A'] = AtoX;
  mutation_tag['C'] = CtoX;

  // Setup cdf
  std::map<int32_t, std::vector<std::vector<double>* > > cdf;
  std::vector<int32_t> mr_keys;

  std::ifstream ifs;
  std::string line;
  std::vector<std::string> words;
  double x;
  ifs.open(inBin, std::ifstream::in);
  if (!ifs.is_open()) {
    error("Info file cannot be opened");
  }
  while(std::getline(ifs, line)) {
    words.clear();
    split(words, "\t", line);
    try {
      x = std::stod(words[0]);
    }
    catch (const std::invalid_argument& ia) {
      continue;
    }
    mr_keys.push_back((int32_t) (x/unit));
  }
  ifs.close();

  std::sort(mr_keys.begin(),mr_keys.end());
  if (mr_keys[0] > 0) {
    mr_keys.insert(mr_keys.begin(), 0);
  }

  for (auto & v : mr_keys) {
    for (int32_t i = 0; i <= 1; ++i) { // to simplify
      std::vector<int64_t> *tmp = new std::vector<double>(0,0);
      cdf[v].push_back(tmp);
    }
    for (int32_t i = 2; i <= maxac; ++i) {
      std::vector<int64_t> *tmp = new std::vector<double>(nbins,0);
      cdf[v].push_back(tmp);
    }
  }
  notice("Finish set up %d mutation rate bins", mr_keys.size() );


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
  float   *val = NULL;
  float   *info_mr = NULL;
  int32_t *info_ac = NULL;
  int32_t *info_an = NULL;
  int32_t nval = 0, n_mr = 0, n_ac = 0, n_an = 0, ac = 0, an = 0;
  int32_t nVariant = 0;

  notice("Started reading site information from VCF file");

  int32_t pre_pos = 0;
  int32_t multi_allele = 0;
  for(int32_t k=0; odr.read(iv); ++k) {  // read marker
    // periodic message to user
    if ( (k+1) % verbose == 0 ) {
      notice("Processing %d markers at %s:%d. Skipped %d filtered markers and %d uninformative markers, %d multi-allelic variants, retaining %d variants", k, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1, nskip, nmono, multi_allele, nVariant);
    }
    if (!bcf_is_snp(iv)) {++nskip; continue;}
    if (ignore_multiallele) {
      if (iv->pos == pre_pos) {multi_allele++; continue;}
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
      int32_t ftest = filter_test(filt, iv, NULL);
      if ( filter_logic == FLT_INCLUDE ) { if ( !ftest)  has_filter = false; }
      else if ( ftest ) { has_filter = false; }
    }
    if ( ! has_filter ) { ++nskip; continue; }

    if (bcf_get_info_int32(odr.hdr, iv, "AC", &info_ac, &n_ac) < 0) {continue;}
    if (bcf_get_info_int32(odr.hdr, iv, "AN", &info_an, &n_an) < 0) {continue;}
    ac = info_ac[0]; an = info_an[0];
    if (ac > an/2) {ac = an-ac;}
    if (ac <= 1 || ac > maxac) {nmono++; continue;}
    if (bcf_get_info_float(odr.hdr, iv, tag.c_str(), &val, &nval) < 0) {
      nmono += 1; continue;}

    int32_t it = val[0]/binsize;
    if (it >= nbins) {it = nbins - 1;}
    if (it < 0) {continue;}

    char ref = iv->d.allele[0][0];
    char alt = iv->d.allele[1][0];
    if (ac > an/2) {
      ac = an - ac;
    }
    if ( ref != 'A' && ref != 'C' ) {
      ref = bpair[ref];
      alt = bpair[ref];
    }
    std::string mtag = "MR_";
    mtag.append(1, ref);
    mtag.append(1, '_');
    mtag.append(1, alt);

    if (bcf_get_info_float(odr.hdr, iv, mtag.c_str(), &info_mr, &n_mr) < 0) {
      nmono += 1; continue; }

    int32_t w = (int32_t) (info_mr[0] / unit);
    int32_t key = binarySearch<int32_t>(mr_keys,0,mr_keys.size()-1,w);
    if (key == -1) {key = 0;}
    key = mr_keys[key];

    auto ptr = cdf.find( key );
    (*ptr->second[ac])[it] += 1.;

    nVariant++;
    pre_pos = iv->pos;
  }

  notice("Finished Processing %d markers. Skipping %d filtered markers and %d uninformative markers", nVariant, nskip, nmono);


  for (int32_t ac = 2; ac <= maxac; ++ac) {
    std::string outf = out + "_" + std::to_string(ac) + ".txt";
    std::ofstream wf(outf.c_str(), std::ios::out);
    for (auto & v : mr_keys) {
      auto ptr = cdf.find(v);
      for (int32_t i = 0; i < nbins; ++i) {
        wf << v * unit << '\t' << (i+1)*binsize << '\t' << (*ptr->second[ac])[i] << '\n';
      }
    }
    wf.close();
  }
  return 0;
}










