#include "cramore.h"
#include "bcf_filter_arg.h"
#include "bcf_ordered_reader.h"
#include "compact_matrix.h"

#include "bp2cm.h"

int32_t cmdVcfIBS0full(int32_t argc, char** argv) {
  std::string inVcf, inMap;
  std::string out, bout, rout;
  int32_t min_hom_gts = 1;
  int32_t verbose = 1000;
  std::string reg;
  int32_t min_variant = 1;
  int32_t n_samples = 0;
  int32_t min_ibs0_bp = 500000;
  double  min_ibs0_cm  = 3.0;
  int32_t unit = 0;

  bcf_vfilter_arg vfilt;
  bcf_gfilter_arg gfilt;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Sites", NULL)
    LONG_STRING_PARAM("in-vcf",&inVcf, "Input VCF/BCF file")
    LONG_STRING_PARAM("region",&reg,"Genomic region to focus on")
    LONG_STRING_PARAM("map",&inMap, "Input linkage map file")

    LONG_PARAM_GROUP("Variant Filtering Options", NULL)
    LONG_MULTI_STRING_PARAM("apply-filter",&vfilt.required_filters, "Require at least one of the listed FILTER strings")
    LONG_STRING_PARAM("include-expr",&vfilt.include_expr, "Include sites for which expression is true")
    LONG_STRING_PARAM("exclude-expr",&vfilt.exclude_expr, "Exclude sites for which expression is true")

    LONG_PARAM_GROUP("Additional Options", NULL)
    LONG_INT_PARAM("min-hom",&min_hom_gts, "Minimum number of homozygous genotypes to be counted for IBS0")
    LONG_INT_PARAM("min-variant",&min_variant, "Minimum number of variants to present to have output file")
    LONG_INT_PARAM("min-bp",&min_ibs0_bp, "Minimum no-ibs0 in bp")
    LONG_DOUBLE_PARAM("min-cm",&min_ibs0_cm, "Minimum no-ibs0 in cm")
    LONG_INT_PARAM("unit",&unit, "Use cm (0) or bp (1) for lower bound")
    // For testing
    LONG_INT_PARAM("num-samples",&n_samples, "Number of samples to test")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output VCF file name")
    LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")


  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  // sanity check of input arguments
  if ( inVcf.empty() || out.empty() || inMap.empty() ) {
    error("[E:%s:%d %s] --in-vcf, --out, --map are required parameters",__FILE__,__LINE__,__FUNCTION__);
  }

  // physical -> genetic position
  bp2cmMap pgmap(inMap, " ");

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

  notice("Started Reading site information from VCF file, identifying %d samples", nsamples);

  if ( nsamples == 0 )
    error("FATAL ERROR: The VCF does not have any samples with genotypes");

  bitmatrix bmatRR(nsamples);
  bitmatrix bmatAA(nsamples);
  std::vector<int32_t> snoopy; // informative position

  int32_t* p_gt = NULL;
  int32_t  n_gt = 0;
  int32_t  nskip = 0, nmono = 0;
  uint8_t* gtRR = (uint8_t*)calloc(nsamples, sizeof(uint8_t));
  uint8_t* gtAA = (uint8_t*)calloc(nsamples, sizeof(uint8_t));
  int32_t *info_ac = NULL;
  int32_t *info_an = NULL;
  int32_t n_ac = 0, n_an = 0;

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

    // extract genotype and apply genotype level filter
    if ( bcf_get_genotypes(odr.hdr, iv, &p_gt, &n_gt) < 0 ) {
      error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);
    }
    if (bcf_get_info_int32(odr.hdr, iv, "AC", &info_ac, &n_ac) < 0) {continue;}
    if (bcf_get_info_int32(odr.hdr, iv, "AN", &info_an, &n_an) < 0) {continue;}
    if (info_ac[0] < 2 || info_an[0]-info_ac[0] < 2) {continue;}
    memset(gtRR, 0, nsamples);
    memset(gtAA, 0, nsamples);
    int32_t gcs[3] = {0,0,0};
    for(int32_t i=0; i < nsamples; ++i) {
      int32_t g1 = p_gt[2*i];
      int32_t g2 = p_gt[2*i+1];
      int32_t geno;
      if ( bcf_gt_is_missing(g1) || bcf_gt_is_missing(g2) ) {
	       //geno = 0;
      }
      else {
      	geno = ((bcf_gt_allele(g1) > 0) ? 1 : 0) + ((bcf_gt_allele(g2) > 0) ? 1 : 0);
      	if ( geno == 0 )    { gtRR[i] = 1; }
      	else if ( geno == 2 ) { gtAA[i] = 1; }
      	++gcs[geno];
      }
    }
    if (  ( gcs[0] < min_hom_gts ) || ( gcs[2] < min_hom_gts ) ) { ++nmono; }
    else {
      bmatRR.add_row_bytes(gtRR);
      bmatAA.add_row_bytes(gtAA);
      snoopy.push_back(iv->pos+1);
      ++nVariant;
    }
  }
  notice("Finished Processing %d markers across %d samples, Skipping %d filtered markers and %d uninformative markers", nVariant, nsamples, nskip, nmono);

  free(gtRR);
  free(gtAA);

  if ( nVariant < min_variant ) {
    notice("Observed only %d informative markers. Skipping IBD segment detection for this chunk...", nVariant);
    return 0;
  }

  bmatRR.transpose();
  bmatAA.transpose();

  // For testing
  if (n_samples > 0 && nsamples > n_samples)
    nsamples = n_samples;

  notice("Searching for potential IBD segments..");
  int32_t k = 0;

  bout = out + ".boundary";
  htsFile* lf = hts_open(bout.c_str(), "w");
  rout = out + ".noibs0";
  htsFile* wf = hts_open(rout.c_str(), "w");

  for(int32_t i=1; i < nsamples; ++i) {
    if (i % 500 == 0)
      notice("Processing pairs including the %dth individual...", i);
    uint8_t* iRR = bmatRR.get_row_bits(i);
    uint8_t* iAA = bmatAA.get_row_bits(i);
    for(int32_t j=0; j < i; ++j) {
      uint8_t* jRR = bmatRR.get_row_bits(j);
      uint8_t* jAA = bmatAA.get_row_bits(j);
      int32_t prepos = 0;
      std::vector<int32_t> bound(2,0);
      for(k=0; k < bmatRR.nbytes_col; ++k) {
        uint8_t byte = ( iRR[k] ^ jRR[k] ) & ( iAA[k] ^ jAA[k] );
        if (byte) {
          for (int32_t bit=0; bit<8 && k*8+bit<bmatRR.ncol; ++bit) {
            int32_t pt = k*8+bit;
          	if ( (byte >> (7-bit)) & 0x01 ) { // IBS0
              if ( prepos > 1 ) {
                double lcm = pgmap.bp2cm(snoopy[pt]) - pgmap.bp2cm(prepos);
                int32_t lbp = snoopy[pt] - prepos;
                if ((unit == 0 && lcm > min_ibs0_cm)||(unit != 0 && lbp > min_ibs0_bp)) {
                  hprintf(wf,"%d\t%d\t%d\t%.3f\t%s\t%s\n", prepos, snoopy[pt], lbp, lcm, odr.hdr->id[BCF_DT_SAMPLE][i].key, odr.hdr->id[BCF_DT_SAMPLE][j].key);
                }
              } else {
                bound[0] = snoopy[pt];
              }
              prepos = snoopy[pt];
          	}
          }
        }
      }
      bound[1] = prepos;
      // Output boundary for this pair of individual
      hprintf(lf,"%s\t%s\t%d\t%d\n",odr.hdr->id[BCF_DT_SAMPLE][i].key, odr.hdr->id[BCF_DT_SAMPLE][j].key, bound[0], bound[1]);
    }
  }

  hts_close(wf);
  hts_close(lf);

  notice("Finished searching for ibs0 segments");

  return 0;
}


int32_t cmdIBS0Bridge(int32_t argc, char** argv) {

  std::string path, pref, suff, out, outf, inMap;
  std::string filename, id1, id2;
  int32_t verbose = 0;
  int32_t chunksize = 1000000, start = 1, nchunk = 65;
  int32_t st, ed, ck=0, nline = 0;
  int32_t pos1, pos2;
  int32_t unit = 0;
  int32_t min_ibs0_bp = 500000;
  double  min_ibs0_cm  = 3.0;
  std::ifstream infile;
  std::vector<int32_t> ibs0vec;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input", NULL)
    LONG_STRING_PARAM("path",&path, "Input file path")
    LONG_STRING_PARAM("pref",&pref, "Input file prefix")
    LONG_STRING_PARAM("suff",&suff, "Input file suffix")
    LONG_STRING_PARAM("map",&inMap, "Input linkage map file")
    LONG_INT_PARAM("chunk-n",&nchunk, "Number of chunks (input files)")
    LONG_INT_PARAM("chunk-size",&chunksize, "Size of chunks in bp")

    LONG_PARAM_GROUP("Additional Options", NULL)
    LONG_INT_PARAM("min-bp",&min_ibs0_bp, "Minimum no-ibs0 in bp")
    LONG_DOUBLE_PARAM("min-cm",&min_ibs0_cm, "Minimum no-ibs0 in cm")
    LONG_INT_PARAM("unit",&unit, "Use cm (0) or bp (1) for lower bound")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output file prefix")
    LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")


  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  // sanity check of input arguments
  if ( path.empty() || pref.empty() || suff.empty() || out.empty() || inMap.empty() ) {
    error("[E:%s:%d %s] --path, --pref, --suff, --out, --map are required parameters",__FILE__,__LINE__,__FUNCTION__);
  }
  if (path.back() != '/' && pref.at(0) != '/')
    path+="/";

  bp2cmMap pgmap(inMap, " ");

  outf = out + ".noibs0";
  FILE * wf;
  wf = fopen(outf.c_str(), "w");

  st = ck*chunksize + start;
  ed = (ck+1)*chunksize;
  filename = path + pref + std::to_string(st) + "-" + std::to_string(ed) + suff;

  infile.open(filename, std::ifstream::in);
  while(infile >> id1 >> id2 >> pos1 >> pos2) {
    ibs0vec.push_back(pos2);
    nline++;
  }
  infile.close();

  notice("Detected %d lines.", nline);

  for (ck = 1; ck < nchunk; ++ck) {
    st = ck*chunksize + start;
    ed = (ck+1)*chunksize;
    filename = path + pref + std::to_string(st) + "-" + std::to_string(ed) + suff;
    infile.open(filename, std::ifstream::in);
    if (!infile.is_open()) {continue;}
    int32_t k = 0;
    while(infile >> id1 >> id2 >> pos1 >> pos2) {
      if (ibs0vec[k] == 0) {
        ibs0vec[k] = pos2;
        k++;
        continue;
      }
      if (pos1 != 0) {
        int32_t lbp = pos1 - ibs0vec[k];
        double  lcm = pgmap.bp2cm(pos1) - pgmap.bp2cm(ibs0vec[k]);
        if ((unit == 0 && lcm > min_ibs0_cm)||(unit != 0 && lbp > min_ibs0_bp)) {
          fprintf(wf, "%d\t%d\t%d\t%.3f\t%s\t%s\n", ibs0vec[k],pos1,lbp,lcm,id1.c_str(),id2.c_str());
        }
        ibs0vec[k] = pos2;
      }
      k++;
      if ( verbose > 0 && k % verbose == 0 )
        notice("Processed %d pairs in region %d-%d.", k, st, ed);
    }
    infile.close();
    notice("Processed %d chunks (at %d-%d).", ck, st, ed);
  }
  fclose(wf);

  return 0;

}
















