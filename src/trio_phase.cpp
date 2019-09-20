#include "bcf_filter_arg.h"
#include "cramore.h"
#include "bcf_ordered_reader.h"
#include "compact_matrix.h"

#include <algorithm>
#include "bp2cm.h"

// goal -- Use parents' genotype to detect switch errors in phased child

int32_t trioSwitchDetect(int32_t argc, char** argv) {
  std::string childVCF, parentsVCF;
  std::string out, outf;
  int32_t verbose = 1000;
  std::string reg;
  int32_t min_variant = 1;
  int32_t n_samples = 0;

  bcf_vfilter_arg vfilt;
  bcf_gfilter_arg gfilt;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Sites", NULL)
    LONG_STRING_PARAM("child-vcf",&childVCF, "Input VCF/BCF file of phased child genotype")
    LONG_STRING_PARAM("parents-vcf",&parentsVCF, "Input VCF/BCF file of parents genotype")
    LONG_STRING_PARAM("region",&reg,"Genomic region to focus on")

    LONG_PARAM_GROUP("Variant Filtering Options", NULL)
    LONG_MULTI_STRING_PARAM("apply-filter",&vfilt.required_filters, "Require at least one of the listed FILTER strings")
    LONG_STRING_PARAM("include-expr",&vfilt.include_expr, "Include sites for which expression is true")
    LONG_STRING_PARAM("exclude-expr",&vfilt.exclude_expr, "Exclude sites for which expression is true")

    LONG_PARAM_GROUP("Additional Options", NULL)
    LONG_INT_PARAM("min-variant",&min_variant, "Minimum number of variants to present to have output file")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output VCF file name")
    LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")

    // For testing
    LONG_INT_PARAM("num-samples",&n_samples, "Number of samples to test")

  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  // sanity check of input arguments
  if ( childVCF.empty() || out.empty() || childVCF.empty() ) {
    error("[E:%s:%d %s] --child-vcf, --parents-vcf and --out are required parameters",__FILE__,__LINE__,__FUNCTION__);
  }

  // bcf reader
  std::vector<GenomeInterval> intervals;
  if ( !reg.empty() ) {
    parse_intervals(intervals, "", reg);
  }
  BCFOrderedReader odr_c(childVCF, intervals);
  bcf1_t* iv_c = bcf_init();
  BCFOrderedReader odr_p(parentsVCF, intervals);
  bcf1_t* iv_p = bcf_init();

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
  if ( filter_logic != 0 ) {
    filter_init(odr_c.hdr, filter_str.c_str());
  }

  // handle --apply-filtrs
  std::vector<int32_t> req_flt_ids;
  if ( !vfilt.required_filters.empty() ) {
    for(int32_t i=0; i < (int32_t)vfilt.required_filters.size(); ++i) {
      req_flt_ids.push_back(bcf_hdr_id2int(odr_c.hdr, BCF_DT_ID, vfilt.required_filters[i].c_str()));
    }
  }

  int32_t nsamples = bcf_hdr_nsamples(odr_c.hdr);


  notice("Started Reading site information from VCF file, identifying %d samples", nsamples);
  if ( nsamples == 0 )
    error("FATAL ERROR: The VCF does not have any samples with genotypes");

  bitmatrix bmatAA(2*nsamples);  // Parents AA
  bitmatrix bmatRR(2*nsamples);  // Parents AA
  bitmatrix bmatHap(2*nsamples); // Child haplotype

  std::vector<int32_t> pos_c, pos_p; // informative position

  int32_t* p_gt = NULL;
  int32_t n_gt = 0;
  uint8_t* gtAA  = (uint8_t*)calloc(2*nsamples, sizeof(uint8_t));
  uint8_t* gtRR  = (uint8_t*)calloc(2*nsamples, sizeof(uint8_t));
  uint8_t* gthap = (uint8_t*)calloc(2*nsamples, sizeof(uint8_t));
  int32_t nVariant = 0;

  for(int32_t k=0; odr_c.read(iv_c); ++k) {  // read haps of child
    // periodic message to user
    if ( k % verbose == 0 )
      notice("Processing child vcf: %d markers at %s:%d", k, bcf_hdr_id2name(odr_c.hdr, iv_c->rid), iv_c->pos+1);
    // unpack FILTER column
    bcf_unpack(iv_c, BCF_UN_FLT);
    // check --apply-filters
    bool has_filter = req_flt_ids.empty() ? true : false;
    if ( ! has_filter ) {
      //notice("%d %d", iv->d.n_flt, (int32_t)req_flt_ids.size());
      for(int32_t i=0; i < iv_c->d.n_flt; ++i) {
      	for(int32_t j=0; j < (int32_t)req_flt_ids.size(); ++j) {
      	  if ( req_flt_ids[j] == iv_c->d.flt[i] )
      	    has_filter = true;
      	}
      }
    }
    if ( ! has_filter ) { continue; }
    // check filter logic
    if ( filt != NULL ) {
      int32_t ret = filter_test(filt, iv_c, NULL);
      if ( filter_logic == FLT_INCLUDE ) { if ( !ret)  has_filter = false; }
      else if ( ret ) { has_filter = false; }
    }
    if ( ! has_filter ) { continue; }

    // extract genotype and apply genotype level filter
    if ( bcf_get_genotypes(odr_c.hdr, iv_c, &p_gt, &n_gt) < 0 ) {
      error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr_c.hdr, iv_c->rid), iv_c->pos+1);
    }
    memset(gthap, 0, 2*nsamples);
    int32_t gcs = 0;
    for(int32_t i=0; i < nsamples; ++i) {
      int32_t g1 = p_gt[2*i];
      int32_t g2 = p_gt[2*i+1];
      if ( bcf_gt_is_missing(g1) || bcf_gt_is_missing(g2)) {
        // Default 0
      }
      else {
        gthap[2*i] = (bcf_gt_allele(g1) > 0) ? 1 : 0;
        gthap[2*i+1] = (bcf_gt_allele(g2) > 0) ? 1 : 0;
        if (gthap[2*i] + gthap[2*i+1] == 1) gcs++;
      }
    }
    if (gcs == 0) {continue;}
    nVariant++;
    pos_c.push_back(iv_c->pos+1);
    bmatHap.add_row_bytes(gthap);
  }

  int32_t k = 0;
  while (k<nVariant) {  // read gt of parents
    memset(gtRR, 0, 2*nsamples);
    memset(gtAA, 0, 2*nsamples);
    if (!odr_p.read(iv_p)) {
      while (k < nVariant) {
        bmatRR.add_row_bytes(gtRR);
        bmatAA.add_row_bytes(gtAA);
        k++;
      }
      break;
    }
    while (k < nVariant && iv_p->pos+1 > pos_c[k]) {
      bmatRR.add_row_bytes(gtRR);
      bmatAA.add_row_bytes(gtAA);
      k++;
    }
    if (k == nVariant) {
      break;
    }
    bool flag = 0;
    while (iv_p->pos+1 < pos_c[k]) {
      if (!odr_p.read(iv_p)) {
        flag = 1;
        break;
      }
    }
    if (flag)
      continue;

    // periodic message to user
    if ( k % verbose == 0 )
      notice("Processing parents vcf: %d markers at %s:%d", k, bcf_hdr_id2name(odr_p.hdr, iv_p->rid), iv_p->pos+1);
    // extract genotype and apply genotype level filter
    if ( bcf_get_genotypes(odr_p.hdr, iv_p, &p_gt, &n_gt) < 0 ) {
      error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr_p.hdr, iv_p->rid), iv_p->pos+1);
    }
    for(int32_t i=0; i < nsamples; ++i) {
      int32_t pi = 2*i;
      int32_t mi = 2*i+1;
      int32_t g11 = p_gt[2*pi], g12 = p_gt[2*pi+1];
      int32_t g21 = p_gt[2*mi], g22 = p_gt[2*mi+1];
      if ( bcf_gt_is_missing(g11) || bcf_gt_is_missing(g12) || bcf_gt_is_missing(g21) || bcf_gt_is_missing(g22) ) {
        // Default 0
      }
      else {
        gtRR[2*i]   = ((bcf_gt_allele(g11) > 0) ? 0 : 1) && ((bcf_gt_allele(g12) > 0) ? 0 : 1);
        gtRR[2*i+1] = ((bcf_gt_allele(g21) > 0) ? 0 : 1) && ((bcf_gt_allele(g22) > 0) ? 0 : 1);
        gtAA[2*i]   = ((bcf_gt_allele(g11) > 0) ? 1 : 0) && ((bcf_gt_allele(g12) > 0) ? 1 : 0);
        gtAA[2*i+1] = ((bcf_gt_allele(g21) > 0) ? 1 : 0) && ((bcf_gt_allele(g22) > 0) ? 1 : 0);
      }
    }
    bmatRR.add_row_bytes(gtRR);
    bmatAA.add_row_bytes(gtAA);
    k++;
  }

  notice("Finished Processing %d/%d markers across %d samples.", nVariant, bmatRR.nrow, nsamples);

  free(gtRR);
  free(gtAA);
  free(gthap);

  if ( nVariant < min_variant ) {
    notice("Observed only %d informative markers. Skipping IBD segment detection for this chunk...", nVariant);
    return 0;
  }

  bmatRR.transpose();
  bmatAA.transpose();
  bmatHap.transpose();

  // For testing
  if (n_samples > 0 && nsamples > n_samples)
    nsamples = n_samples;

  outf = out;
  if ( !reg.empty() ) {
    outf += ("_" + reg);
  }
  outf += "_switch.pos";
  htsFile* wf = hts_open(outf.c_str(), "w");
  for(int32_t i=0; i < nsamples; ++i) {
    int32_t pre_config = -1, config = 0;
    int32_t st = 0, ed = 0;
    std::vector<int32_t> fvec;
    if (i % 200 == 0)
      notice("Processing the %dth individual...", i);
    uint8_t* pRR = bmatRR.get_row_bits(2*i);
    uint8_t* pAA = bmatAA.get_row_bits(2*i);
    uint8_t* mRR = bmatRR.get_row_bits(2*i+1);
    uint8_t* mAA = bmatAA.get_row_bits(2*i+1);
    uint8_t* c1 = bmatHap.get_row_bits(2*i);
    uint8_t* c2 = bmatHap.get_row_bits(2*i+1);
    for(k=0; k < bmatHap.nbytes_col; ++k) {
      uint8_t byte = (c1[k] ^ c2[k]) & ((~c1[k]) ^ (~c2[k])) & ((pAA[k] ^ mAA[k]) ^ (pRR[k] ^ mRR[k]));
      if (!byte) {continue;}
      config = (c1[k] & pAA[k]) ^ (c2[k] & mAA[k]) ^ (c2[k] & pRR[k]) ^ (c1[k] & mRR[k]);
      for (int32_t bit=0; bit<8 && k*8+bit<bmatHap.ncol; ++bit) {
        int32_t pt = k*8+bit;
        if ((byte >> (7-bit)) & 0x01) { // Informative
          if (pre_config < 0) {
            pre_config = (config>>(7-bit) & 0x01);
            st = pre_config;
          } else {
            if ( (config>>(7-bit) & 0x01) != pre_config ) {
              fvec.push_back(pos_c[pt]);
              pre_config = 1 - pre_config;
            }
          }
        }
      }
    }
    ed = pre_config;
    std::string nline = "";
    for (auto & v : fvec) {
      nline += std::to_string(v) + ",";
    }
    std::string trio_id(odr_c.hdr->id[BCF_DT_SAMPLE][i].key);
    trio_id += '\t';
    trio_id.append(odr_p.hdr->id[BCF_DT_SAMPLE][2*i].key);
    trio_id += '\t';
    trio_id.append(odr_p.hdr->id[BCF_DT_SAMPLE][2*i+1].key);
    hprintf(wf,"%s\t%d\t%s\t%d\n", trio_id.c_str(), st, nline.c_str(), ed);
  }
  hts_close(wf);
  notice("Finished searching for switch errors\n");

  return 0;
}


















// To find potential blips using trio

int32_t trioSwitchDetect_onepass(int32_t argc, char** argv) {
  std::string childVCF, parentsVCF;
  std::string out, outf, outb;
  int32_t verbose = 100000;
  std::string reg;
  int32_t min_variant = 1;
  int32_t n_samples = 0;

  bcf_vfilter_arg vfilt;
  bcf_gfilter_arg gfilt;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Sites", NULL)
    LONG_STRING_PARAM("child-vcf",&childVCF, "Input VCF/BCF file of phased child genotype")
    LONG_STRING_PARAM("parents-vcf",&parentsVCF, "Input VCF/BCF file of parents genotype")
    LONG_STRING_PARAM("region",&reg,"Genomic region to focus on")

    LONG_PARAM_GROUP("Variant Filtering Options", NULL)
    LONG_MULTI_STRING_PARAM("apply-filter",&vfilt.required_filters, "Require at least one of the listed FILTER strings")
    LONG_STRING_PARAM("include-expr",&vfilt.include_expr, "Include sites for which expression is true")
    LONG_STRING_PARAM("exclude-expr",&vfilt.exclude_expr, "Exclude sites for which expression is true")

    LONG_PARAM_GROUP("Additional Options", NULL)
    LONG_INT_PARAM("min-variant",&min_variant, "Minimum number of variants to present to have output file")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output VCF file name")
    LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")

    // For testing
    LONG_INT_PARAM("num-samples",&n_samples, "Number of samples to test")

  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  // sanity check of input arguments
  if ( childVCF.empty() || out.empty() || childVCF.empty() ) {
    error("[E:%s:%d %s] --child-vcf, --parents-vcf and --out are required parameters",__FILE__,__LINE__,__FUNCTION__);
  }

  // bcf reader
  std::vector<GenomeInterval> intervals;
  if ( !reg.empty() ) {
    parse_intervals(intervals, "", reg);
  }
  BCFOrderedReader odr_c(childVCF, intervals);
  bcf1_t* iv_c = bcf_init();
  BCFOrderedReader odr_p(parentsVCF, intervals);
  bcf1_t* iv_p = bcf_init();

  // Output
  outf = out;
  if ( !reg.empty() ) {
    outf += ("_" + reg);
  }
  outb = outf;
  outf += "_switch.pos";
  outb += "_blip.pos";
  std::ofstream wf,wfb;
  wf.open(outf, std::ios::trunc);
  wfb.open(outb, std::ios::trunc);

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
  if ( filter_logic != 0 ) {
    filter_init(odr_c.hdr, filter_str.c_str());
  }

  // handle --apply-filtrs
  std::vector<int32_t> req_flt_ids;
  if ( !vfilt.required_filters.empty() ) {
    for(int32_t i=0; i < (int32_t)vfilt.required_filters.size(); ++i) {
      req_flt_ids.push_back(bcf_hdr_id2int(odr_c.hdr, BCF_DT_ID, vfilt.required_filters[i].c_str()));
    }
  }

  int32_t nsamples = bcf_hdr_nsamples(odr_c.hdr);


  notice("Started Reading site information from VCF file, identifying %d samples", nsamples);
  if ( nsamples == 0 )
    error("FATAL ERROR: The VCF does not have any samples with genotypes");

  std::vector<bool> preconfig(nsamples,0);
  std::vector<int32_t> prehet(nsamples,1);
  std::vector<std::vector<int32_t> > switchpos(nsamples, std::vector<int32_t>(1,0));
  int32_t* p_gt = NULL;
  int32_t  n_gt = 0;
  uint8_t* gthap = (uint8_t*)calloc(2*nsamples, sizeof(uint8_t));
  int32_t  nVariant = 0;
  odr_p.read(iv_p);

  for(int32_t k=0; odr_c.read(iv_c); ++k) {
    // periodic message to user
    if ( k % verbose == 0 )
      notice("Processing child vcf: %d markers at %s:%d", k, bcf_hdr_id2name(odr_c.hdr, iv_c->rid), iv_c->pos+1);
    // unpack FILTER column
    bcf_unpack(iv_c, BCF_UN_FLT);
    // check --apply-filters
    bool has_filter = req_flt_ids.empty() ? true : false;
    if ( ! has_filter ) {
      //notice("%d %d", iv_c->d.n_flt, (int32_t)req_flt_ids.size());
      for(int32_t i=0; i < iv_c->d.n_flt; ++i) {
        for(int32_t j=0; j < (int32_t)req_flt_ids.size(); ++j) {
          if ( req_flt_ids[j] == iv_c->d.flt[i] )
            has_filter = true;
        }
      }
    }
    if ( ! has_filter ) { continue; }
    // check filter logic
    if ( filt != NULL ) {
      int32_t ret = filter_test(filt, iv_c, NULL);
      if ( filter_logic == FLT_INCLUDE ) { if ( !ret)  has_filter = false; }
      else if ( ret ) { has_filter = false; }
    }
    if ( ! has_filter ) { continue; }

    // extract genotype and apply genotype level filter
    if ( bcf_get_genotypes(odr_c.hdr, iv_c, &p_gt, &n_gt) < 0 ) {
      error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr_c.hdr, iv_c->rid), iv_c->pos+1);
    }
    memset(gthap, 0, 2*nsamples);
    int32_t gcs = 0;
    for(int32_t i=0; i < nsamples; ++i) {
      int32_t g1 = p_gt[2*i];
      int32_t g2 = p_gt[2*i+1];
      if ( bcf_gt_is_missing(g1) || bcf_gt_is_missing(g2)) {
        // If any allele is missing, don't make decision
      }
      else {
        gthap[2*i] = (bcf_gt_allele(g1) > 0) ? 1 : 0;
        gthap[2*i+1] = (bcf_gt_allele(g2) > 0) ? 1 : 0;
        if (gthap[2*i] + gthap[2*i+1] == 1) gcs++;
      }
    }
    if (gcs == 0) {continue;}

    // Move the reader of parents genotype to this position
    while (iv_p->pos < iv_c->pos) {odr_p.read(iv_p);}
    if (iv_p->pos != iv_c->pos) {continue;}

    // Extract prarents genotypes
    if ( bcf_get_genotypes(odr_p.hdr, iv_p, &p_gt, &n_gt) < 0 ) {
      error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr_p.hdr, iv_p->rid), iv_p->pos+1);
    }
    for(int32_t i=0; i < nsamples; ++i) {
      bool c1 = gthap[2*i], c2 = gthap[2*i+1];
      if (c1 == c2) {continue;}
      int32_t pi = 2*i, mi = 2*i+1;
      int32_t g11 = p_gt[2*pi], g12 = p_gt[2*pi+1];
      int32_t g21 = p_gt[2*mi], g22 = p_gt[2*mi+1];
      bool gtRR1, gtRR2, gtAA1, gtAA2;
      if ( bcf_gt_is_missing(g11) || bcf_gt_is_missing(g12) || bcf_gt_is_missing(g21) || bcf_gt_is_missing(g22) ) {
        // If any is missing, don't make decision at this position
        continue;
      }
      else {
        gtRR1 = ((bcf_gt_allele(g11) > 0) ? 0 : 1) && ((bcf_gt_allele(g12) > 0) ? 0 : 1);
        gtRR2 = ((bcf_gt_allele(g21) > 0) ? 0 : 1) && ((bcf_gt_allele(g22) > 0) ? 0 : 1);
        gtAA1 = ((bcf_gt_allele(g11) > 0) ? 1 : 0) && ((bcf_gt_allele(g12) > 0) ? 1 : 0);
        gtAA2 = ((bcf_gt_allele(g21) > 0) ? 1 : 0) && ((bcf_gt_allele(g22) > 0) ? 1 : 0);
      }
      if (!(gtRR1 || gtRR2 || gtAA1 || gtAA2)) {continue;}
      if (gtRR1 == gtRR2 && gtAA1 == gtAA2) {continue;}
      // Informative
      bool config = ((c1 && gtAA1) || (c2 && gtAA2) || (c2 && gtRR1) || (c1 && gtRR2));
      if (preconfig[i]!=config) {
// if (i == 0) {
//   std::cout << prehet[i] << '\t' << preconfig[i] << '\n';
//   std::cout << iv_c->pos << '\t' << c1 << c2 << '\t' << bcf_gt_allele(g11) << bcf_gt_allele(g12) << '\t' << bcf_gt_allele(g21) << bcf_gt_allele(g22) << '\n';
// }
        if (prehet[i] == switchpos[i].back()) { // Blip
          wfb << iv_c->pos << '\t' << prehet[i] << '\t' << i << '\n';
        }
        switchpos[i].push_back(iv_c->pos);
        preconfig[i] = config;
      }
      prehet[i] = iv_c->pos;
    }
    nVariant++;
  }
  free(gthap);
  wfb.close();

  for (int32_t i = 0; i < nsamples; ++i) {
    std::stringstream recline;
    recline << i << '\t' << switchpos[i].size() << '\t';
    uint32_t j = 1;
    while (j < switchpos[i].size()) {
      recline << switchpos[i][j] << ',';
      j++;
    }
    recline.seekp(-1, std::ios_base::end);
    recline << '\n';
    wf << recline.str();
  }
  wf.close();
  notice("Finished searching for switch errors\n");

  return 0;
}


