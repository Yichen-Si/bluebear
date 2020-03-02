#include "cramore.h"
#include "utils.h"
#include "bcf_filter_arg.h"
#include "bcf_ordered_reader.h"
#include "bcf_ordered_writer.h"
#include "ibs0.h"
#include "rare_variant_ibs0.h"
#include <stdio.h>

// goal -- For a given list of positions, random sample major allele carriers
//         to get no-ibs0 length for null distribution
int32_t RandomPairIBS0(int32_t argc, char** argv) {
  std::string inVcf, inMap, inPos, chrom, reg;
  std::string out;
  int32_t min_hom_gts = 1;
  int32_t verbose = 1000;
  int32_t start, end;
  int32_t window_size = 500000;
  int32_t cst = -1, ced = -1;
  int32_t ck_len = 500000;
  int32_t bp_limit = 1000000;
  // double  cm_limit = -1.0;

  bcf_vfilter_arg vfilt;
  bcf_gfilter_arg gfilt;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Sites", NULL)
    LONG_STRING_PARAM("in-vcf",&inVcf, "Input VCF/BCF file")
    LONG_STRING_PARAM("position",&inPos, "Genomic position (sorted) to sample ibs0 length")
    LONG_STRING_PARAM("chrom",&chrom, "Chromosome to sample ibs0 length")
    LONG_STRING_PARAM("map",&inMap, "Genetic map for genetic distance")

    LONG_PARAM_GROUP("Variant Filtering Options", NULL)
    LONG_INT_PARAM("centromere-st",&cst, "Start position of centromere")
    LONG_INT_PARAM("centromere-ed",&ced, "End position of centromere")

    LONG_PARAM_GROUP("Additional Options", NULL)
    LONG_INT_PARAM("min-hom",&min_hom_gts, "Minimum number of homozygous genotypes to be counted for IBS0")
    LONG_INT_PARAM("bp-limit",&bp_limit, "Max distance to look for IBS0 (bp)")
    LONG_INT_PARAM("window-size",&window_size, "Length of the sliding window to procee (bp)")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output file name")
    LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")

  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  // sanity check of input arguments
  if ( inVcf.empty() || inPos.empty() || inMap.empty() || out.empty() ) {
    error("[E:%s:%d %s] --in-vcf, --position, --inMap, --out are required parameters",__FILE__,__LINE__,__FUNCTION__);
  }

  // Genetic map
  bp2cmMap pgmap(inMap, " ", "", cst, ced);
  notice("Read map file and detect centromere: %d-%d", pgmap.centromere_st, pgmap.centromere_ed);

  // Input vcf
  std::vector<GenomeInterval> intervals;
  parse_intervals(intervals, "", "");
  BCFOrderedReader odr(inVcf, intervals);
  int32_t nsamples = bcf_hdr_nsamples(odr.hdr);

  // Output
  FILE * wf;
  wf = fopen(out.c_str(), "w");

  // To randomly sample individual pairs
  std::random_device rd;
  std::mt19937 rng(rd());
  std::set<int32_t> chosen;
  std::vector<int32_t> pair(2);
  int32_t pos, leftibs0, rightibs0, lengthbp;
  double lengthcm;

  // Parse positions. Assume sorted
  std::ifstream rf (inPos);

  rf >> start;
  end = start + window_size;
  reg = chrom + ":" + std::to_string(start) + "-" + std::to_string(end);
  rf.clear();
  rf.seekg(0, std::ios::beg);

  // IBS0 reader
  IBS0lookup ibs0finder(inVcf, reg, pgmap, bp_limit, ck_len, min_hom_gts);
  notice("Set up ibs0finder for region %s", reg.c_str());

  int32_t ct = 0;
  while (rf >> pos) {
    if (ct % verbose == 0) {
      notice("Processing %d variants at %s:%d.", ct, chrom.c_str(), pos);
    }
    if (abs(pos-pgmap.centromere_st) < bp_limit || abs(pos-pgmap.centromere_ed) < bp_limit || pos < start) {continue;}
    NchooseK(nsamples, 2, chosen, rng, 0);
    std::copy(chosen.begin(),chosen.end(),pair.begin());
    chosen.clear();
    leftibs0 = ibs0finder.FindIBS0(pair[0],pair[1],pos,1);
    rightibs0= ibs0finder.FindIBS0(pair[0],pair[1],pos,0);
    if (leftibs0 <= 0 || rightibs0 <= 0) {
      lengthbp = bp_limit*2+window_size;
      lengthcm = pgmap.bpinterval2cm(pos-bp_limit,pos+bp_limit);
    } else {
      lengthbp = rightibs0 - leftibs0;
      lengthcm = pgmap.bpinterval2cm(leftibs0,rightibs0);
    }
    fprintf(wf, "%d\t%d\t%d\t%d\t%d\t%d\t%.5f\n",pos,pair[0],pair[1],leftibs0,rightibs0,lengthbp,lengthcm);
    ct++;
    if (pos > end) { // Update
      start = end + 1;
      end = start + window_size;
      while (pgmap.overlap_centro(start,end) > 0 && end < pgmap.maxpos) {
        start = end + 1;
        end = start + window_size;
      }
      if (start > pgmap.maxpos) {break;}
      reg = chrom + ":" + std::to_string(start) + "-" + std::to_string(end);
      ibs0finder.Update(reg);
      notice("Update ibs0finder for region %s", reg.c_str());
    }
  }
  notice("Finished processing %d variants", ct);

  pgmap.clear();
  rf.close();
  fclose(wf);

  return 0;
}

