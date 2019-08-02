#include "cramore.h"
#include "bcf_filter_arg.h"
#include "bcf_ordered_reader.h"

#include "pbwt_build.h"
#include "cm2bp.h"
#include "bp2cm.h"

int32_t IBS0Phase(int32_t argc, char** argv) {
  std::string inVcf, inMap, out, reg;
  int32_t verbose = 1000;
  int32_t min_variant = 1;
  int32_t min_hom_gts = 1;

  inMap = "/net/wonderland/home/ycsi/Data/Map/plink.chr20.GRCh38.map";
  cm2bpMap cbmap(inMap, " ");

  // Build prefix pbwt from begining. TODO: start from a close snapshot


  // Process & store \delta cM blocks for ibs0 look up

  //




  return 0;
}
















