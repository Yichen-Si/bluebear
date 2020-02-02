#include "cramore.h"
#include "bcf_filter_arg.h"
#include "bcf_ordered_reader.h"
#include "bcf_ordered_writer.h"
#include "hts_utils.h"
#include "bp2cm.h"
#include "ibs0.h"

int32_t test(int32_t argc, char** argv) {

  std::string inVcf, reg, inMap, chrom;
  int32_t pos = 1000000, st = 15000000, ed = 16000000;
  // std::string inVcf, inFa, chrom, out;
  // int32_t kmer = 7;
  // int64_t verbose = 50000;
  // bool withGT = false;

  paramList pl;
  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input", NULL)
    LONG_STRING_PARAM("map",&inMap, "Input map file")
    LONG_INT_PARAM("pos",&pos, "POS")
    LONG_INT_PARAM("st",&st, "ST")
    LONG_INT_PARAM("ed",&ed, "ED")

    // LONG_PARAM_GROUP("Additional Options", NULL)
    LONG_STRING_PARAM("in-vcf",&inVcf, "Input VCF file to annotate")
    LONG_STRING_PARAM("region",&reg, "Input VCF file to annotate")
  //   LONG_STRING_PARAM("in-fasta",&inFa, "Input Fasta file")
    LONG_STRING_PARAM("chr",&chrom,"Chromosome to process")

  //   LONG_PARAM_GROUP("Additional Options", NULL)
  //   LONG_INT_PARAM("kmer",&kmer, "Kmer length around each variant")

  //   LONG_PARAM_GROUP("Output Options", NULL)
  //   LONG_STRING_PARAM("out", &out, "Output VCF/BCF file")
  //   LONG_PARAM("with-GT", &withGT, "If GT field is included in input VCF/BCF and should be written to output VCF/BCF")

  END_LONG_PARAMS();
  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  bp2cmMap pgmap(inMap, " ", chrom);
  IBS0lookup ibs0finder(inVcf, reg, pgmap, 5000, 5000, 1);
  std::cout << "Tested constructor\n";
  int32_t win = 50000;
  int32_t i = 0;
  for (int32_t i = 0; i < 2; ++i) {
    st = ed + 1;
    ed = st + win;
    reg = "chr" + chrom + ":" + std::to_string(st) + "-" + std::to_string(ed);
    ibs0finder.Update(reg);
    std::cout << ibs0finder.bmatRR_que.size() << ' ' << ibs0finder.bmatAA_que.size() << ' ' << ibs0finder.start_que[0] << '-' << ibs0finder.posvec_que.back()->back() << '\n';
  }
  std::cout << "Tested update\n";
  reg = "chr" + chrom + ":" + std::to_string(pos) + "-" + std::to_string(pos+win);
  ibs0finder.Update_Fixed(reg);
  std::cout << "Tested update_fixed\n";

  return 0;
}




