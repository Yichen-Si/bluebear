#include "cramore.h"
#include "hts_utils.h"
#include "utils.h"
#include "fa_reader.h"
#include "seq_basics.h"
#include "fstream"
#include "sstream"

int32_t test(int32_t argc, char** argv) {
    // std::string inFa, out, chrom;
    // int32_t start = 0, end = -1;
    // int32_t kmer = 9, kmer0 = 3;
    // int32_t verbose = 1000000;
    // int32_t step = 50000;
    //
    // paramList pl;
    // BEGIN_LONG_PARAMS(longParameters)
    //     LONG_PARAM_GROUP("Input", NULL)
    //     LONG_STRING_PARAM("in-fasta",&inFa, "Input Fasta file")
    //     LONG_STRING_PARAM("chr",&chrom,"Chromosome to process. Should exist in the fasta file")
    //
    //     LONG_PARAM_GROUP("Additional Options", NULL)
    //     LONG_INT_PARAM("kmer",&kmer, "k-mer")
    //     LONG_INT_PARAM("chunk-size",&step, "Process by chunk")
    //     LONG_INT_PARAM("start",&start, "Start position")
    //     LONG_INT_PARAM("end",&end, "End position")
    //
    //     LONG_PARAM_GROUP("Output Options", NULL)
    //     LONG_STRING_PARAM("out", &out, "Output file")
    //     LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (bp)")
    //
    // END_LONG_PARAMS();
    // pl.Add(new longParams("Available Options", longParameters));
    // pl.Read(argc, argv);
    // pl.Status();
    //
    // std::string inFai = inFa + ".fai";
    // if (!file_exists(inFa) || !file_exists(inFai)) {
    // error("[E:%s:%d %s] input --in-fasta or its .fai index file does not exist",__FILE__,__LINE__,__FUNCTION__);
    // }
    //
    // FaReader findx(inFa);
    // if (!findx.FindFaiVal(chrom)) {
    //     error("Cannot find fai index for this chromosome");
    // }
    // notice("Initialized reader of fasta file. Index info: %ld\t%ld\t%ld\t%ld\n ", findx.val->len, findx.val->offset, findx.val->line_blen, findx.val->line_len);
    //
    // htsFile *wf = hts_open(out.c_str(), "w");
	// notice("Output saves to \n%s", out.c_str());
    //
    // std::string word0, word;
	// std::set<char> focal_base_set {'A', 'C', 'N'};
	// std::map<char, char> bpair = { {'A','T'}, {'T','A'}, {'C','G'}, {'G','C'}, {'N','N'} };
    //
    // int32_t verb = 0, pos = 0, nvb = 0;
	// int32_t st = start;
	// int32_t ed = st + step - 1;
    // int32_t kmer_pad_left = (int32_t) (kmer/2);
    // int32_t kmer_pad_right= kmer - kmer_pad_left - 1;
    // if (st < kmer_pad_left) {
    //     st = kmer_pad_left;
    // }
    //
    // st=60000;
    // ed=70000;
    // std::string seq = findx.fa_get_seq(st, ed);
    // for (int32_t i = 0; i < ed-st; ++i ) {
    //     if (i % 50 == 0) {
    //         std::cout << '\n';
    //     }
    //     std::cout << seq[i];
    // }
  return 0;
}
