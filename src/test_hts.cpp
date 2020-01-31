#include "cramore.h"
#include "bcf_filter_arg.h"
#include "bcf_ordered_reader.h"
#include "bcf_ordered_writer.h"
#include "hts_utils.h"
#include "gmap.h"

// struct Faival {
//   int64_t len, offset, line_blen, line_len;
//   Faival(int64_t _l, int64_t _os, int64_t _bl, int64_t _ll) :
//   len(_l), offset(_os), line_blen(_bl), line_len(_ll) {}
// };

// class FaReader {

//   public:
//     std::map<std::string, Faival*> lookup;

//     FaReader(std::string &fnfai) {
//       std::ifstream mfile(fnfai);
//       std::string line;
//       std::vector<std::string> words;
//       int64_t len=0, offset=0, line_blen=0, line_len=0;
//       const char* sep = "\t";
//       if (mfile.is_open()) {
//         while(std::getline(mfile, line)) {
//           split(words, sep, line);
//           len = std::stoll(words[1]);
//           offset = std::stoll(words[2]);
//           line_blen = std::stoll(words[3]);
//           line_len = std::stoll(words[4]);
//           Faival *tmp = new Faival(len, offset, line_blen, line_len);
//           lookup[words[0]] = tmp;
//         }
//       } else {
//         error("Unable to open .fai file");
//       }
//     }

//     Faival* FindFaiVal(std::string &chrom) {
//       if (lookup.find(chrom) != lookup.end()) {
//         return lookup[chrom];
//       } else {
//         return NULL;
//       }
//     }

// };

// std::string fa_kmer(faidx_t *fai, Faival *val, int32_t pos, int32_t k) {
//   int32_t st = pos - k/2;
//   int32_t ret = bgzf_useek(fai->bgzf, val->offset + st / val->line_blen * val->line_len + st % val->line_blen, SEEK_SET);
//   if ( ret<0 ) {
//     error("Error: bgzf_useek failed.");
//   }
//   int32_t l = 0;
//   char c;
//   std::string seq = "";
//   while ( (c=bgzf_getc(fai->bgzf))>=0 && l < k ) {
//     if (isgraph(c)) {
//       c = toupper(c);
//       seq += c;
//       l++;
//     }
//   }
//   return(seq);
// };


int32_t test(int32_t argc, char** argv) {

  std::string inMap;
  int32_t pos = 1000000, st = 15000000, ed = 16000000;
  // std::string inVcf, inFa, chrom, out;
  // int32_t kmer = 7;
  // int64_t verbose = 50000;
  // bool withGT = false;

  paramList pl;
  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input", NULL)
    LONG_STRING_PARAM("in-map",&inMap, "Input map file")
    LONG_INT_PARAM("pos",&pos, "POS")
    LONG_INT_PARAM("st",&st, "ST")
    LONG_INT_PARAM("ed",&ed, "ED")

    // LONG_PARAM_GROUP("Additional Options", NULL)
    // LONG_STRING_PARAM("in-vcf",&inVcf, "Input VCF file to annotate")
  //   LONG_STRING_PARAM("in-fasta",&inFa, "Input Fasta file")
  //   LONG_STRING_PARAM("chr",&chrom,"Chromosome to process")

  //   LONG_PARAM_GROUP("Additional Options", NULL)
  //   LONG_INT_PARAM("kmer",&kmer, "Kmer length around each variant")

  //   LONG_PARAM_GROUP("Output Options", NULL)
  //   LONG_STRING_PARAM("out", &out, "Output VCF/BCF file")
  //   LONG_PARAM("with-GT", &withGT, "If GT field is included in input VCF/BCF and should be written to output VCF/BCF")

  END_LONG_PARAMS();
  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  gMap pgmap(inMap, " ", "19");
  std::cout << pgmap.minpos << '\t' << pgmap.maxpos << '\n';
  std::cout << pgmap.maxcm  << '\t' << pgmap.ctrcm  << '\n';
  std::cout << pgmap.centromere_st << '\t' << pgmap.centromere_ed << '\n';
  std::cout << pgmap.bp2cm(pos) << '\n';
  std::cout << pgmap.bpinterval2cm(st,ed) << '\n';

  return 0;
}




