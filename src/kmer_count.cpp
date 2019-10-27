#include "cramore.h"
#include "hts_utils.h"

// Goal: count the occurence of kmers

int32_t KmerCount(int32_t argc, char** argv) {
  std::string inFa, out, outf;
  std::string chrom="chr20";
  int32_t kmer = 7;
  int64_t verbose = 5000000;

  paramList pl;
  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input", NULL)
    LONG_STRING_PARAM("in-fasta",&inFa, "Input Fasta file")
    LONG_STRING_PARAM("chr",&chrom,"Chromosome to process")

    LONG_PARAM_GROUP("Additional Options", NULL)
    LONG_INT_PARAM("kmer",&kmer, "k-mer")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output file prefix")

  END_LONG_PARAMS();
  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  faidx_t *fai = fai_load(inFa.c_str());

  // Read .fai index
  std::string fnfai = inFa + ".fai";
  std::ifstream mfile(fnfai);
  std::string line;
  std::vector<std::string> words;
  int64_t len=0, offset=0, line_blen=0, line_len=0;
  const char* sep = "\t";
  if (mfile.is_open()) {
    while(std::getline(mfile, line)) {
      split(words, sep, line);
      if (words[0] == chrom) {
        len = std::stoll(words[1]);
        offset = std::stoll(words[2]);
        line_blen = std::stoll(words[3]);
        line_len = std::stoll(words[4]);
        break;
      }
    }
  }

  if (len == 0) {
    error("Cannot read index.");
  }
  std::cout << "Initialized reader\n";
  std::cout << len << '\t' << offset << '\t' << line_blen << '\t' << line_len << '\n';

  std::string seq = "";
  std::set<char> alphabet {'A','T','C','G'};
  std::map<std::string, int32_t> kmerct;
  char c;
  int32_t l = 0, p_beg_i = 0;

  // Retrive the first kmer
  int ret = bgzf_useek(fai->bgzf, offset + p_beg_i / line_blen * line_len + p_beg_i % line_blen, SEEK_SET);
  std::cout << ret << '\n';

  if ( ret<0 ) {
    error("Error: bgzf_useek failed.");
  }
  l = 0;
  while ( (c=bgzf_getc(fai->bgzf))>=0 && l < kmer ) {
    if (isgraph(c)) {
      c = toupper(c);
      if (alphabet.find(c) == alphabet.end()) {
        l = 0;
        seq = "";
        continue;
      }
      seq += c;
      l++;
    }
  }
  kmerct[seq]++;

  int32_t verb = 0;
  // Read one base at a time
  while ((c=bgzf_getc(fai->bgzf))>=0) {
    if (isgraph(c)) {
      if (c == '>') {
        break;
      }
      c = toupper(c);
      if (alphabet.find(c) == alphabet.end()) {
        seq = "";
        continue;
      }
      if ((int32_t) seq.length() < kmer) {
        seq += c;
      }
      else if ((int32_t) seq.length() == kmer) {
        seq = seq.substr(1) + c;
        kmerct[seq]++;
      } else {
        notice("Potential incorrect shift in reading by kmer.");
        seq = "";
      }
      verb++;
    if ( verb % verbose == 0 )
      notice("Processed %d Mb", verb/1000000);
    }
  }

  outf = out + "_" + chrom + "_kmer_" + std::to_string(kmer) + "_count.txt";
  htsFile *wf = hts_open(outf.c_str(), "w");
  for (auto const& v : kmerct) {
    hprintf(wf, "%s\t%d\t%s\t%d\n", chrom.c_str(), len, v.first.c_str(), v.second);
  }

  hts_close(wf);
  return 0;
}


