#include "cramore.h"
#include "bcf_filter_arg.h"
#include "bcf_ordered_reader.h"
#include "bcf_ordered_writer.h"
#include "hts_utils.h"

struct Faival {
  int64_t len, offset, line_blen, line_len;
  Faival(int64_t _l, int64_t _os, int64_t _bl, int64_t _ll) :
  len(_l), offset(_os), line_blen(_bl), line_len(_ll) {}
};

class FaReader {

  public:
    std::map<std::string, Faival*> lookup;

    FaReader(std::string &fnfai) {
      std::ifstream mfile(fnfai);
      std::string line;
      std::vector<std::string> words;
      int64_t len=0, offset=0, line_blen=0, line_len=0;
      const char* sep = "\t";
      if (mfile.is_open()) {
        while(std::getline(mfile, line)) {
          split(words, sep, line);
          len = std::stoll(words[1]);
          offset = std::stoll(words[2]);
          line_blen = std::stoll(words[3]);
          line_len = std::stoll(words[4]);
          Faival *tmp = new Faival(len, offset, line_blen, line_len);
          lookup[words[0]] = tmp;
        }
      } else {
        error("Unable to open .fai file");
      }
    }

    Faival* FindFaiVal(std::string &chrom) {
      if (lookup.find(chrom) != lookup.end()) {
        return lookup[chrom];
      } else {
        return NULL;
      }
    }

};

std::string fa_kmer(faidx_t *fai, Faival *val, int32_t pos, int32_t k) {
  int32_t st = pos - k/2;
  int32_t ret = bgzf_useek(fai->bgzf, val->offset + st / val->line_blen * val->line_len + st % val->line_blen, SEEK_SET);
  if ( ret<0 ) {
    error("Error: bgzf_useek failed.");
  }
  int32_t l = 0;
  char c;
  std::string seq = "";
  while ( (c=bgzf_getc(fai->bgzf))>=0 && l < k ) {
    if (isgraph(c)) {
      c = toupper(c);
      seq += c;
      l++;
    }
  }
  return(seq);
};


int32_t cmdVcfAddContexte(int32_t argc, char** argv) {
  std::string inVcf, inFa, chrom, out;
  int32_t kmer = 7;
  int64_t verbose = 50000;
  bool withGT = false;

  paramList pl;
  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input", NULL)
    LONG_STRING_PARAM("in-vcf",&inVcf, "Input VCF file to annotate")
    LONG_STRING_PARAM("in-fasta",&inFa, "Input Fasta file")
    LONG_STRING_PARAM("chr",&chrom,"Chromosome to process")

    LONG_PARAM_GROUP("Additional Options", NULL)
    LONG_INT_PARAM("kmer",&kmer, "Kmer length around each variant")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output VCF/BCF file")
    LONG_PARAM("with-GT", &withGT, "If GT field is included in input VCF/BCF and should be written to output VCF/BCF")

  END_LONG_PARAMS();
  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  if ( inVcf.empty() || inFa.empty() || chrom.empty() || out.empty() ) {
    error("[E:%s:%d %s] --in-vcf, --in-fasta, --chr, --out are required parameters",__FILE__,__LINE__,__FUNCTION__);
  }

  faidx_t *fai = fai_load(inFa.c_str());
  // Read .fai index
  std::string fnfai = inFa + ".fai";
  FaReader findx(fnfai);
  Faival *val = findx.FindFaiVal(chrom);
  if (val == NULL) {
    error("Cannot find fai index for this chromosome");
  }
  notice("Initialized reader of fasta file. Index info: %ld\t%ld\t%ld\t%ld\n ", val->len, val->offset, val->line_blen, val->line_len);


  // setup input & output BCF/VCF file
  std::vector<GenomeInterval> intervals;
  BCFOrderedReader* odr = new BCFOrderedReader(inVcf, intervals);
  bcf1_t* iv = bcf_init();

  BCFOrderedWriter odw(out.c_str(),0);
  if (withGT) {
    odw.set_hdr(odr->hdr);
  } else {
    bcf_hdr_t* hnull = bcf_hdr_subset(odr->hdr, 0, 0, 0);
    bcf_hdr_remove(hnull, BCF_HL_FMT, NULL);
    odw.set_hdr(hnull);
  }

  char buffer[65536];
  // check the existence of header and create one if needed
  if ( bcf_hdr_id2int(odr->hdr, BCF_DT_ID, "CpG") < 0 ) {
    sprintf(buffer,"##INFO=<ID=CpG,Number=0,Type=String,Description=\"If the variant is at a CpG site with the ref/major allele being the C\">\n");
    bcf_hdr_append(odw.hdr, buffer);
  }
  if ( bcf_hdr_id2int(odr->hdr, BCF_DT_ID, "CONTEXT") < 0 ) {
    sprintf(buffer,"##INFO=<ID=CONTEXT,Number=1,Type=String,Description=\"Nucleotide context from 5\'to 3\': %d bases motif around a variant\">\n", kmer);
    bcf_hdr_append(odw.hdr, buffer);
  }
  odw.write_hdr();

  // Process each variant and retrieve kmer context
  int32_t k = 0;
  int32_t mpt = kmer/2;
  int32_t *info_ac = NULL, *info_an = NULL;
  int32_t n_ac = 0, n_an = 0, ac, an;
  std::string ctx, cpg, refcpg, majcpg;
  for(; odr->read(iv); ++k) {
    if ( k % verbose == 0 )
      notice("Processing %d markers at %s:%d.", k, bcf_hdr_id2name(odr->hdr, iv->rid), iv->pos+1);
    if (bcf_get_info_int32(odr->hdr, iv, "AC", &info_ac, &n_ac) < 0) {continue;}
    if (bcf_get_info_int32(odr->hdr, iv, "AN", &info_an, &n_an) < 0) {continue;}

    bcf1_t* nv = bcf_dup(iv);
    ac = info_ac[0]; an = info_an[0];
    ctx = fa_kmer(fai, val, iv->pos, kmer);
    refcpg = "F";
    majcpg = "F";

    if (bcf_is_snp(iv)) { // SNP
      if (ctx[mpt]=='C' && ctx[mpt+1]=='G') { // Ref CpG
        refcpg = "T";
      }
      if (ac <= an/2) { // Ref == Major
        majcpg = refcpg;
      } else {
        if ((iv->d).allele[1][0]=='C' && ctx[mpt+1]=='G') { // Major CpG
          majcpg="T";
        }
      }
    }

    cpg = refcpg+majcpg;
    bcf_update_info_string(odw.hdr, nv, "CpG", cpg.c_str());
    bcf_update_info_string(odw.hdr, nv, "CONTEXT", ctx.c_str());
    odw.write(nv);
    bcf_destroy(nv);
  }

  odw.close();
  return 0;
}




