#include "cramore.h"
#include "bcf_filter_arg.h"
#include "bcf_ordered_reader.h"
#include "bcf_ordered_writer.h"
#include "hts_utils.h"
#include "fa_reader.h"

int32_t cmdVcfAddContext(int32_t argc, char** argv) {
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

  // Locate target chr
  FaReader findx(inFa);
  if (! findx.FindFaiVal(chrom) ) {
    error("Cannot find fai index for this chromosome");
  }
  notice("Initialized reader of fasta file. Index info: %ld\t%ld\t%ld\t%ld\n ", findx.val->len, findx.val->offset, findx.val->line_blen, findx.val->line_len);


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
  std::string ctx, cpg, refcpg, majcpg;
  for(; odr->read(iv); ++k) {
    if ( k % verbose == 0 )
      notice("Processing %d markers at %s:%d.", k, bcf_hdr_id2name(odr->hdr, iv->rid), iv->pos+1);

    bcf1_t* nv = bcf_dup(iv);
    bcf_unpack(nv, BCF_UN_ALL);
    if (!withGT)
      bcf_subset(odw.hdr, nv, 0, 0);
    ctx = findx.fa_kmer(iv->pos, kmer);
    cpg = "F";

    if (bcf_is_snp(iv)) { // SNP
      if ((ctx[mpt]=='C' && ctx[mpt+1]=='G')||(ctx[mpt]=='G' && ctx[mpt-1]=='C')) { // CpG
        cpg="T";
      }
    }

    bcf_update_info_string(odw.hdr, nv, "CpG", cpg.c_str());
    bcf_update_info_string(odw.hdr, nv, "CONTEXT", ctx.c_str());
    odw.write(nv);
    bcf_destroy(nv);
  }

  odw.close();
  return 0;
}
