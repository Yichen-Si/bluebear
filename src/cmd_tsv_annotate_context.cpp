#include "cramore.h"
#include "tsv_reader.h"
#include "hts_utils.h"
#include "utils.h"
#include "fa_reader.h"
#include "seq_basics.h"

// Goal: annotate tsv file with sequence context
// notes: 1 currently it does not process more than one chromosome
int32_t cmdTsvAddContext(int32_t argc, char** argv) {
    std::string inTsv, inFa, chrom, reg, out, filter;
    int32_t kmer = 9, fold = 0;
    int64_t verbose = 50000;
    int32_t has_tbi = 0;
    int32_t header_row = 0, pos_column = -1, filter_column = -1, filter_soft = 0;
    bool if_filter = 0;

  paramList pl;
    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input", NULL)
    LONG_STRING_PARAM("in-tsv",&inTsv, "Input VCF file to annotate")
    LONG_STRING_PARAM("in-fasta",&inFa, "Input Fasta file")
    LONG_STRING_PARAM("chr",&chrom,"Chromosome to process")
    LONG_STRING_PARAM("region",&reg,"Region to process (require tabix)")
    LONG_INT_PARAM("tbi",&has_tbi,"If input has tabix index")

    LONG_PARAM_GROUP("Additional Options", NULL)
    LONG_INT_PARAM("kmer",&kmer, "Kmer length around each variant")
    LONG_INT_PARAM("fold",&fold, "If reverse complemented sequences are collapsed")
    LONG_INT_PARAM("header-row",&header_row, "Which row contains column names (1-based, 0 for no header)")
    LONG_INT_PARAM("pos-column",&pos_column, "Which column contains variant position (0-based)")
    LONG_STRING_PARAM("filter-match",&filter, "To enforce a pattern matching. --filter-column is required")
    LONG_INT_PARAM("filter-column",&filter_column, "Which column contains string to match with the pattern for filtering (0-based)")
    LONG_INT_PARAM("filter-soft",&filter_soft, "If 'containing pattern' is enough (1) or 'exact match' (0)")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output file name")

  END_LONG_PARAMS();
  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

if ( inTsv.empty() || inFa.empty() || out.empty() ) {
    error("[E:%s:%d %s] --in-tsv, --in-fasta, --out are required parameters",__FILE__,__LINE__,__FUNCTION__);
}
if (!filter.empty()) {
	if_filter = 1;
	if (filter_column < 0) {
		notice("The default column containing the source of filtering is the 6-th (0-based)");
		filter_column = 6;
	}
}
if (reg.empty() && chrom.empty()) {
		error("[E:%s:%d %s] At least one of --chrom & --region is required",__FILE__,__LINE__,__FUNCTION__);
} else {
	if (chrom.empty()) {
		uint32_t i = 0;
		while (i < reg.size()) {
			if (reg[i] == ':') {break;}
			i++;
		}
		chrom = reg.substr(0, i);
	}
}

  // input & output

tsv_reader tr = (inTsv.c_str());
htsFile *wf = hts_open(out.c_str(), "w");
std::string line;
int32_t nVariants = 0;
if (header_row > 0) {
	while (nVariants < header_row && tr.read_line(line) > 0) {
		nVariants += 1;
	}
	hprintf(wf, "%s\t%s\t%s\n", line.c_str(), "CONTEXT", "CpG");
	nVariants = 0;
	if (pos_column < 0 && tr.nfields > 1) {
		pos_column = 0;
		while( pos_column < tr.nfields) {
			std::string cname(tr.str_field_at(pos_column));
			if (cname.compare("POS") == 0) {break;}
			pos_column++;
		}
		if (pos_column >= tr.nfields) {
			error("Which column to read as genome position is unknown");
		}
	}
}
if (!reg.empty()) {
	tr.jump_to(reg.c_str());
}

// Locate target chr
FaReader findx(inFa);
if (! findx.FindFaiVal(chrom) ) {
error("Cannot find fai index for this chromosome");
}
notice("Initialized reader of fasta file. Index info: %ld\t%ld\t%ld\t%ld\n ", findx.val->len, findx.val->offset, findx.val->line_blen, findx.val->line_len);

std::set<char> focal_base_set {'A', 'C', 'N'};
std::map<char, char> bpair = { {'A','T'}, {'T','A'}, {'C','G'}, {'G','C'}, {'N','N'} };

// Process each variant and retrieve kmer context
int32_t pos;
int32_t mid_p = (int32_t) (kmer/2);
std::string ctx;
while(tr.read_line(line) > 0) {
	if (if_filter) {
		std::string field = tr.str_field_at(filter_column);
		if (filter_soft) {
			if (field.find(filter) == std::string::npos) {
				continue;
			}
		} else {
			if (field.compare(filter) != 0) {
				continue;
			}
		}
	}
	pos = tr.int_field_at(pos_column);
    ctx = findx.fa_kmer(pos, kmer);
	bool cpg = if_cpg(ctx, mid_p, 0);
	if ( fold && focal_base_set.find(ctx[mid_p]) == focal_base_set.end() ) {
		ctx = reverse_complement(ctx, bpair);
	}
	hprintf(wf, "%s\t%s\t%d\n", line.c_str(), ctx.c_str(), cpg);
	if ( nVariants % verbose == 0 )
    notice("Processing %d markers at %s:%d.", nVariants, chrom.c_str(), pos);
    nVariants++;
}

	tr.close();
	hts_close(wf);

  return 0;
}
