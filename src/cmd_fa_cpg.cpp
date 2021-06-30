#include "cramore.h"
#include "hts_utils.h"
#include "utils.h"
#include "fa_reader.h"
#include "seq_basics.h"
#include "fstream"
#include "sstream"



// Output all CpG sites
int32_t cmdFaCpG(int32_t argc, char** argv) {
	std::string inFa, out, chrom;
	int32_t start = 0, end = -1;
	int32_t kmer = 9, kmer0 = 3;
	int32_t verbose = 1000000;
	int32_t step = 50000;

	paramList pl;
	BEGIN_LONG_PARAMS(longParameters)
		LONG_PARAM_GROUP("Input", NULL)
		LONG_STRING_PARAM("in-fasta",&inFa, "Input Fasta file")
		LONG_STRING_PARAM("chr",&chrom,"Chromosome to process. Should exist in the fasta file")

		LONG_PARAM_GROUP("Additional Options", NULL)
		LONG_INT_PARAM("kmer",&kmer, "k-mer")
		LONG_INT_PARAM("chunk-size",&step, "Process by chunk")
		LONG_INT_PARAM("start",&start, "Start position")
		LONG_INT_PARAM("end",&end, "End position")

		LONG_PARAM_GROUP("Output Options", NULL)
		LONG_STRING_PARAM("out", &out, "Output file")
		LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (bp)")

	END_LONG_PARAMS();
	pl.Add(new longParams("Available Options", longParameters));
	pl.Read(argc, argv);
	pl.Status();

	std::string inFai = inFa + ".fai";
	// sanity check of input arguments
	if ( inFa.empty() || out.empty() ) {
	error("[E:%s:%d %s] --in-fasta, --out are required parameters",__FILE__,__LINE__,__FUNCTION__);
	}
	if (!file_exists(inFa) || !file_exists(inFai)) {
	error("[E:%s:%d %s] input --in-fasta or its .fai index file does not exist",__FILE__,__LINE__,__FUNCTION__);
	}

	FaReader findx(inFa);
	if (!findx.FindFaiVal(chrom)) {
		error("Cannot find fai index for this chromosome");
	}
	notice("Initialized reader of fasta file. Index info: %ld\t%ld\t%ld\t%ld\n ", findx.val->len, findx.val->offset, findx.val->line_blen, findx.val->line_len);

	htsFile *wf = hts_open(out.c_str(), "w");
	notice("Output saves to \n%s", out.c_str());

	std::string word0, word;
	std::set<char> focal_base_set {'A', 'C', 'N'};
	std::map<char, char> bpair = { {'A','T'}, {'T','A'}, {'C','G'}, {'G','C'}, {'N','N'} };

	int32_t verb = 0, pos = 0, nvb = 0;
	int32_t st = start;
	int32_t ed = st + step - 1;
	int32_t kmer_pad_left = (int32_t) (kmer/2);
	int32_t kmer_pad_right= kmer - kmer_pad_left - 1;
	if (st < kmer_pad_left) {
		st = kmer_pad_left;
	}

	while (1) {
		if ( st > nvb * verbose ) {
			notice("Processed %d CpG sites, covering %.2f Mb", verb, st*1e-6);
			nvb++;
		}
		// Read sequence by window
		std::string seq = findx.fa_get_seq(st-kmer_pad_left, ed+kmer_pad_right);
		if (seq.size() == 0) {
			break;
		}
		std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
		// Count kmer occurence
		pos = st-1;
		for(int32_t i=kmer_pad_left; i < (int32_t)seq.size()-kmer_pad_right; ++i) {
			pos++;
			word0=seq.substr(i-1,kmer0);
			if (!if_cpg(word0,1,0)) {continue;}
				word=seq.substr(i-kmer_pad_left,kmer);
			if ( focal_base_set.find(word[kmer_pad_left]) == focal_base_set.end() ) {
				word = reverse_complement(word, bpair);
				word0= word.substr( kmer_pad_left-1, kmer0 );
			}
			hprintf(wf, "%s\t%d\t%s\t%s\n", chrom.c_str(), pos, word0.c_str(), word.c_str() );
			verb += 1;
		}
		st += step;
		ed = st + step - 1;
		if (end > 0 && st > end) {break;}
	}

  hts_close(wf);
  return 0;
}










// Output all CpG sites into bed files with strand info
// chrom start end strand kmer
int32_t cmdFaCpG_bed(int32_t argc, char** argv) {
	std::string inFa, out, chrom;
	int32_t start = 0, end = -1;
	int32_t kmer = 9, kmer0 = 3;
	int32_t verbose = 1000000;
	int32_t step = 50000;

	paramList pl;
	BEGIN_LONG_PARAMS(longParameters)
		LONG_PARAM_GROUP("Input", NULL)
		LONG_STRING_PARAM("in-fasta",&inFa, "Input Fasta file")
		LONG_STRING_PARAM("chr",&chrom,"Chromosome to process. Should exist in the fasta file")

		LONG_PARAM_GROUP("Additional Options", NULL)
		LONG_INT_PARAM("kmer",&kmer, "k-mer")
		LONG_INT_PARAM("chunk-size",&step, "Process by chunk")
		LONG_INT_PARAM("start",&start, "Start position")
		LONG_INT_PARAM("end",&end, "End position")

		LONG_PARAM_GROUP("Output Options", NULL)
		LONG_STRING_PARAM("out", &out, "Output file")
		LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (bp)")

	END_LONG_PARAMS();
	pl.Add(new longParams("Available Options", longParameters));
	pl.Read(argc, argv);
	pl.Status();

	std::string inFai = inFa + ".fai";
	// sanity check of input arguments
	if ( inFa.empty() || out.empty() ) {
	error("[E:%s:%d %s] --in-fasta, --out are required parameters",__FILE__,__LINE__,__FUNCTION__);
	}
	if (!file_exists(inFa) || !file_exists(inFai)) {
	error("[E:%s:%d %s] input --in-fasta or its .fai index file does not exist",__FILE__,__LINE__,__FUNCTION__);
	}

	FaReader findx(inFa);
	if (!findx.FindFaiVal(chrom)) {
		error("Cannot find fai index for this chromosome");
	}
	notice("Initialized reader of fasta file. Index info: %ld\t%ld\t%ld\t%ld\n ", findx.val->len, findx.val->offset, findx.val->line_blen, findx.val->line_len);

	htsFile *wf = hts_open(out.c_str(), "w");
	notice("Output saves to \n%s", out.c_str());

	std::string word0, word;
	std::set<char> focal_base_set {'A', 'C', 'N'};
	std::map<char, char> bpair = { {'A','T'}, {'T','A'}, {'C','G'}, {'G','C'}, {'N','N'} };
	std::map<int32_t, char> strand_label = { { 1, '+' }, { -1, '-' } };

	int32_t verb = 0, pos = 0, nvb = 0;
	int32_t st = start;
	int32_t ed = st + step - 1;
	int32_t kmer_pad_left = (int32_t) (kmer/2);
	int32_t kmer_pad_right= kmer - kmer_pad_left - 1;
	if (st < kmer_pad_left) {
		st = kmer_pad_left;
	}

	while (1) {
		if ( st > nvb * verbose ) {
			notice("Processed %d CpG sites, covering %.2f Mb", verb, st*1e-6);
			nvb++;
		}
		// Read sequence by window
		std::string seq = findx.fa_get_seq(st-kmer_pad_left, ed+kmer_pad_right);
		if (seq.size() == 0) {
			break;
		}
		std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
		// Count kmer occurence
		pos = st-1;
		for(int32_t i=kmer_pad_left; i < (int32_t)seq.size()-kmer_pad_right; ++i) {
			pos++;
			word0=seq.substr(i-1,kmer0);
			int32_t strand = which_cpg(word0,1);
			if (strand == 0) {continue;}
				word=seq.substr(i-kmer_pad_left,kmer);
			hprintf(wf, "%s\t%d\t%d\t%c\t%s\t%s\n", chrom.c_str(), pos, pos+1, strand_label[strand], word0.c_str(), word.c_str() );
			verb += 1;
		}
		st += step;
		ed = st + step - 1;
		if (end > 0 && st > end) {break;}
	}
  hts_close(wf);
  return 0;
}
