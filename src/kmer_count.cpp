#include "cramore.h"
#include "hts_utils.h"
#include "utils.h"
#include "fa_reader.h"
#include "seq_basics.h"
#include "fstream"
#include "sstream"

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

  // Read .fa
	FaReader findx(inFa);
  if (!findx.FindFaiVal(chrom)) {
    error("Cannot find fai index for this chromosome");
  }
	notice("Initialized .fa reader");

  std::string seq = "";
  std::set<char> alphabet {'A','T','C','G'};
  std::map<std::string, int32_t> kmerct;
  char c;
  int32_t l = 0;

  l = 0;
  while ( (c=bgzf_getc(findx.fai->bgzf))>=0 && l < kmer ) {
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
  while ((c=bgzf_getc(findx.fai->bgzf))>=0) {
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
    hprintf(wf, "%s\t%d\t%s\t%d\n", chrom.c_str(), v.first.c_str(), v.second);
  }

  hts_close(wf);
  return 0;
}









// Goal: count the occurence of kmers by window
// Output one window per line or one window x kmer per line

int32_t KmerCount_window(int32_t argc, char** argv) {
  std::string inFa, out, outf;
  std::string chrom="chr20";
  int32_t kmer = 3;
  int32_t long_format = 0, cpg_only = 1;
  int32_t start = 0, end = -1;
  int32_t verbose = 1000;
  int32_t w_size = 1000, ovlp = 500;

  paramList pl;
  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input", NULL)
    LONG_STRING_PARAM("in-fasta",&inFa, "Input Fasta file")
    LONG_STRING_PARAM("chr",&chrom,"Chromosome to process. Should exist in the fasta file")

    LONG_PARAM_GROUP("Additional Options", NULL)
    LONG_INT_PARAM("kmer",&kmer, "k-mer")
    LONG_INT_PARAM("window-size",&w_size, "Window size (bp) to count kmer occurence & output")
    LONG_INT_PARAM("overlap",&ovlp, "Overlap between adjacent windows")
    LONG_INT_PARAM("cpg-only",&cpg_only, "Whether to only include kmers with CpG in the middle")
    LONG_INT_PARAM("start",&start, "Start position")
    LONG_INT_PARAM("end",&end, "End position")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output file prefix")
    LONG_INT_PARAM("long",&long_format, "Whether to output in the long (1) or wide (0) format")
		LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n window)")

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
	if (end < start + w_size) {
		error("End point --end is required");
	}

  int32_t kmer_pad_left = (int32_t) (kmer/2);
  int32_t kmer_pad_right= kmer - kmer_pad_left - 1;
  int32_t step = w_size - ovlp;
  if (cpg_only && (kmer % 2 == 0) )
    notice("Kmer size is even (%d), CpG C site is define as the %d-th position", kmer, (int) (kmer/2+1));
  if ( w_size % (w_size-ovlp) != 0 )
    notice("Window size & overlapping length privided does not support a partition as a subset of the output windows" );

  FaReader findx(inFa);
  if (!findx.FindFaiVal(chrom)) {
    error("Cannot find fai index for this chromosome");
  }
  notice("Initialized reader of fasta file. Index info: %ld\t%ld\t%ld\t%ld\n ", findx.val->len, findx.val->offset, findx.val->line_blen, findx.val->line_len);

  outf = out + "." + chrom + ".kmer_" + std::to_string(kmer) + "_wsize_"+std::to_string((int32_t) (w_size*1e-3))+"kb_ovlp_"+std::to_string((ovlp))+"bp.count";
	if (cpg_only) {
		outf += ".cpg_only.bed";
	} else {
		outf += ".bed";
	}
  htsFile *wf = hts_open(outf.c_str(), "w");
	notice("Output saves to \n%s", outf.c_str());

	std::map<std::string, int32_t> kmerct;
	std::map<std::string, int32_t> kmerct_n;
  std::string word;

	std::vector<char> alphabet {'A','T','C','G'};
	std::set<char> alphabet_set {'A','T','C','G'};
	std::vector<char> alphabet_n {'A','T','C','G','N'};
	std::vector<char> focal_base {'A', 'C', 'N'};
	std::set<char> focal_base_set {'A', 'C', 'N'};
	std::map<char, char> bpair = { {'A','T'}, {'T','A'}, {'C','G'}, {'G','C'}, {'N','N'} };
	std::vector<std::string > motifs;

	// Define an output kmer sequence if CpG && k=3
	// Otherwise, output non-zero counts for each window
	std::vector<std::string> output_required {"ACG","CCG","GCG","TCG"};

	if (cpg_only) {
		// Generate all possible k-mers (folded) with C in CG being the middle base
		int32_t n_comb = AllConfig(alphabet_n,kmer-2,motifs);
		notice("There are %d possible CpG kmer configurations", n_comb);
		for (int32_t i = 0; i < n_comb; ++i) {
			word = motifs[i].substr(0,kmer_pad_left) + "CG" + motifs[i].substr(kmer_pad_left);
			if (word.find('N') == std::string::npos)
				kmerct[word] = 0;
			else
				kmerct_n[word] = 0;
		}
	} else {
		// Generate all possible k-mers (folded)
		// Set middle (focal) bp to be A or C (or N), flip to the negative strand otherwise
		int32_t n_comb = AllConfig(alphabet_n,kmer-1,motifs);
		notice("There are %d possible folded kmer configurations", n_comb*focal_base.size());
		for (int32_t i = 0; i < n_comb; ++i) {
			for (auto x : focal_base) {
				word = motifs[i].substr(0,kmer_pad_left);
				word += x;
				word += motifs[i].substr(kmer_pad_left);
				if (word.find('N') == std::string::npos)
					kmerct[word] = 0;
				else
					kmerct_n[word] = 0;
			}
		}
	}

for (auto& kv : kmerct)	{
	std::cout << kv.first << '\t';
}
std::cout << '\n';
for (auto& kv : kmerct_n)	{
	std::cout << kv.first << '\t';
}
std::cout << '\n';

	int32_t verb = 0;
	int32_t st = start;
	int32_t ed = st + w_size - 1;

	while ( ed < end ) {

		if ( verb % verbose == 0 )
			notice("Processing %d-th window in %s at %d kb (%.2f Mb)", verb, chrom.c_str(), st/1000, st*1e-6);

		// Read sequence by window
		std::string seq = findx.fa_get_seq(st-kmer_pad_left, ed+kmer_pad_right);
		if (seq.size() == 0) {
			st += step;
			ed = st + w_size - 1;
			continue;
		}
		std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
		// Count kmer occurence
		for(int32_t i=kmer_pad_left; i < (int32_t)seq.size()-kmer_pad_right; ++i) {
			word=seq.substr(i-kmer_pad_left,kmer);
			// bool if_cpg (const std::string& seq, int32_t indx, bool one_side=1)
			if (cpg_only && !if_cpg(word,kmer_pad_left,0)) {
				continue;
			}
			if ( focal_base_set.find(word[kmer_pad_left]) == focal_base_set.end() ) {
				word = reverse_complement(word, bpair);
			}
			auto ptr = kmerct.find(word);
			if (ptr == kmerct.end()) {
				auto ptr_n = kmerct_n.find(word);
				if (ptr_n == kmerct_n.end()) { // There are other letters than ACTGN
					std::stringstream word_replace;
					for (uint32_t i = 0; i < word.size(); ++i) {
						if (alphabet_set.find(word[i]) == alphabet_set.end()) {
							word_replace << 'N';
						} else {
							word_replace << word[i];
						}
					}
					word = word_replace.str();
					auto ptr_rd = kmerct_n.find(word);
					if (ptr_rd == kmerct_n.end()) { // Should never happen
						error("sth is wrong: incountered kmer %s not in pre-enumerated legal set", word.c_str());
					} else {
						ptr_rd->second++;
					}
				} else {
					ptr_n->second++;
				}
			} else {
				ptr->second++;
			}
		}
		// Output counts for this window
		if (long_format==1) {
			for (auto const& kv : kmerct) {
				if (kv.second > 0)
					hprintf(wf, "%s\t%d\t%d\t%s\t%d\n", chrom.c_str(), st, st + w_size, kv.first.c_str(), kv.second);
		  }
			for (auto const& kv : kmerct_n) {
				if (kv.second > 0)
					hprintf(wf, "%s\t%d\t%d\t%s\t%d\n", chrom.c_str(), st, st + w_size, kv.first.c_str(), kv.second);
		  }
		} else {
			std::stringstream ss;
			if (cpg_only && kmer == 3) {
				for (auto k : output_required) {
					ss << k << ":" << kmerct[k] << ";";
				}
			} else {
				for (auto const& kv : kmerct) {
					if (kv.second > 0)
						ss << kv.first << ":" << kv.second << ";";
			  }
			}
			for (auto const& kv : kmerct_n) {
				if (kv.second > 0)
					ss << kv.first << ":" << kv.second << ";";
			}
			hprintf(wf, "%s\t%d\t%d\t%s\n", chrom.c_str(), st, st + w_size, ss.str().c_str());
		}
		// Clear up for next window
		for (auto& kv : kmerct) {
			kv.second = 0;
		}
		for (auto& kv : kmerct_n) {
			kv.second = 0;
		}
		st += step;
		ed = st + w_size - 1;
		verb += 1;
	}

  hts_close(wf);
  return 0;
}
