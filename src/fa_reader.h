#ifndef FA_READER_H
#define FA_READER_H

#include "utils.h"
#include "hts_utils.h"

struct Faival {
  int64_t len, offset, line_blen, line_len;
  Faival(int64_t _l, int64_t _os, int64_t _bl, int64_t _ll) : len(_l), offset(_os), line_blen(_bl), line_len(_ll) {}
};

class FaReader {

    public:
        std::map<std::string, Faival*> lookup;
        faidx_t* fai;
        Faival*  val;
        std::string inFa, inFai;

        FaReader(std::string _inFa, std::string _inFai="") : inFa(_inFa), inFai(_inFai) {
            fai = fai_load(inFa.c_str());
            if (inFai == "") {
            	inFai = inFa + ".fai";
            }
        std::ifstream mfile(inFai);
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

    bool FindFaiVal(std::string &chrom) {
        if (lookup.find(chrom) != lookup.end()) {
    		val = lookup[chrom];
    		return 1;
        } else {
            return 0;
        }
    }

		std::string fa_get_seq(int32_t st, int32_t ed) {
            if (ed > val->len) {return "";}
            if (st < 0) {
                notice("[FaReader::fa_get_seq] Received negative query position, forced to zero");
                st = 0;
            }
            if (ed < st) {
                error("[FaReader::fa_get_seq] Invalid query region");
            }
            int32_t ret = bgzf_useek(fai->bgzf, val->offset + st / val->line_blen * val->line_len + st % val->line_blen, SEEK_SET);
            if ( ret<0 ) {
                error("Error: bgzf_useek failed at %d", st);
            }
            int32_t seq_len = ed-st+1;
            int32_t l = 0;
            char c;
            std::string seq = "";
            while ( (c=bgzf_getc(fai->bgzf))>=0 && l < seq_len ) {
                if (isgraph(c)) {
                    c = toupper(c);
                    seq += c;
                    l++;
                }
            }
            return(seq);
		}

		std::string fa_kmer(int32_t pos, int32_t k) {
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
		}
};

#endif
