#ifndef SEQ_BASICS_H
#define SEQ_BASICS_H

#include <string>
#include <vector>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <map>
#include <cmath>
#include <cfloat>

inline bool if_cpg (const std::string& seq, int32_t indx, bool one_side=1) {
	if (one_side) {
		return seq.compare(indx,2,"CG") == 0;
	} else {
		return (seq.compare(indx,2,"CG") == 0) || (seq.compare(indx-1,2,"CG") == 0);
	}
}

inline int32_t which_cpg (const std::string& seq, int32_t indx) {
	if (seq.compare(indx,2,"CG") == 0) {
		return 1;
	} else if (seq.compare(indx-1,2,"CG") == 0) {
		return -1;
	} else {
		return 0;
	}
}

inline std::string reverse_complement (std::string& seq, std::map<char, char>& bpair, char default_missing = 'N') {
	std::stringstream reversed;
	int32_t seq_len = seq.size();
	for (int32_t i = 1; i <= seq_len; i++) {
		auto ptr = bpair.find(seq[seq_len-i]);
		if (ptr == bpair.end()) {
			reversed << default_missing;
		} else {
			reversed << bpair[seq[seq_len-i]];
		}
	}
	return reversed.str();
}

#endif
