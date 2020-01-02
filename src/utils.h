#ifndef UTILS_H
#define UTILS_H

#include <sys/stat.h>
#include <cstdio>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
//#include <cstdint>
#include <cstring>
#include <cmath>
#include <cfloat>
#include <climits>
#include <vector>
#include <map>
#include <set>
#include <queue>

#include <random>
#include "htslib/kstring.h"
#include "htslib/khash.h"
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/vcf.h"
#include "htslib/vcfutils.h"

/**
 * Splits a line into a vector - PERL style
 */
void split(std::vector<std::string>& vec, const char* delims, std::string& str, uint32_t limit=UINT_MAX, bool clear=true, bool collapse=true);

/**
 * Splits a line into a vector - PERL style
 */
void split(std::vector<std::string>& vec, const char* delims, const char* str, uint32_t limit=UINT_MAX, bool clear=true, bool collapse=true);


/**
 * Casts a string into int32.  Returns true if successful.
 */
bool str2int32(std::string& s, int32_t& i);

/**
 * Casts a string into uint32.  Returns true if successful.
 */
bool str2uint32(std::string& s, uint32_t& i);

/**
 * Appends cuurent working directoy to a path.
 * Returns true if successful.
 */
bool append_cwd(std::string& path);

unsigned int str_hash(const char* s, unsigned int seed = 0);

// New

/**
 * Get the last line of a file .
 */
std::string GetLastLine(const std::string& file);

/**
 * Generate all possible k-mer from an alphabet.
 */
void EnumerateMotif(std::vector<char>& alphabet, int32_t k, std::vector<std::string >& res);
int32_t AllConfig(std::vector<char>& alphabet, int32_t k, std::vector<std::string >& res);

/**
 * Sample without replacement.
 */
void NchooseK(int32_t N, int32_t k, std::set<int32_t> & chosen, std::mt19937& rng, int32_t base = 1);
void NchooseK(int32_t N, int32_t k, std::set<int32_t> & chosen, std::vector<int32_t> & avoid, std::mt19937& rng, int32_t base = 1);

#endif



