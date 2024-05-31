#include "utils.h"

/**
 * Splits a line into a vector - PERL style
 */
void split(std::vector<std::string>& vec, const char *delims, std::string& str, uint32_t limit, bool clear, bool collapse)
{
    std::map<char, int32_t> delim_set;

    for (uint32_t i=0; i<strlen(delims); ++i)
    {
        delim_set[delims[i]] = 1;
    }

    if (clear)
    {
        vec.clear();
    }
    const char* tempStr = str.c_str();
    int32_t i=0, lastIndex = str.size()-1;
    std::stringstream token;

    if (lastIndex<0) return;

    uint32_t noTokens = 0;
    bool isDelim = false;
    while (i<=lastIndex)
    {
        isDelim = (delim_set.find(tempStr[i])!=delim_set.end());

        if (!isDelim || noTokens>=limit-1)
        {
            token << tempStr[i];
        }

        if ((isDelim && noTokens<limit-1) || i==lastIndex)
        {
            if (collapse && token.str()!="")
            {
                vec.push_back(token.str());
                ++noTokens;
                token.str("");
            }
        }

        ++i;
    }
};

/**
 * Splits a line into a vector - PERL style
 */
void split(std::vector<std::string>& vec, const char *delims, const char* str, uint32_t limit, bool clear, bool collapse)
{
    std::map<char, int32_t> delim_set;

    for (uint32_t i=0; i<strlen(delims); ++i)
    {
        delim_set[delims[i]] = 1;
    }

    if (clear)
    {
        vec.clear();
    }
    const char* tempStr = str;
    int32_t i=0, lastIndex = strlen(str)-1;
    std::stringstream token;

    if (lastIndex<0) return;

    uint32_t noTokens = 0;
    bool isDelim = false;
    while (i<=lastIndex)
    {
        isDelim = (delim_set.find(tempStr[i])!=delim_set.end());

        if (!isDelim || noTokens>=limit-1)
        {
            token << tempStr[i];
        }

        if ((isDelim && noTokens<limit-1) || i==lastIndex)
        {
            if (collapse && token.str()!="")
            {
                vec.push_back(token.str());
                ++noTokens;
                token.str("");
            }
        }

        ++i;
    }
};

/**
 * Casts a string into int32.  Returns true if successful.
 */
bool str2int32(std::string& s, int32_t& i)
{
    const char* start = s.c_str();
    char *end = 0;
    i = std::strtol(s.c_str(), &end, 10);
    return (end!=start);
};

/**
 * Casts a string into uint32.  Returns true if successful.
 */
bool str2uint32(std::string& s, uint32_t& i)
{
    const char* start = s.c_str();
    char *end = 0;
    i = std::strtoul(s.c_str(), &end, 10);
    return (end!=start);
};

/**
 * Appends cuurent working directoy to a path.
 * Returns true if successful.
 */
bool append_cwd(std::string& path)
{
    if (path.size()>0 && path.c_str()[0]!='/')
    {
        char cwd[1024];
        if (getcwd(cwd, sizeof(cwd))!=NULL)
        {
            std::string cwd_path(cwd);
            path = cwd_path + "/" + path;

            return true;
        }
    }

    return false;
};


unsigned int str_hash(const char* s, unsigned int seed)
{
  unsigned int hash = seed;
  while (*s) {
    hash = hash * 101  +  *s++;
  }
  return hash;
};


/**
 * Get the last line of a file .
 */
std::string GetLastLine(const std::string& f) {
    std::ifstream file(f);
    file.seekg(-1, std::ios_base::end);
    char c;
    file.get(c);
    while (c == '\n') {
      file.seekg(-2,std::ios_base::cur);
      file.get(c);
    }
    while (c != '\n') {
      file.seekg(-2,std::ios_base::cur);
      file.get(c);
    }
    std::string line;
    std::getline(file, line);
    return line;
};


/**
 * Generate all possible k-mer from an alphabet.
 */
void EnumerateMotif(std::vector<char>& alphabet, int32_t k, std::vector<std::string >& res) {
  uint32_t n = alphabet.size();
  if (k == 1) {
    for (uint32_t i = 0; i < n; ++i) {
      std::string tmp(1,alphabet[i]);
      res.push_back(tmp);
    }
    return;
  }
  std::vector<std::string > sres;
  EnumerateMotif(alphabet, k-1, sres);
  for (uint32_t i = 0; i < n; ++i) {
    std::vector<std::string > nres = sres;
    for (uint32_t j = 0; j < nres.size(); j++) {
      nres[j] += alphabet[i];
    }
    res.insert(res.end(), nres.begin(), nres.end());
  }
  return;
};

int32_t AllConfig(std::vector<char>& alphabet, int32_t k, std::vector<std::string >& res) {
  EnumerateMotif(alphabet,k,res);
  return res.size();
};

/**
 * Sample without replacement.
 */
void NchooseK(int32_t N, int32_t k, std::set<int32_t> & chosen, std::mt19937& rng, int32_t base) {
    if (N < 3 * k) {
      if (N == k) {
        for (int32_t r = base; r < N; ++r) {
          chosen.emplace_hint(chosen.end(), r);
        }
      }
      std::vector<int32_t> v(N);
      std::iota(v.begin(), v.end(), base);
      std::shuffle(v.begin(), v.end(), rng);
      for (int32_t r = 0; r < k; ++r) {
        chosen.insert(v[r]);
      }
    } else {
      for (int32_t r = N - k - 1 + base ; r < N - 1 + base; ++r) {
        int32_t v = std::uniform_int_distribution<>(base, r)(rng);
        if (!chosen.insert(v).second) {
          chosen.insert(r+1);
        }
      }
    }
};

void NchooseK(int32_t N, int32_t k, std::set<int32_t> & chosen, std::vector<int32_t> & avoid, std::mt19937& rng, int32_t base) {
  int32_t n = N, k0 = 0, i = 0;
  std::vector<int32_t> v(N);
  std::iota(v.begin(), v.end(), base);
  for (auto a : avoid) {
    if (a < base || a >= N + base) {
      continue;
    }
    v[a] = -1;
    n--;
  }
  if (n < k) {
    error("NchooseK: not enough elements to choose from");
  }
  std::shuffle(v.begin(), v.end(), rng);
  while (k0 < k) {
    if (v[i] < 0) {
      continue;
    }
    chosen.emplace_hint(chosen.end(), v[i]);
    k0++;
  }
};

/**
 * Generate random & even partitions of indices.
*/
void generatePartitions(int32_t totalIndices, int32_t nPartitions, std::vector<std::vector<size_t>>& partitions, std::mt19937& rng) {
  std::vector<size_t> indices(totalIndices);
  for (size_t i = 0; i < totalIndices; ++i) {
      indices[i] = i;
  }
  std::shuffle(indices.begin(), indices.end(), rng);
  // Divide shuffled indices into nearly even partitions
  partitions.clear();
  partitions.resize(nPartitions);
  size_t perPartition = totalIndices /nPartitions;
  size_t extra = totalIndices % nPartitions;
  size_t start = 0, end = 0;
  for (int i = 0; i < nPartitions; ++i) {
      end += perPartition + (i < extra ? 1 : 0);
      partitions[i] = std::vector<size_t>(indices.begin() + start, indices.begin() + end);
      start = end;
  }
}

template<typename T>
std::vector<int32_t> argsort(std::vector<T> &array) {
    std::vector<int32_t> indices(array.size());
    std::iota(indices.begin(), indices.end(), 0);
    sort(indices.begin(), indices.end(),
      [&array](int32_t left, int32_t right) -> bool {
          return array[left] <= array[right]; } );
    return indices;
};

template<typename T, typename U>
std::vector<std::pair<T, U> > sort_map_val(std::map<T,  U>& M)
{
    std::vector<std::pair<T, U> > A;
    for (auto& it : M) {
        A.push_back(it);
    }
    sort(A.begin(), A.end(), cmp_pair2<T, U>);
    return A;
};
