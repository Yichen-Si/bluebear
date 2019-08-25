#include "cramore.h"
#include "utils.h"

struct interval_t {
  int32_t st, ed;
  int32_t lbp, min_ac;
  double  lcm;
  // std::vector<int32_t> ac;
  // std::vector<string> idlist;
};

class ptowner {
public:
  std::string id;
  int32_t tot = 0, inside = 0, outside = 0;
  std::vector<int32_t> pos;
  std::vector<interval_t*> regions;

  ptowner(std::string _id, std::vector<int32_t> _pos) : id(_id), pos(_pos) {}

  void count_pt() {
    tot = pos.size();
    int32_t lower = 0;
    int32_t k = 0;
    while (k < (int32_t) pos.size()) {
      bool flag = 0;
      if (regions[lower]->st > pos[k]) {
        k++;
        continue;
      }
      for (int32_t i = lower; i < (int32_t) regions.size(); ++i) {
        if (regions[i]->ed < pos[k]) {
          lower = i+1;
          break;
        }
      }
      if (lower >= (int32_t) regions.size()) {
        break;
      }
      for (int32_t i = lower; i < (int32_t) regions.size(); ++i) {
        if (regions[i]->st <= pos[k] && regions[i]->ed >= pos[k]) {
          flag = 1;
        }
      }
      if (flag) {
        inside++;
      }
      k++;
    }
    outside = tot - inside;
  }

};

int32_t mapPt2Interval(int32_t argc, char** argv) {

  std::string rfile, pfile, out, line;
  std::string tmp="4,5";
  int32_t st_col = 0, ed_col = 1, id_col = 0, pt_col = 4;
  std::vector<int32_t> id_col_list;
  std::vector<std::string> words;
  std::map<std::string, ptowner*> owners;

  paramList pl;
  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Sites", NULL)
    LONG_STRING_PARAM("rfile",&rfile, "Input interval file")
    LONG_STRING_PARAM("pfile",&pfile, "Input position file")
    LONG_INT_PARAM("st-col",&st_col, "Col. of start point in interval file, 0-based.")
    LONG_INT_PARAM("ed-col",&ed_col, "Col. of end point in interval file, 0-based.")
    LONG_STRING_PARAM("id-list-col",&tmp, "Col. of involved individuals in interval file, 0-based.")
    LONG_INT_PARAM("id-col",&id_col, "Col. of individual id in position file, 0-based.")
    LONG_INT_PARAM("pt-col",&pt_col, "Col. of points in position file, 0-based.")

    // LONG_PARAM_GROUP("Additional Options", NULL)

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output file")

  END_LONG_PARAMS();
  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  split(words, ",", tmp);
  for (auto & v : words)
    id_col_list.push_back(stoi(v));

  std::ifstream infile;

  infile.open(pfile, std::ifstream::in);
  while (std::getline(infile, line)) {
    split(words, "\t", line);
    std::vector<std::string> posline;
    std::vector<int32_t> posvec;
    split(posline, ",", words[pt_col]);
    for (auto & v: posline) {
      int32_t pt;
      if (!str2int32(v, pt)) {
        error("Incorrect point lists.");
      }
      posvec.push_back(pt);
    }
    owners[words[id_col]] = new ptowner(words[id_col],posvec);
  }
  infile.close();

notice("Finish parsing point file. Found %d individuals", owners.size());

  infile.open(rfile, std::ifstream::in);
  while (std::getline(infile, line)) {
    split(words, "\t", line);
    interval_t * reg = new interval_t;
    if (!str2int32(words[st_col], reg->st) || !str2int32(words[ed_col], reg->ed)) {
      error("Incorrect interval input.");
    }
    for (auto & it : id_col_list) {
      auto mit = owners.find(words[it]);
      if (mit != owners.end())
        (owners[words[it]]->regions).push_back(reg);
    }
  }
  infile.close();
notice("Finish parsing interval file.");

  FILE * wf;
  wf = fopen(out.c_str(), "w");
  for (auto const& ptr : owners) {
    ptr.second->count_pt();
    fprintf(wf,"%s\t%d\t%d\t%d\n", (ptr.first).c_str(), (ptr.second)->tot, (ptr.second)->inside, (ptr.second)->outside);
  }
  fclose(wf);

  return 0;
}














