#include "cramore.h"
#include "utils.h"
#include "bcf_filter_arg.h"
#include "bcf_ordered_reader.h"
#include "compact_matrix.h"
#include <iomanip>
#include "bp2cm.h"
#include "ibs0.h"
#include <numeric>

struct ibs0interval {
public:
  double leftcm, rightcm, totcm;
  ibs0interval(double _l, double _r, double _t) :
           leftcm(_l), rightcm(_r), totcm(_t) {}
};


class Doubleton {
public:
  int32_t pos;
  std::string oldrec;
  double cm0;
  std::vector<int32_t> carriers;   // IDs of the two carriers
  std::vector<ibs0interval> ibs0pairs; // Non-carriers
  std::vector<double> refcm;       // IBS0 length between control pairs
  double mean=0.0, median=0.0, stdev=0.0;
  int32_t rank = 0;

Doubleton(int32_t _p, std::string _r, double _l) : pos(_p), oldrec(_r), cm0(_l) {}

void Add(double left, double right, double tot) {
  ibs0interval rec(left, right, tot);
  ibs0pairs.push_back(rec);
  refcm.push_back(tot);
}

void AddID(int32_t _id) {
  carriers.push_back(_id);
}

void Organize() {
  // Sort IBS0 length, get summary statistics

  std::sort(refcm.begin(), refcm.end());
  int32_t n = ibs0pairs.size();
  while(rank < n && cm0 > refcm[rank]) {
    rank ++;
  }
  median = refcm[n/2-1];
  if (n % 2 == 0) {
    median = (median + refcm[n/2]) / 2.0;
  }
  mean = std::accumulate(refcm.begin(), refcm.end(), 0.0) / n;
  double sq_sum = 0.0;
  for (auto & v : refcm) {
    sq_sum += pow(v-mean, 2);
  }
  stdev = std::sqrt(sq_sum / n);

}

};


// goal -- For each doubleton, add non-carrier pairs as control
int32_t IBS0AddNoncarrierControl(int32_t argc, char** argv) {
  std::string inVcf, inMap, inRec, chrom, reg;
  std::string out, outf;
  int32_t cst = -1, ced = -1;
  int32_t min_hom_gts = 1;
  int32_t verbose = 5000;
  int32_t min_variant = 1;
  int32_t n_samples = 0;
  int32_t start, end;
  int32_t ck_len = 1000000;
  double  cm_limit = 0.5;

  int32_t AC = 2;
  int32_t npairs = 100;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Files", NULL)
    LONG_STRING_PARAM("in-rec",&inRec, "Input variant record")
    LONG_STRING_PARAM("in-vcf",&inVcf, "Input VCF/BCF file")
    LONG_STRING_PARAM("region",&reg,"Genomic region to focus on")
    LONG_STRING_PARAM("map",&inMap, "Map file for genetic distance")

    LONG_PARAM_GROUP("Additional Options", NULL)
    LONG_INT_PARAM("min-hom",&min_hom_gts, "Minimum number of homozygous genotypes to be counted for IBS0")
    LONG_INT_PARAM("min-variant",&min_variant, "Minimum number of variants to present to have output file")
    LONG_DOUBLE_PARAM("mix-ibs0",&cm_limit, "Maximal ibs0 in cM to search for.")

    LONG_INT_PARAM("centromere-st",&cst, "Start position of centromere")
    LONG_INT_PARAM("centromere-ed",&ced, "End position of centromere")

    // For testing
    LONG_INT_PARAM("num-samples",&n_samples, "Number of samples to test")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output file prefix")
    LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")

  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  // sanity check of input arguments
  if ( inRec.empty() || inVcf.empty() || reg.empty() || inMap.empty() || out.empty() ) {
    error("[E:%s:%d %s] --in-rec, --in-vcf, --region, --inMap, --out are required parameters",__FILE__,__LINE__,__FUNCTION__);
  }

  // Genetic map
  bp2cmMap pgmap(inMap, " ", cst, ced);

  // Output
  std::ofstream wf;
  outf = out + "_with_control_" + reg + ".txt";
  wf.open(outf, std::ios::trunc);

  // Region to process
  std::vector<std::string> v;
  split(v, ":-", reg);
  chrom = v[0];
  if (!str2int32(v[1], start) || !str2int32(v[2], end)) {
    error("Invalid region.");
  }

  // Setup record of results
  std::vector<Doubleton*> Doubts;

  // IBS0 lookup
  IBS0lookup ibs0finder(inVcf, reg, pgmap, cm_limit, ck_len, 1);
  notice("Set up IBS0 finder");

  // Read info file, find corresponding region
  // *** Input is assumed to be sorted by position
  std::string line, word;
  int32_t ac, pos;
  std::ifstream rf(inRec);
  while(std::getline(rf, line)) {
    std::istringstream ss(line);
    ss >> chrom >> ac >> pos;
    if (pos >= start) {break;} // Find starting point
  }
  while(std::getline(rf, line) && pos < end) {
    std::istringstream ss(line);
    ss >> chrom >> ac >> pos;
    while(!ss.eof()) {
      ss >> word;
    }
    double l0 = std::stod(word);
    line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
    Doubleton* rec = new Doubleton(pos, line, l0);
    Doubts.push_back(rec);
  }
  notice("Read input records");

  // bcf reader
  std::vector<GenomeInterval> intervals;
  parse_intervals(intervals, "", reg);
  BCFOrderedReader odr(inVcf, intervals);
  bcf1_t* iv = bcf_init();

  int32_t nsamples = bcf_hdr_nsamples(odr.hdr);

  notice("Identifying %d samples from Vcf file", nsamples);
  if ( nsamples == 0 )
    error("FATAL ERROR: The VCF does not have any samples with genotypes");

  // For testing
  if (n_samples > 0 && nsamples > n_samples)
    nsamples = n_samples;

  // First pass: add carrier ids to record to avoid
  int32_t* p_gt = NULL;
  int32_t  n_gt = 0;
  int32_t *info_ac = NULL;
  int32_t n_ac = 0;

  int32_t ptr = 0, matchedct = 0;
  for(int32_t k=0; odr.read(iv); ++k) {

    if (bcf_get_info_int32(odr.hdr, iv, "AC", &info_ac, &n_ac) < 0) {continue;}
    if (info_ac[0] != AC) {continue;}
    if (iv->pos+1 < Doubts[ptr]->pos) {continue;}
    if (iv->pos+1 > Doubts[ptr]->pos) {ptr ++; continue;}
    matchedct ++;

    // Extract genotype and identify carriers
    if ( bcf_get_genotypes(odr.hdr, iv, &p_gt, &n_gt) < 0 ) {
      error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);
    }
    for (int32_t it = 0; it < nsamples; ++it) {
      if (bcf_gt_allele(p_gt[2*it]) > 0 || bcf_gt_allele(p_gt[2*it+1]) > 0) {
        Doubts[ptr]->AddID(it);
      }
    }
    if ( (int32_t) (Doubts[ptr]->carriers).size() < AC) { // Never
      notice("Error at pos %s:%d: wrong number of carriers", chrom, iv->pos+1);
      Doubts[ptr]->carriers[1] = Doubts[ptr]->carriers[0];
    }
  }
  notice("Found pairs of carriers for %d among %d variants", matchedct, Doubts.size());

  // Process doubletons one by one
  int32_t ndone = 0;
  std::vector<int32_t> candy(nsamples);
  std::iota(candy.begin(), candy.end(), 0);
  std::set<int32_t> chosen;
  std::random_device rd;
  std::mt19937 rng(rd());
  std::uniform_int_distribution<int> uni(0,npairs-1);
  for (auto & rec : Doubts) {
    pos = rec->pos;
    if ( (int32_t) (rec->carriers).size() < AC) {continue;}
    if (ndone % verbose == 0) {
      notice("Processing %d variants at %s:%d.", ndone, chrom.c_str(), pos);
    }
// if (ndone > 100) { // Test
//   break;
// }
    // Randomly sample npairs of non-carriers
    NchooseK(nsamples, npairs*2, chosen, rng, 0);
    while (chosen.find(rec->carriers[0]) != chosen.end() || chosen.find(rec->carriers[1]) != chosen.end()) { // Inefficient way
      chosen.clear();
      NchooseK(nsamples, npairs*2, chosen, rng, 0);
    }
    // Find no-IBS0 length for each pair
    std::vector<int32_t> chosen_vec(chosen.begin(), chosen.end());
    for (int32_t it = 0; it < npairs; ++it) {
      int32_t r = ibs0finder.FindIBS0(chosen_vec[it*2],chosen_vec[it*2+1],pos,0);
      int32_t l = ibs0finder.FindIBS0(chosen_vec[it*2],chosen_vec[it*2+1],pos,1);
      if (l < 0) {l = ibs0finder.start_que[0];}
      if (r < 0) {r = ibs0finder.posvec_que.back()->back();}
      double rcm = pgmap.bp2cm(r) - pgmap.bp2cm(pos);
      double lcm = pgmap.bp2cm(pos) - pgmap.bp2cm(l);
      double tcm = pgmap.bp2cm(r) - pgmap.bp2cm(l);
      rec->Add(rcm,lcm,tcm);
    }
    // Summarize the control pairs and output for this variant
    rec->Organize();
    int32_t eg = uni(rng);
// std::cout << rec->oldrec << '\t' << rec->ibs0pairs[eg].leftcm << '\t' << rec->ibs0pairs[eg].rightcm << '\t' << rec->ibs0pairs[eg].totcm << '\t' << rec->rank << '\t' << rec->mean << '\t' << rec->median << '\t' << rec->stdev << '\n';
    wf << rec->oldrec << '\t' << rec->ibs0pairs[eg].leftcm << '\t' << rec->ibs0pairs[eg].rightcm << '\t' << rec->ibs0pairs[eg].totcm << '\t' << rec->rank << '\t' << rec->mean << '\t' << rec->median << '\t' << rec->stdev << '\n';
    ndone ++;
  }

  wf.close();

  notice("Finish. Processed %d variants outof %d.", ndone, Doubts.size());

  return 0;
}

