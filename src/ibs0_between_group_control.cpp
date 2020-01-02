#include "cramore.h"
#include "utils.h"
#include "bcf_filter_arg.h"
#include "bcf_ordered_reader.h"
#include "compact_matrix.h"
#include <iomanip>
#include "bp2cm.h"
#include "ibs0.h"
#include <numeric>
#include <random>


struct RareVar {

  int32_t AC, pos;
  std::vector<int32_t> carriers;
  RareVar(int32_t _ac, int32_t _p, std::vector<int32_t> & vec) : AC(_ac), pos(_p) {
    for (auto & v : vec)
      carriers.push_back(v);
  }
  RareVar(int32_t _ac, int32_t _p) : AC(_ac), pos(_p) {}
  void AddID(int32_t _id) {carriers.push_back(_id);}
};

struct RareFocal {

  int32_t AC, pos;
  std::vector< std::pair<RareVar*,RareVar*> > pair_list;
  RareFocal(int32_t _ac, int32_t _p) : AC(_ac), pos(_p) {}
  void AddPair(RareVar* ptr1, RareVar* ptr2) {
    pair_list.push_back( std::pair<RareVar*,RareVar*>(ptr1, ptr2) );
  }
};

template<class T>
std::vector<double> VecSummary(std::vector<T>& vec) {
  // Basic summary statistics
  std::sort(vec.begin(), vec.end());
  int32_t n = vec.size();
  double median;
  if (n < 2)
    median = vec[0];
  else
    median = vec[n/2-1];
  if (n % 2 == 0) {
    median = (median + vec[n/2]) / 2;
  }
  double mean = std::accumulate(vec.begin(), vec.end(), 0.0) / n;
  double min = *std::min_element(vec.begin(), vec.end());
  double max = *std::max_element(vec.begin(), vec.end());
  std::vector<double> res{mean, median, min, max};
  return res;
}

// goal -- For a focal rare variant, find a pair of vairanst nearby with allele counts summed to its AC, then get stats on between-cluster distance
int32_t IBS0AddPMDistr(int32_t argc, char** argv) {
  std::string largeVcf, smallVcf, inMap, reg;
  std::string out, outf;
  int32_t cst = -1, ced = -1;
  int32_t min_hom_gts = 1;
  int32_t verbose = 5000;
  int32_t min_variant = 1;
  int32_t start, end;
  int32_t ck_len = 500000;
  double  cm_limit = 2.0;
  int32_t MAXAC = 10;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Files", NULL)
    LONG_STRING_PARAM("in-vcf-full",&largeVcf, "Input VCF/BCF file, including singletons")
    LONG_STRING_PARAM("in-vcf-small",&smallVcf, "Input VCF/BCF file, including variants for identify IBS0")
    LONG_STRING_PARAM("region",&reg,"Genomic region to focus on")
    LONG_STRING_PARAM("map",&inMap, "Map file for genetic distance")

    LONG_PARAM_GROUP("Additional Options", NULL)
    LONG_INT_PARAM("min-hom",&min_hom_gts, "Minimum number of homozygous genotypes to be counted for IBS0")
    LONG_INT_PARAM("min-variant",&min_variant, "Minimum number of variants to present to have output file")
    LONG_DOUBLE_PARAM("max-ibs0",&cm_limit, "Maximal ibs0 in cM to search for.")
    LONG_INT_PARAM("max-ac",&MAXAC, "Maximal allele count to consider as focal rare variants")

    LONG_INT_PARAM("centromere-st",&cst, "Start position of centromere")
    LONG_INT_PARAM("centromere-ed",&ced, "End position of centromere")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output file prefix")
    LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")

  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  // sanity check of input arguments
  if ( largeVcf.empty() || reg.empty() || inMap.empty() || out.empty() ) {
    error("[E:%s:%d %s] --in-vcf-full, --region, --inMap, --out are required parameters",__FILE__,__LINE__,__FUNCTION__);
  }
  if (smallVcf.empty())
    smallVcf = largeVcf;

  std::random_device rd;
  std::mt19937 rng(rd());

  // Genetic map
  bp2cmMap pgmap(inMap, " ", cst, ced);

  // Output
  std::ofstream wf;
  outf = out + "_pm_subdivide_null_" + reg + ".txt";
  wf.open(outf, std::ios::trunc);

  // Region to process
  std::vector<std::string> v;
  split(v, ":-", reg);
  if (!str2int32(v[1], start) || !str2int32(v[2], end)) {
    error("Invalid region.");
  }

  // bcf reader
  std::vector<GenomeInterval> intervals;
  parse_intervals(intervals, "", reg);
  BCFOrderedReader odr(largeVcf, intervals);
  bcf1_t* iv = bcf_init();

  int32_t nsamples = bcf_hdr_nsamples(odr.hdr);
  notice("Identifying %d samples from Vcf file", nsamples);
  if ( nsamples == 0 )
    error("FATAL ERROR: The VCF does not have any samples with genotypes");

  // Store all rare variants
  std::vector<RareVar*> AllRare;
  // First pass: list all rare variants & their carriers
  int32_t* p_gt = NULL;
  int32_t  n_gt = 0;
  int32_t *info_ac = NULL;
  int32_t  n_ac = 0;
  for(int32_t k=0; odr.read(iv); ++k) {
    // Notice: Assuming no need for QC filter
    if (bcf_get_info_int32(odr.hdr, iv, "AC", &info_ac, &n_ac) < 0) {continue;}
    if (info_ac[0] < 1 || info_ac[0] > MAXAC) {continue;}

    // Extract genotype and identify carriers
    if ( bcf_get_genotypes(odr.hdr, iv, &p_gt, &n_gt) < 0 ) {
      error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);
    }
    RareVar* ptr = new RareVar(info_ac[0], iv->pos+1);
    for (int32_t it = 0; it < nsamples; ++it) {
      if (bcf_gt_allele(p_gt[2*it]) > 0 || bcf_gt_allele(p_gt[2*it+1]) > 0) {
        ptr->AddID(it);
      }
    }
    if ( (int32_t) (ptr->carriers).size() != info_ac[0]) // Ignore
      delete ptr;
    else
      AllRare.push_back(ptr);
  }
  int32_t totrare = AllRare.size();
  notice("Found %d rare variants with AC <= %d.", totrare, MAXAC);

  // Identify nearby pairs of vairanst
  // Currently: find a pair of flanking variants
  std::vector<RareFocal> FocalSites;
  std::bernoulli_distribution ber(0.5);
  for (int32_t i = 1; i < totrare-1; ++i) {
    if (AllRare[i]->AC >= 2) {
      int32_t ac = AllRare[i]->AC;
      RareFocal tmp(ac, AllRare[i]->pos);
      // FocalSites.push_back(RareFocal(ac, AllRare[i]->pos));
      for (int32_t m = 1; m <= ac/2; ++m) {
        RareVar *ptr1, *ptr2;
        if (ber(rng)) {
          int32_t j = i - 1;
          while (AllRare[j]->AC != m && j > 0)
            j--;
          if (AllRare[j]->AC != m) {
            continue;
          }
          ptr1 = AllRare[j];
          j = i + 1;
          while (AllRare[j]->AC != ac - m && j < totrare-1)
            j++;
          if (AllRare[j]->AC != ac - m) {
            continue;
          }
          ptr2 = AllRare[j];
        } else {
          int32_t j = i + 1;
          while (AllRare[j]->AC != m && j < totrare-1)
            j++;
          if (AllRare[j]->AC != m) {
            continue;
          }
          ptr1 = AllRare[j];
          j = i - 1;
          while (AllRare[j]->AC != ac - m && j > 0)
            j--;
          if (AllRare[j]->AC != ac - m) {
            continue;
          }
          ptr2 = AllRare[j];
        }
        bool flag = 0;
        for (auto & c : ptr1->carriers) {
          if (std::find((ptr2->carriers).begin(), (ptr2->carriers).end(), c) != (ptr2->carriers).end()) {
            FocalSites.pop_back();
            flag = 1;
            break;
          }
        } // Need to check the carrier sets do not overlap
        if (flag)
          continue;
        FocalSites.push_back(tmp);
        FocalSites.back().AddPair(ptr1, ptr2);
      }
    }
  }
  int32_t totfocal = (int32_t) FocalSites.size();
  notice("Found %d focal rare variants with 2 <= AC <= %d", totfocal, MAXAC);

  // IBS0 lookup
  IBS0lookup ibs0finder(smallVcf, reg, pgmap, cm_limit, ck_len, 1);
  notice("Set up IBS0 finder");

  // Find pairwise ibs0
  int32_t ndone = 0;
  for (auto & v : FocalSites) {
    if (ndone % verbose == 0) {
      notice("Processed %d variants.", ndone);
    }
    for (auto & w : v.pair_list) {
      auto ptr1 = w.first;
      auto ptr2 = w.second;
      std::vector<int32_t> tvec;
      for (auto id1 : ptr1->carriers) {
        for (auto id2 : ptr2->carriers) {
          int32_t r = ibs0finder.FindIBS0(id1,id2,v.pos,0);
          int32_t l = ibs0finder.FindIBS0(id1,id2,v.pos,1);
          if (l < 0) {l = ibs0finder.start_que[0];}
          if (r < 0) {r = ibs0finder.posvec_que.back()->back();}
          int32_t offset = pgmap.overlap_centro(l,r);
          tvec.push_back(r-l-offset);
        }
      }
      std::vector<double> svec1 = VecSummary<int32_t>(tvec);
      // Random pairs for comparison
      std::set<int32_t> chosen1, chosen2;
      NchooseK(nsamples, ptr1->AC, chosen1, rng, 0);
      NchooseK(nsamples, ptr2->AC, chosen2, rng, 0);
      std::vector<int32_t> nvec;
      for (auto id1 : chosen1) {
        for (auto id2 : chosen2) {
          int32_t r = ibs0finder.FindIBS0(id1,id2,v.pos,0);
          int32_t l = ibs0finder.FindIBS0(id1,id2,v.pos,1);
          if (l < 0) {l = ibs0finder.start_que[0];}
          if (r < 0) {r = ibs0finder.posvec_que.back()->back();}
          int32_t offset = pgmap.overlap_centro(l,r);
          nvec.push_back(r-l-offset);
        }
      }
      std::vector<double> svec2 = VecSummary<int32_t>(nvec);
      // output: pos, ac1, ac2, (mean, med, min, max) x 2
      std::stringstream ss;
      ss << v.pos << '\t' << ptr1->AC << '\t' << ptr2->AC;
      for (auto v : svec1)
        ss << '\t' << std::round(v);
      for (auto v : svec2)
        ss << '\t' << std::round(v);
      ss << '\n';
      wf << ss.str();
    }
    ndone++;
  }
  wf.close();
  notice("Finish. Processed %d variants.", ndone);

  return 0;
}





