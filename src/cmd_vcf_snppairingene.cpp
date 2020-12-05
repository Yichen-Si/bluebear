#include "cramore.h"
#include "bcf_filter_arg.h"
#include "bcf_ordered_reader.h"
// #include "bcf_ordered_writer.h"
#include "hts_utils.h"
#include "bp2cm.h"
#include <bits/stdc++.h>
// #include "ibs0.h"
#include <sstream>

struct snp {
  int32_t pos, ac;
  double cm;
  // std::vector<int32_t> gt;
  // std::map<int32_t, std::vector<int32_t> > carry;
  std::set<int32_t> carriers;
  int32_t max_impact = 0;
  std::string all_impact;
  // std::map<std::string, std::vector<std::string> > impact;
  // transcript -> (transcript type, snp function, impact category)

  snp(int32_t _p, int32_t _ac = 0) : pos(_p), ac(_ac) {}
  snp(int32_t _p, int32_t _ac, double _c) : pos(_p), ac(_ac), cm(_c) {}
  void add_gpos(double _c) {cm=_c;}
  void add_carrier(int32_t _id) {
    carriers.insert(_id);
  }
  void add_carrier_ct(int32_t _id, int32_t _g = 1) {
    carriers.insert(_id);
    ac += _g;
  }
  void set_ac(int32_t _c) {ac = _c;}
  void add_function(char* _ann, std::string & name,
               std::map<std::string, int32_t> & impact_category) {
    std::string fun(_ann);
    std::vector<std::string> ts;
    split(ts, ";", fun);
    std::stringstream ss;
    for (auto & v : ts) {
      std::vector<std::string> wd;
      split(wd, "|", v);
      if (wd[4] == name && wd[6] != "") {
        // impact[wd[6]] = std::vector<std::string>{wd[7],wd[2], wd[1]};
        const auto & ptr = impact_category.find(wd[2]);
        if ( ptr != impact_category.end() ) {
          max_impact = std::max(max_impact, ptr->second);
          ss<<wd[6]<<":"<<wd[1]<<"|"<<wd[2]<<";";
        }
      }
    }
    all_impact = ss.str();
    if (all_impact == "") {
      all_impact = ".";
    }
  }
};

int32_t cmdVcfSnpPairInGeneCount(int32_t argc, char** argv) {

  std::string inVcf, inMap, inGene, out;
  int32_t max_ac=10, min_ac=2, min_bp=-1, max_bp=INT_MAX;
  int32_t verbose = 10000;
  bool coding_only = 0;

  paramList pl;
  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input", NULL)
    LONG_STRING_PARAM("in-vcf",&inVcf, "Input VCF/BCF file")
    LONG_STRING_PARAM("map",&inMap, "Genetic map file")
    LONG_STRING_PARAM("gene",&inGene, "Gene region bed")

    LONG_PARAM_GROUP("Additional Options", NULL)
    LONG_INT_PARAM("max-ac",&max_ac, "Maximum sample allele count to consider")
    LONG_INT_PARAM("min-ac",&min_ac, "Maximum sample allele count to consider")
    LONG_INT_PARAM("coding-only",&coding_only, "If only include coding regions")
    LONG_INT_PARAM("min-inter",&min_bp, "Minimum distance between a pair of SNPs")
    LONG_INT_PARAM("max-inter",&max_bp, "Maximum distance between a pair of SNPs")
    LONG_INT_PARAM("verbose",&verbose, "Periodic report when scanning VCF")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output file")

  END_LONG_PARAMS();
  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  // bp2cmMap pgmap(inMap, " ");

  // Read gene region file
  std::vector<std::vector<std::string> > gene_info;
  std::ifstream ifs;
  std::string line, reg;
  std::vector<std::string> words;
  ifs.open(inGene, std::ifstream::in);
  if (!ifs.is_open()) {
    error("Gene region file cannot be opened");
  }
  while(std::getline(ifs, line)) {
    words.clear();
    split(words, "\t", line);
    reg = words[0] + ":" + words[1] + "-" + words[2];
    gene_info.push_back( std::vector<std::string>{reg, words[3], words[4], words[0]} );
  }
  ifs.close();

  std::map<std::string, int32_t> impact_category = { {"HIGH",4}, {"MODERATE",3}, {"LOW",2}, {"MODIFIER",1} };

  // output
  FILE * wf;
  wf = fopen(out.c_str(), "w");
  // std::ofstream wf(out.c_str(), std::ios::out);
  fprintf(wf, "%s\n", "CHR\tGeneID\tGeneName\tnSNP\tnSyn\tnMis\tnLoF\tNoShare\tShareOne\tShareMore\tNested\tMinorTwo\tMinorMore\tFunTwo\tFunMore");

for (auto & one_gene : gene_info) {

  std::vector<GenomeInterval> intervals;
  parse_intervals(intervals, "", one_gene[0]);
  BCFOrderedReader odr(inVcf, intervals);
  bcf1_t* iv = bcf_init();
  int32_t nsamples = bcf_hdr_nsamples(odr.hdr);

  // Store minor alleles
  std::vector<snp> snps;
  // Count snps in gene by function
  std::vector<int32_t> ct_fun(5, 0);
  std::vector<std::vector<int32_t> > ct_indiv;
  for (int32_t it = 0; it < 2; ++it) {
    std::vector<int32_t> tmp(nsamples,0);
    ct_indiv.push_back(tmp);
  }
  notice("Start region %s for gene %s", one_gene[0].c_str(), one_gene[1].c_str() );

  // Scan the gene region for rare alleles
  int32_t nsnps = 0;
  int32_t *p_gt = NULL;
  int32_t  n_gt = 0;
  int32_t *info_ac = NULL;
  int32_t  n_ac = 0;
  char  *ann = NULL;
  int32_t  n_ann = 0;
  int32_t  h1, h2;

  for (int32_t k = 0; odr.read(iv); ++k) {
    if (!bcf_is_snp(iv)) {continue;}
    if (bcf_get_info_int32(odr.hdr, iv, "AC", &info_ac, &n_ac) < 0) {continue;}
    int32_t ac = info_ac[0];
    int32_t minor = 1;
    if (ac > nsamples) {ac = 2*nsamples - ac; minor = 0;}
    if (ac < min_ac || ac > max_ac ) {continue;}

    if ( k % verbose == 0 )
      notice("Processing %d variants at %s:%d. %d kept", k, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1, nsnps);

    if ( bcf_get_genotypes(odr.hdr, iv, &p_gt, &n_gt) < 0 ) {
      error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);
    }

    snp rec(iv->pos+1, ac);
    if (bcf_get_info_string(odr.hdr, iv, "ANN", &ann, &n_ann) >= 0) {
      rec.add_function(ann, one_gene[1], impact_category);
    } else {
      rec.max_impact = 0;
      rec.all_impact = ".";
    }
    if (coding_only && rec.max_impact < 2) {continue;}

    ct_fun[rec.max_impact] += 1;
    int32_t anno_fun = (int32_t) (rec.max_impact > 2);

    for (int32_t it = 0; it < nsamples; ++it) {
      h1 = bcf_gt_allele(p_gt[2*it]);
      h2 = bcf_gt_allele(p_gt[2*it+1]);
      if (h1 == minor || h2 == minor) {
        rec.add_carrier(it);
        ct_indiv[anno_fun][it] += 1;
      }
    }
    snps.emplace_back(rec);
    nsnps += 1;
  }
  notice("Finished reading genotypes, recorded %d rare variants", nsnps);
  if (nsnps < 2) {
    continue;
  }

  std::map<std::string, int64_t> ct_pair = {
    {"NoSharing",0}, {"Share_one", 0}, {"Share_more", 0}, {"Nested", 0}
  };
  std::map<std::string, int32_t> ct_carry = {
    {"Minor_two",0}, {"Minor_more", 0}, {"Fun_two", 0}, {"Fun_more", 0}
  };
  for (int32_t it = 0; it < nsamples; ++it) {
    int32_t tot_minor = ct_indiv[0][it] + ct_indiv[1][it];
    if (tot_minor == 2) {
      ct_carry["Minor_two"] += 1;
    } else if (tot_minor > 2) {
      ct_carry["Minor_more"] += 1;
    }
    if (ct_indiv[1][it] == 2) {
      ct_carry["Fun_two"] += 1;
    } else if (ct_indiv[1][it] > 2) {
      ct_carry["Fun_more"] += 1;
    }
  }

  // Get pairwise genotype information
  for (uint32_t i = 1; i < snps.size(); i++) {
    for (uint32_t j = 0; j < i; j++) {
      if ( snps[i].pos - snps[j].pos > max_bp ) {continue;}
      if ( snps[i].pos - snps[j].pos < min_bp ) {break;}

      std::set<int32_t> id_both;
      std::set_intersection(snps[i].carriers.begin(), snps[i].carriers.end(),
                            snps[j].carriers.begin(), snps[j].carriers.end(),
                            std::inserter(id_both, id_both.begin()));
      if (id_both.size() > 0) {
        if (id_both.size() == 1) {
          ct_pair["Share_one"] += 1;
        } else {
          ct_pair["Share_more"] += 1;
        }
        if (id_both.size() < snps[i].carriers.size() &&
            id_both.size() < snps[j].carriers.size()) {
          continue;
        }
        std::vector<int32_t> gtpair(3,0);
        for (auto & u : snps[j].carriers) {
          if (id_both.find(u) == id_both.end()) {
            gtpair[1] += 1;
          }
        }
        for (auto & u : snps[i].carriers) {
          if (id_both.find(u) == id_both.end()) {
            gtpair[2] += 1;
          }
        }
        if (gtpair[1] == 0 || gtpair[1] == 0) {
          ct_pair["Nested"] += 1;
        }
      }
      else {
        ct_pair["NoSharing"] += 1;
      }
    }
  }
  bcf_destroy(iv);

  fprintf(wf, "%s\t%s\t%s\t", one_gene[3].c_str(), one_gene[1].c_str(), one_gene[2].c_str() );
  fprintf(wf, "%d\t%d\t%d\t%d\t", nsnps, ct_fun[2], ct_fun[3], ct_fun[4]);
  fprintf(wf, "%d\t%d\t%d\t%d\t", ct_pair["NoSharing"],ct_pair["Share_one"],ct_pair["Share_more"],ct_pair["Nested"]);
  fprintf(wf, "%d\t%d\t%d\t%d\n", ct_carry["Minor_two"], ct_carry["Minor_more"], ct_carry["Fun_two"], ct_carry["Fun_more"]);
}
  fclose(wf);
  return 0;
}











int32_t cmdVcfSnpPairInGeneSummary(int32_t argc, char** argv) {

  std::string inVcf, inMap, inGene, inSample, out;
  int32_t max_ac, min_ac, min_bp=-1, max_bp=INT_MAX;
  int32_t nsamples = 0;
  int32_t verbose = 10000;
  std::vector<int32_t> acThres;
  std::vector<double> afThres;
  std::vector<double> cmThres{0.01, 0.05, 0.1};
  std::vector<std::string> cmLab{"<0.01","0.01-0.05","0.05-0.1",">0.1"};
  std::map<std::string, int32_t> impact_category = { {"HIGH",4}, {"MODERATE",3}, {"LOW",2}, {"MODIFIER",1} };
  bool coding_only = 1;
  int32_t n_cm = cmThres.size();

  paramList pl;
  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input", NULL)
    LONG_STRING_PARAM("in-vcf",&inVcf, "Input VCF/BCF file")
    LONG_STRING_PARAM("map",&inMap, "Genetic map file")
    LONG_STRING_PARAM("gene",&inGene, "Gene region bed")

    LONG_PARAM_GROUP("Additional Options", NULL)
    // LONG_INT_PARAM("max-ac",&max_ac, "Maximum sample allele count to consider")
    // LONG_INT_PARAM("min-ac",&min_ac, "Maximum sample allele count to consider")
    LONG_INT_PARAM("sample-size",&nsamples, "Sample size")
    LONG_STRING_PARAM("subset-sample",&inSample, "Index of subset of samples")
    LONG_INT_PARAM("coding-only",&coding_only, "If only include coding regions")
    LONG_MULTI_INT_PARAM("ac",&acThres,"Allele count threshold to count rare/common variants")
    LONG_MULTI_DOUBLE_PARAM("af",&afThres,"Allele frequency threshold to count rare/common variants")
    LONG_INT_PARAM("min-inter",&min_bp, "Minimum distance between a pair of SNPs")
    LONG_INT_PARAM("max-inter",&max_bp, "Maximum distance between a pair of SNPs")
    LONG_INT_PARAM("verbose",&verbose, "Periodic report when scanning VCF")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output file")

  END_LONG_PARAMS();
  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  if ( inVcf.empty() || out.empty() || inMap.empty() || inGene.empty() || acThres.size() == 0 ) {
    error("[E:%s:%d %s] --in-vcf, --gene, --map, --ac, --out are required parameters",__FILE__,__LINE__,__FUNCTION__);
  }

  if (nsamples == 0 && afThres.size() == 0) {
    error("[E:%s:%d %s] If afThres is specified, sample-size is required",__FILE__,__LINE__,__FUNCTION__);
  }

  std::vector<int32_t> ac_lower, ac_upper;
  std::sort(acThres.begin(), acThres.end());
  for(int32_t i=0; i < (int32_t)acThres.size()-1; ++i) {
    if (acThres[i] > 0) {
      ac_lower.push_back(acThres[i]);
      ac_upper.push_back(acThres[i+1]-1);
    }
  }
  min_ac = ac_lower[0];
  max_ac = ac_upper.back();
  if (afThres.size() > 0) {
    ac_lower.push_back(acThres.back());
    std::sort(afThres.begin(), afThres.end());
    for (int32_t i=0; i < (int32_t)afThres.size()-1; ++i) {
      int32_t tmp_ac = (int32_t) (afThres[i] * nsamples * 2);
      if ( tmp_ac > max_ac ) {
        ac_upper.push_back(tmp_ac-1);
        ac_lower.push_back(tmp_ac);
      }
    }
    max_ac = (int32_t) (afThres.back() * nsamples * 2);
    ac_upper.push_back(max_ac);
  }

  std::cout << ac_lower.size() << '\t' << ac_upper.size() << '\n';
  for (uint32_t i = 0; i < ac_lower.size(); ++ i) {
    std::cout << ac_lower[i] << '\t' << ac_upper[i] << '\n';
  }

  bp2cmMap pgmap(inMap, " ");

  // Read gene region file
  std::vector<std::vector<std::string> > gene_info;
  std::ifstream ifs;
  std::string line, reg;
  std::vector<std::string> words;
  ifs.open(inGene, std::ifstream::in);
  if (!ifs.is_open()) {
    error("Gene region file cannot be opened");
  }
  while(std::getline(ifs, line)) {
    words.clear();
    split(words, "\t", line);
    reg = words[0] + ":" + words[1] + "-" + words[2];
    gene_info.push_back( std::vector<std::string>{reg, words[3], words[4], words[0]} );
  }
  ifs.close();

  std::vector<int32_t> sample_id;
  if (! inSample.empty()) {
    ifs.open(inSample, std::ifstream::in);
    if (!ifs.is_open()) {
      error("Sample id file cannot be opened");
    }
    while(std::getline(ifs, line)) {
      words.clear();
      split(words, "\t", line);
      sample_id.push_back(std::stoi(words[0]));
    }
    ifs.close();
  }


  std::map<int32_t, std::map<int32_t, std::map<int32_t, std::vector<double> > > > ct_pair;
  for (uint32_t i = 0; i < ac_lower.size(); ++i) {
    for (uint32_t j = i; j < ac_lower.size(); ++j) {
      for (int32_t k = 0; k < n_cm + 1; ++k) {
        std::vector<double> tmp(4, 0.0);
        ct_pair[i][j][k] = tmp;
      }
    }
  }

for (auto & one_gene : gene_info) {

  std::vector<GenomeInterval> intervals;
  parse_intervals(intervals, "", one_gene[0]);
  BCFOrderedReader odr(inVcf, intervals);
  bcf1_t* iv = bcf_init();
  int32_t n_samples = bcf_hdr_nsamples(odr.hdr);
  if (sample_id.size() == 0) {
    for (int32_t it = 0; it < n_samples; ++it) {
      sample_id.push_back(it);
    }
    nsamples = n_samples;
  }

  // Store minor alleles
  std::vector<snp> snps;
  // Count snps in gene by function
  std::vector<int32_t> ct_fun(5, 0);
  std::vector<std::vector<int32_t> > ct_indiv;
  for (int32_t it = 0; it < 2; ++it) {
    std::vector<int32_t> tmp(nsamples,0);
    ct_indiv.push_back(tmp);
  }
  notice("Start region %s for gene %s", one_gene[0].c_str(), one_gene[1].c_str() );

  // Scan the gene region for rare alleles
  int32_t nsnps = 0;
  int32_t *p_gt = NULL;
  int32_t  n_gt = 0;
  int32_t *info_ac = NULL;
  int32_t  n_ac = 0;
  int32_t  h1, h2;
  char  *ann = NULL;
  int32_t  n_ann = 0;

  for (int32_t k = 0; odr.read(iv); ++k) {
    if (!bcf_is_snp(iv)) {continue;}
    if (bcf_get_info_int32(odr.hdr, iv, "AC", &info_ac, &n_ac) < 0) {continue;}
    int32_t ac = info_ac[0];
    int32_t minor = 1;
    if (ac > nsamples) {ac = 2*nsamples - ac; minor = 0;}
    if (ac < min_ac ) {continue;}

    if ( k % verbose == 0 )
      notice("Processing %d variants at %s:%d. %d kept", k, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1, nsnps);

    if ( bcf_get_genotypes(odr.hdr, iv, &p_gt, &n_gt) < 0 ) {
      error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);
    }

    double gp = pgmap.bp2cm(iv->pos);
    snp rec(iv->pos+1, 0, gp);

    if (bcf_get_info_string(odr.hdr, iv, "ANN", &ann, &n_ann) >= 0) {
      rec.add_function(ann, one_gene[1], impact_category);
    } else {
      continue;
    }
    if (coding_only && rec.max_impact < 2) {continue;}

    for ( auto & it : sample_id ) {
      h1 = bcf_gt_allele(p_gt[2*it]);
      h2 = bcf_gt_allele(p_gt[2*it+1]);
      if (h1 == minor || h2 == minor) {
        rec.add_carrier_ct(it, (int32_t) (h1==minor)+ (int32_t) (h2==minor) );
      }
    }
    if ( rec.ac < min_ac || rec.ac > max_ac ) {continue;}
    int32_t ac_cat = -1;
    for (int32_t i = 0; i < (int32_t) ac_lower.size(); ++i ) {
      if (rec.ac >= ac_lower[i] && rec.ac <= ac_upper[i]) {
        ac_cat = i;
        break;
      }
    }
    if (ac_cat < 0) {continue;}
    rec.set_ac(ac_cat);
    snps.emplace_back(rec);
    nsnps += 1;
  }
  notice("Finished reading genotypes, recorded %d rare variants", nsnps);
  if (nsnps < 2) {
    bcf_destroy(iv);
    continue;
  }

  // Get pairwise information
  for (uint32_t i = 0; i < snps.size() - 1; i++) {
    for (uint32_t j = i+1; j < snps.size(); j++) {
      if ( snps[j].pos - snps[i].pos > max_bp ) {break;}
      if ( snps[j].pos - snps[i].pos < min_bp ) {continue;}

      std::set<int32_t> id_both;
      std::set_intersection(snps[i].carriers.begin(), snps[i].carriers.end(),
                            snps[j].carriers.begin(), snps[j].carriers.end(),
                            std::inserter(id_both, id_both.begin()));
      int32_t ac1 = std::min(snps[i].ac, snps[j].ac);
      int32_t ac2 = std::max(snps[i].ac, snps[j].ac);
      double g_dist = snps[j].cm - snps[i].cm;
      int32_t k_cm = 0;
      while ( k_cm < n_cm && g_dist > cmThres[k_cm] ) {
        k_cm += 1;
      }

      ct_pair[ac1][ac2][k_cm][0] += 1;
      if (id_both.size() > 0) {
        if (id_both.size() == 1) {
          ct_pair[ac1][ac2][k_cm][1] += 1;
        } else {
          ct_pair[ac1][ac2][k_cm][2] += 1;
        }
        if (id_both.size() == snps[i].carriers.size() ||
            id_both.size() == snps[j].carriers.size()) {
          std::vector<int32_t> gtpair(3,0);
          for (auto & u : snps[j].carriers) {
            if (id_both.find(u) == id_both.end()) {
              gtpair[1] += 1;
            }
          }
          for (auto & u : snps[i].carriers) {
            if (id_both.find(u) == id_both.end()) {
              gtpair[2] += 1;
            }
          }
          if (gtpair[1] == 0 || gtpair[1] == 0) {
            ct_pair[ac1][ac2][k_cm][3] += 1;
          }
        }
      }
    }
  }
  bcf_destroy(iv);
}

  // output
  FILE * wf;
  wf = fopen(out.c_str(), "w");
  // std::ofstream wf(out.c_str(), std::ios::out);
  fprintf(wf, "%s\n", "AC1_bin\tAC2_bin\tcM_bin\tcM_range\tTotalPair\tShareOne\tShareMore\tNested");

  for ( auto p1 : ct_pair ) {
    for ( auto p2 : p1.second ) {
      for ( auto p3 : p2.second ) {
        fprintf(wf, "%d\t%d\t%d\t%s\t", p1.first, p2.first, p3.first, cmLab[p3.first].c_str() );
        fprintf(wf, "%.2f\t%.2f\t%.2f\t%.2f\n", p3.second[0], p3.second[1], p3.second[2], p3.second[3]  );
      }
    }
  }

  fclose(wf);
  return 0;
}


















