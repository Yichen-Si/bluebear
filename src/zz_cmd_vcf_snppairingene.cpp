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
  int32_t pos, ac=0;
  std::vector<int32_t> gt;
  std::map<int32_t, std::vector<int32_t> > carry;
  std::set<int32_t> carriers;
  int32_t max_impact = 0;
  std::string all_impact;
  // std::map<std::string, std::vector<std::string> > impact;
  // transcript -> (transcript type, snp function, impact category)

  snp(int32_t _p) : pos(_p) {
    for (int32_t i = 0; i < 3; ++i)
      gt.push_back(0);
  }
  void add_carrier(int32_t _id, int32_t _h1, int32_t _h2) {
    int32_t x = _h1+_h2;
    std::vector<int32_t> indiv{x, _h1, _h2};
    carry[_id] = indiv;
    ac+=x;
    gt[x] += 1;
    carriers.insert(_id);
  }
  void add_function(char* _ann, std::string name,
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
  int32_t max_ac=10, min_bp=-1, max_bp=INT_MAX;
  int32_t verbose = 10000;

  paramList pl;
  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input", NULL)
    LONG_STRING_PARAM("in-vcf",&inVcf, "Input VCF/BCF file")
    LONG_STRING_PARAM("map",&inMap, "Genetic map file")
    LONG_STRING_PARAM("gene",&inGene, "Gene region bed")

    LONG_PARAM_GROUP("Additional Options", NULL)
    LONG_INT_PARAM("max-ac",&max_ac, "Maximum sample allele count to consider")
    LONG_INT_PARAM("min-inter",&min_bp, "Minimum distance between a pair of SNPs")
    LONG_INT_PARAM("max-inter",&max_bp, "Maximum distance between a pair of SNPs")
    LONG_INT_PARAM("verbose",&verbose, "Periodic report when scanning VCF")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output file")

  END_LONG_PARAMS();
  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

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

  std::map<std::string, int32_t> impact_category = { {"HIGH",4}, {"MODERATE",3}, {"LOW",2}, {"MODIFIER",1} };

  // output
  FILE * wf;
  wf = fopen(out.c_str(), "w");
  // std::ofstream wf(out.c_str(), std::ios::out);
  fprintf(wf, "%s\n", "CHR\tGeneID\tGeneName\tPOS1\tPOS2\tDist_bp\tDist_cM\tAC1\tAC2\tGT10\tGT01\tGT11\tHTcis\tHTtrans\tAnno1\tAnno2\tMaxImpact1\tMaxImpact2");

for (auto & one_gene : gene_info) {

  notice("Start region %s for gene %s", one_gene[0].c_str(), one_gene[1].c_str() );
  std::vector<GenomeInterval> intervals;
  parse_intervals(intervals, "", one_gene[0]);
  BCFOrderedReader odr(inVcf, intervals);
  bcf1_t* iv = bcf_init();
  int32_t nsamples = bcf_hdr_nsamples(odr.hdr);

  // Store minor alleles
  std::vector<snp> snps;

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
    if (ac < 1 || ac > max_ac ) {continue;}

    if ( k % verbose == 0 )
      notice("Processing %d variants at %s:%d. %d kept", k, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1, nsnps);

    if ( bcf_get_genotypes(odr.hdr, iv, &p_gt, &n_gt) < 0 ) {
      error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);
    }

    snp rec(iv->pos+1);
    if (bcf_get_info_string(odr.hdr, iv, "ANN", &ann, &n_ann) >= 0) {
      rec.add_function(ann, one_gene[1], impact_category);
    }
    if (minor == 1) {
      for (int32_t it = 0; it < nsamples; ++it) {
        h1 = bcf_gt_allele(p_gt[2*it]);
        h2 = bcf_gt_allele(p_gt[2*it+1]);
        if (h1 == 1 || h2 == 1) {
          rec.add_carrier(it, h1, h2);
        }
      }
    } else {
      for (int32_t it = 0; it < nsamples; ++it) {
        h1 = 1-bcf_gt_allele(p_gt[2*it]);
        h2 = 1-bcf_gt_allele(p_gt[2*it+1]);
        if (h1 == 1 || h2 == 1) {
          rec.add_carrier(it, h1, h2);
        }
      }
    }

    snps.push_back(rec);
    nsnps += 1;
  }
  notice("Finished reading genotypes, recorded %d rare variants", nsnps);
  if (nsnps < 2) {
    continue;
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
        std::vector<int32_t> gtpair(3,0);
        int32_t cis = 0, trans = 0;
        int32_t dist_dp = snps[i].pos - snps[j].pos;
        double  dist_cm = pgmap.bpinterval2cm(snps[j].pos, snps[i].pos);
        for (auto & u : snps[j].carry) {
          if (id_both.find(u.first) == id_both.end()) {
            gtpair[1] += 1;
          }
        }
        for (auto & u : snps[i].carry) {
          if (id_both.find(u.first) == id_both.end()) {
            gtpair[2] += 1;
          }
        }
        for (auto & u : id_both) { // Count genotype config
          if ( snps[i].carry[u][0]==1 && snps[j].carry[u][0]==1 ) {
            gtpair[0]+=1;
            if ( snps[i].carry[u][1]==snps[j].carry[u][1] ) {cis += 1;}
            else {trans += 1;}
          }
        }
        fprintf(wf, "%s\t%s\t%s\t%d\t%d\t%d\t%.5f\t%d\t%d\t", one_gene[3].c_str(), one_gene[1].c_str(), one_gene[2].c_str(), snps[j].pos, snps[i].pos, dist_dp, dist_cm, snps[j].ac, snps[i].ac );
        fprintf(wf, "%d\t%d\t%d\t%d\t%d\t", gtpair[1], gtpair[2], gtpair[0], cis, trans);
        fprintf(wf, "%s\t%s\t%d\t%d\n", snps[j].all_impact.c_str(),snps[i].all_impact.c_str(),snps[j].max_impact, snps[i].max_impact);
      }
    }
  }
  bcf_destroy(iv);
}

  fclose(wf);
  return 0;
}





























