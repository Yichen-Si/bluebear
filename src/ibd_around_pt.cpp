#include "cramore.h"
#include "bcf_filter_arg.h"
#include "bcf_ordered_reader.h"
#include "compact_matrix.h"

#include "ibs0.h"
#include "pbwt_build.h"
#include "bp2cm.h"

#include <map>

struct refhap {
  int32_t id, hapid;
  int32_t hap, ibs;
  bool trunc = false;
};

struct switchpt {
  int32_t pos;
  int32_t id;
  std::map<int32_t, refhap*> reflist;
};

// Input bcf file & a sorted list of positions
// Matched haplotype length v.s. no-ibs0 length

int32_t HapibdVSnoibs0(int32_t argc, char** argv) {
  std::string inVcf, inMap, path_pbwt, inPos, out, reg, chrom;
  std::string outf,outf_s;
  int32_t detailed = 0;
  int32_t verbose = 1000;
  int32_t min_variant = 1;
  int32_t min_hom_gts = 1;
  int32_t min_ac = 1; // If we want to ignore rare variants in haplotype matching
  int32_t nsamples=0, M=0;
  int32_t ibs0_chunk = 1000000;  // process suffix pbwt & ibs0 by chunk
  int32_t pbwt_chunk = 1000000;
  int32_t delta = 2000000;  // Maximal no-ibs0 to look for on each side

  int32_t* p_gt = NULL;
  int32_t  n_gt = 0;

  paramList pl;
  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input", NULL)
    LONG_STRING_PARAM("in-vcf",&inVcf, "Input VCF/BCF file")
    LONG_STRING_PARAM("in-pos",&inPos, "Input list of position")
    LONG_STRING_PARAM("region",&reg,"Genomic region to focus on")
    LONG_STRING_PARAM("chr",&chrom,"Chromosome to process (if region is not provided)")
    LONG_STRING_PARAM("map",&inMap, "Map file for genetic distance")
    LONG_STRING_PARAM("pbwt-path",&path_pbwt,"Pre-computed snapshot of suffix pbwt")

    LONG_PARAM_GROUP("Additional Options", NULL)
    LONG_INT_PARAM("delta",&delta, "Maximal no-ibs0 to look for on each side")
    LONG_INT_PARAM("min-ac",&min_ac, "Minimal AC to be considered")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output file prefix")
    LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")
    LONG_INT_PARAM("detailed",&detailed,"If detailed info for each pair is to be written (0/1)")

  END_LONG_PARAMS();
  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  if (inVcf.empty() || inPos.empty() || inMap.empty() || path_pbwt.empty() || out.empty()) {
    error("[E:%s:%d %s] --in-vcf, --in-pos, --map, --pbwt-path, --out are required parameters",__FILE__,__LINE__,__FUNCTION__);
  }

  // Read the position - id file // Assume sorted by position
  std::fstream rf(inPos,std::ios::in);
  std::vector<std::vector<int32_t> > poslist;
  int32_t pos, id;
  int32_t start, end;
  while ( rf >> pos >> id ) {
    poslist.push_back(std::vector<int32_t> {pos,id});
  }
  notice("Read %d position-id pairs. From %d to %d.", poslist.size(), poslist[0][0], poslist.back()[0]);

  // Translate between genetic & physical position.
  bp2cmMap pgmap(inMap, " ");

  // Output
  outf = out + "_long.list";
  outf_s = out + "_short.list";
  FILE* wf;
  FILE* wfs;
  wfs = fopen(outf_s.c_str(), "w");
  if (detailed) {
    wf = fopen(outf.c_str(), "w");
  }

  // get sample size & chrom region
  htsFile *fp = hts_open(inVcf.c_str(),"r");
  bcf_hdr_t *hdr = bcf_hdr_read(fp);
  nsamples = hdr->n[BCF_DT_SAMPLE];
  M = nsamples * 2;
  hts_close(fp);
  if (reg.empty()) {
    start = poslist[0][0];
    end   = poslist.back()[0];
  } else {
    std::vector<std::string> v;
    split(v, ":-", reg);
    if (!str2int32(v[1], start) || !str2int32(v[2], end)) {
      error("Invalid input region.");
    }
    chrom = v[0];
  }

  // Initialize forward pbwt
  pbwtCursor prepc(M, 1);

  // Can read saved check points, but build it from start for now
  reg = chrom + ":0-" + std::to_string(start-1);
  std::vector<GenomeInterval> intervals;
  parse_intervals(intervals, "", reg);
  BCFOrderedReader odr(inVcf, intervals);
  bcf1_t* iv = bcf_init();
  for (int32_t k=0; odr.read(iv); ++k) { // Read haplotypes
    if (bcf_get_genotypes(odr.hdr, iv, &p_gt, &n_gt) < 0) {
      error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);
    }
    bool y[M];
    int32_t ac = 0;
    for (int32_t i = 0; i <  nsamples; ++i) {
      int32_t g1 = p_gt[2*i];
      int32_t g2 = p_gt[2*i+1];
      if (bcf_gt_is_missing(g1)) {
        y[2*i] = 0;
      } else {
        y[2*i] = (bcf_gt_allele(g1) > 0);
      }
      if (bcf_gt_is_missing(g2)) {
        y[2*i+1] = 0;
      } else {
        y[2*i+1] = (bcf_gt_allele(g2) > 0);
      }
      ac += (y[2*i] + y[2*i+1]);
    } // Finish reading one SNP
    if (ac >= min_ac)
      prepc.ForwardsAD_prefix(y, iv->pos+1);
  } // Finish reading all haplotypes before the first position


  // Start finding ibd around input position-id pairs.
  int32_t st = start;
  int32_t ed = (start/pbwt_chunk + 1)*pbwt_chunk;

  // initialize ibs0 lookup
  reg = chrom + ":" + std::to_string(st) + "-" + std::to_string(ed);
  IBS0lookup ibs0finder(inVcf, reg, pgmap, delta, ibs0_chunk, min_hom_gts);

  uint32_t iter = 0;
  bool flag = 0;
  pos = poslist[iter][0];

  while (st < pgmap.maxpos) {
    // Read and build suffix pbwt by physical chunk
    // TODO: enable efficient random access &/ customizable chunk size

    reg = chrom + ":" + std::to_string(st) + "-" + std::to_string(ed);
    ibs0finder.Update(reg);
    // Read genotype matrix
    parse_intervals(intervals, "", reg);
    odr.jump_to_interval(intervals[0]);
    bcf_clear(iv);

    std::vector<bool*> gtmat;
    std::vector<int32_t> positions;

    for (int32_t k=0; odr.read(iv); ++k) { // Read haplotypes
      if (bcf_get_genotypes(odr.hdr, iv, &p_gt, &n_gt) < 0) {
        error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);
      }
      bool *y = new bool[M];
      int32_t ac = 0;
      for (int32_t i = 0; i <  nsamples; ++i) {
        int32_t g1 = p_gt[2*i];
        int32_t g2 = p_gt[2*i+1];
        if (bcf_gt_is_missing(g1)) {
          y[2*i] = 0;
        } else {
          y[2*i] = (bcf_gt_allele(g1) > 0);
        }
        if (bcf_gt_is_missing(g2)) {
          y[2*i+1] = 0;
        } else {
          y[2*i+1] = (bcf_gt_allele(g2) > 0);
        }
        ac += (y[2*i] + y[2*i+1]);
      } // Finish reading one SNP
      if (ac < min_ac)
        continue;
      gtmat.push_back(y);
      positions.push_back(iv->pos+1);
    } // Finish reading all haplotypes in this chunk
    int32_t N = positions.size();
    if ( N < min_variant ) {
      if (N > 0) {
        for (int32_t k = 0; k < N; ++ k)
          prepc.ForwardsAD_prefix(gtmat[k], positions[k]);
        for (auto& pt : gtmat)
          delete [] pt;
      }
      st = ed + 1;
      ed += pbwt_chunk;
      notice("Observed only %d informative markers. Skipping the region %s", N, reg.c_str());
      continue;
    } else {
      notice("Read %d markers, start processing %s.", N, reg.c_str());
    }

    // Read the last snp in this chunk. should be stored as a checkpoint
    std::string sreg = chrom+":"+std::to_string(ed-pbwt_chunk+1)+"-"+std::to_string(ed);
    std::string sfile, line;
    std::ifstream sf;
    std::vector<std::string> wvec;
    int32_t a[M], d[M];

    sfile = path_pbwt + "_" + sreg + "_suffix.amat";
    sf.open(sfile);
    std::getline(sf, line);
    sf.close();
    split(wvec, "\t", line);
    uint32_t OFFSET = wvec.size() - M;
    for (uint32_t i =  OFFSET; i < wvec.size(); ++i)
      a[i- OFFSET] = std::stoi(wvec[i]);

    sfile = path_pbwt + "_" + sreg + "_suffix.dmat";
    sf.open(sfile);
    std::getline(sf, line);
    sf.close();
    split(wvec, "\t", line);
    for (uint32_t i =  OFFSET; i < wvec.size(); ++i)
      d[i- OFFSET] = std::stoi(wvec[i]);

    pbwtCursor sufpc(M, a, d);

    // Build & store pbwt suffix matricies
    std::vector<int32_t*> dmat(N,NULL);
    std::vector<int32_t*> amat(N,NULL);
    for (int32_t i = 0; i < N; ++i) {
      dmat[i] = new int32_t[M];
      amat[i] = new int32_t[M];
    }
    for (int32_t k = N-1; k >= 0; --k) {
      sufpc.ForwardsAD_suffix(gtmat[k], positions[k]);
      memcpy(dmat[k], sufpc.d, M*sizeof(int32_t));
      memcpy(amat[k], sufpc.a, M*sizeof(int32_t));
    }

    // Build prefix pbwt
    int32_t id1, id2; // individual id
    int32_t i_s, j_s, i_p, j_p, s_i, s_j;
    int32_t l0, r0;
    int32_t dijp, dijs;
    int32_t k = 0;
    while (k <= N-1) {
      // Sorted upto and include position k
      prepc.ForwardsAD_prefix(gtmat[k], positions[k]);
      delete [] gtmat[k];
      if (positions[k] < pos) {
        k++;
        continue;
      }
      while (positions[k] > pos) {
        iter++;
        if (iter >= poslist.size())
          break;
        pos = poslist[iter][0];
      }
      if (iter >= poslist.size()) {
        break;
      }
      int32_t rvec[M]; // Map haplotype id to row in suffix matrix
      for (int32_t it = 0; it < M; ++it) {
        rvec[amat[k][it]] = it;
      }

      while (pos == positions[k]) {
        int32_t ibd = 0, ibs = 0, besthap = -1, bestibs = -1;
// std::cout << k << '\t' << iter << '\t' << positions[k] << '\t' << pos << '\n';
        id1 = poslist[iter][1];
        switchpt pt;
        pt.pos=pos;
        pt.id =id1;
        for (int32_t h = id1*2; h <= id1*2+1; ++h) {
          i_p = prepc.FindIndex(h);
          i_s = rvec[h];

          // Left
          if (i_p == M-1 || prepc.d[i_p] < prepc.d[i_p+1]) {
            j_p = i_p-1;
            dijp = prepc.d[i_p];
          } else {
            j_p = i_p+1;
            dijp = prepc.d[i_p+1];
          }
          id2 = prepc.a[j_p]/2;
          if (id2 != id1) {
            j_s = rvec[prepc.a[j_p]];
            refhap* lh = new refhap;
            lh->hapid = prepc.a[j_p];
            s_i = std::min(i_s,j_s);
            s_j = std::max(i_s,j_s);
            dijs = dmat[k][s_i+1];
            for (int32_t it = s_i+1; it <= s_j; ++it) {
              dijs = std::min(dijs,dmat[k][it]);
            }
            lh->hap = dijs - dijp - pgmap.overlap_centro(dijp,dijs);
  // std::cout << dijp << ' ' << dijs << ' ' << pgmap.overlap_centro(dijp,dijs) << '\n';
            lh->id=id2;
            r0 = ibs0finder.FindIBS0(id1,id2,pos,0);
            l0 = ibs0finder.FindIBS0(id1,id2,pos,1);
            lh->trunc = (r0 < 0 || l0 < 0);
            r0 = (r0 < 0) ? (ibs0finder.posvec_que.back())->back() : r0;
            l0 = (l0 < 0) ? ibs0finder.start_que[0] : l0;
            lh->ibs = r0 - l0 - pgmap.overlap_centro(l0,r0);
  // std::cout << l0 << ' ' << r0 << ' ' << lh->ibs << '\n';
            if (ibd < lh->hap) {
              besthap = lh->hapid;
              ibd = lh->hap;
            } else if (besthap >= 0 && ibd == lh->hap && lh->ibs > pt.reflist[besthap]->ibs) {
              besthap = lh->hapid;
            }
            if (ibs < lh->ibs) {
              bestibs = lh->hapid;
              ibs = lh->ibs;
            } else if (bestibs >= 0 && ibs == lh->ibs && lh->hap > pt.reflist[bestibs]->hap) {
              bestibs = lh->hapid;
            }
            pt.reflist[lh->hapid]=lh;
            if (detailed) {
              fprintf(wf, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",pos,id1,h%2,id2,lh->hapid%2, dijp,dijs,l0,r0);
            }
          }
          // Right
          if (i_s == M-1 || dmat[k][i_s] > dmat[k][i_s+1]) {
            j_s = i_s-1;
            dijs = dmat[k][i_s];
          } else {
            j_s = i_s+1;
            dijs = dmat[k][i_s+1];
          }
          id2 = amat[k][j_s]/2;
          if (pt.reflist.find(amat[k][j_s]) != pt.reflist.end() || id2 == id1) {
            continue;
          }
          j_p = prepc.FindIndex(amat[k][j_s]);
          refhap* rh = new refhap;
          rh->hapid = amat[k][j_s];
          int32_t p_i = std::min(i_p,j_p);
          int32_t p_j = std::max(i_p,j_p);
          dijp = prepc.d[p_i+1];
          for (int32_t it = p_i+1; it <= p_j; ++it) {
            dijp = std::max(dijp,prepc.d[it]);
          }
          rh->hap = dijs - dijp - pgmap.overlap_centro(dijp,dijs);
          rh->id=id2;
          r0 = ibs0finder.FindIBS0(id1,id2,pos,0);
          l0 = ibs0finder.FindIBS0(id1,id2,pos,1);
          rh->trunc = (r0 < 0 || l0 < 0);
          r0 = (r0 < 0) ? (ibs0finder.posvec_que.back())->back() : r0;
          l0 = (l0 < 0) ? ibs0finder.start_que[0] : l0;
          rh->ibs = r0 - l0 - pgmap.overlap_centro(l0,r0);
// std::cout << l0 << ' ' << r0 << ' ' << rh->ibs << '\n';
          if (ibd < rh->hap) {
            besthap = rh->hapid;
            ibd = rh->hap;
          } else if (besthap >=0 && ibd == rh->hap && rh->ibs > pt.reflist[besthap]->ibs) {
            besthap = rh->hapid;
          }
          if (ibs < rh->ibs) {
            bestibs = rh->hapid;
            ibs = rh->ibs;
          } else if (bestibs >=0 && ibs == rh->ibs && rh->hap > pt.reflist[bestibs]->hap) {
            bestibs = rh->hapid;
          }
          pt.reflist[rh->hapid]=rh;
          if (detailed) {
            fprintf(wf, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",pos,id1,h%2,id2,rh->hapid%2, dijp,dijs,l0,r0);
          }
        }
        if (besthap >= 0 && bestibs >= 0) {
          fprintf(wfs, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", pos, id1,
                  pt.reflist[besthap]->hap, pt.reflist[besthap]->ibs,pt.reflist[besthap]->trunc,
                  pt.reflist[bestibs]->hap, pt.reflist[bestibs]->ibs,pt.reflist[bestibs]->trunc);
        }
        for (auto it : pt.reflist)
          delete it.second;
        iter++;
        if (iter >= poslist.size()) {
          flag = 1;
          break;
        }
        pos = poslist[iter][0];
      }
      k++;
      if (flag)
        break;
    } // Finish processing one chunk
    // Free memory used by the suffix pbwt
    for (auto& pt : dmat)
      delete [] pt;
    for (auto& pt : amat)
      delete [] pt;
    if (flag) {
      break;
    }
    st = ed + 1;
    ed += pbwt_chunk;
  } // Finish processing all positions
  fclose(wfs);
  if (detailed)
    fclose(wf);
  return 0;
}

















struct snp_query_simple {
  int32_t pos1, pos2;
  std::vector<int32_t> carry_id;
  std::vector<int32_t> carry_ht;
  std::vector<int32_t> left;
  std::vector<int32_t> right;
  snp_query_simple(int32_t _p1, int32_t _p2) : pos1(_p1), pos2(_p2) {}
};

bool snp_compare3 (snp_query_simple* v, snp_query_simple* u) {return (v->pos1 < u->pos1);};



int32_t hapIBDpairwise(int32_t argc, char** argv) {

  std::string inVcf, inPos, out, path_pbwt;
  int32_t verbose = 10000;
  int32_t nsamples=0, M=0;
  int32_t min_ac = 10;
  int32_t pbwt_chunk = 1000000, OFFSET = 2;

  paramList pl;
  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Sites", NULL)
    LONG_STRING_PARAM("in-vcf",&inVcf, "Input VCF/BCF file")
    LONG_STRING_PARAM("pos",&inPos, "List of positions to query")
    LONG_STRING_PARAM("pbwt-path",&path_pbwt, "Gene region bed")
    LONG_INT_PARAM("min-ac",&min_ac, "Minimum sample allele count to include in haplotype matching")

    // LONG_PARAM_GROUP("Additional Options", NULL)
    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output file prefix")
    LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")

  END_LONG_PARAMS();
  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  // bcf reader
  std::vector<GenomeInterval> intervals;
  BCFOrderedReader odr(inVcf, intervals);
  nsamples = bcf_hdr_nsamples(odr.hdr);
  M = nsamples * 2;

  // Read query position file
  std::vector<snp_query_simple* > query_pos;
  // std::map<int32_t, std::vector<int32_t> > query_id;
  // std::map<int32_t, std::vector<int32_t> > query_ht;
  std::ifstream ifs;
  std::string line, preline, reg, chrom;
  std::vector<std::string> words, ids, hts;
  ifs.open(inPos, std::ifstream::in);
  int32_t pos, p1, p2, start, end, elt;
  if (!ifs.is_open()) {
    error("Query position file cannot be opened");
  }
  while(std::getline(ifs, line)) {
    words.clear();
    split(words, "\t", line);
    if (str2int32(words[1], p1) && str2int32(words[2], p2) ) {
      ids.clear(); hts.clear();
      snp_query_simple* tmp = new snp_query_simple(p1,p2);
      split(ids, ",", words[3]);
      split(hts, ",", words[4]);
      for (auto & v : ids) {
        if (v != "" &&  str2int32(v, elt)) {
          tmp->carry_id.push_back(elt);
        }
      }
      for (auto & v : hts) {
        if (v != "" &&  str2int32(v, elt)) {
          tmp->carry_ht.push_back(elt);
        }
      }
      ids.clear();
      split(ids, ",", words[5]);
      for (auto & v : ids) {
        if (v != "" &&  str2int32(v, elt)) {
          tmp->left.push_back(elt);
        }
      }
      ids.clear();
      split(ids, ",", words[6]);
      for (auto & v : ids) {
        if (v != "" &&  str2int32(v, elt)) {
          tmp->right.push_back(elt);
        }
      }
      query_pos.push_back(tmp);
    }
  }
  chrom = words[0];
  words.clear();
  line = "";
  ifs.close();
  std::sort(query_pos.begin(), query_pos.end(), snp_compare3);
  start = query_pos[0]->pos1;
  end   = query_pos.back()->pos1;

  notice("Finish reading query positions. Range: %d - %d", start, end);

  // Initialize prefix pbwt from saved check points
  int32_t st, ed, st_read, ed_read;
  int32_t a[M], d[M];
  std::string afile, dfile;
  st = start / pbwt_chunk * pbwt_chunk + 1;
  ed = st - 1 + pbwt_chunk;
  preline = "";
  while (preline == "") {
    reg = chrom + ":" + std::to_string(st) + "-" + std::to_string(ed);
    afile = path_pbwt + "_" + reg + "_prefix.amat";
    dfile = path_pbwt + "_" + reg + "_prefix.dmat";
    ifs.open(afile, std::ifstream::in);
    pos = start - 1;
    while (pos < start && std::getline(ifs, line)) {
      words.clear();
      split(words, "\t", line);
      pos = std::stoi(words[1]);
      preline = line;
    }
    ifs.close();
    if (pos < start) {
      break;
    }
    preline = "";
    st -= pbwt_chunk;
    ed -= pbwt_chunk;
  }
  words.clear();
  split(words, "\t", preline);

  for (uint32_t i =  OFFSET; i < words.size(); ++i)
    a[i- OFFSET] = std::stoi(words[i]);

  line = "";
  ifs.open(dfile, std::ifstream::in);
  pos = start - 1;
  while (pos < start && std::getline(ifs, line)) {
    words.clear();
    split(words, "\t", line);
    pos = std::stoi(words[1]);
    preline = line;
  }
  ifs.close();

  words.clear();
  split(words, "\t", preline);

  for (uint32_t i =  OFFSET; i < words.size(); ++i)
    d[i- OFFSET] = std::stoi(words[i]);

  pbwtCursor prepc(M, a, d);
  st_read = std::stoi(words[1]) + 1;

  notice("Finish finitilize pbwt prefix");
// Initialize suffix pbwt from saved check points

  st = end / pbwt_chunk * pbwt_chunk + 1;
  ed = st - 1 + pbwt_chunk;
  preline = "";
  while (preline == "") {
    reg = chrom + ":" + std::to_string(st) + "-" + std::to_string(ed);
    afile = path_pbwt + "_" + reg + "_suffix.amat";
    dfile = path_pbwt + "_" + reg + "_suffix.dmat";
    ifs.open(afile, std::ifstream::in);
    pos = end + 1;
    while (pos > end && std::getline(ifs, line)) {
      words.clear();
      split(words, "\t", line);
      pos = std::stoi(words[1]);
      preline = line;
    }
    ifs.close();
    if (pos > end) {
      break;
    }
    preline = "";
    st += pbwt_chunk;
    ed += pbwt_chunk;
  }
  words.clear();
  split(words, "\t", preline);

  for (uint32_t i =  OFFSET; i < words.size(); ++i)
    a[i- OFFSET] = std::stoi(words[i]);

  line = "";
  ifs.open(dfile, std::ifstream::in);
  pos = end + 1;
  while (pos > end && std::getline(ifs, line)) {
    words.clear();
    split(words, "\t", line);
    pos = std::stoi(words[1]);
    preline = line;
  }
  ifs.close();

  words.clear();
  split(words, "\t", preline);

  for (uint32_t i =  OFFSET; i < words.size(); ++i)
    d[i- OFFSET] = std::stoi(words[i]);

  pbwtCursor sufpc(M, a, d);
  ed_read = std::stoi(words[1]) - 1;
  notice("Finish finitilize pbwt suffix");

  reg = chrom + ":" + std::to_string(st_read) + "-" + std::to_string(ed_read);

  // reset bcf reader range
  parse_intervals(intervals, "", reg);
  odr.jump_to_interval(intervals[0]);
  bcf1_t* iv = bcf_init();

  // Output
  FILE * wf;
  wf = fopen(out.c_str(), "w");
  fprintf(wf, "%s\n", "POS\tRefPOS\tID\tHT\tIBD1\tIBD2\tDirection");


  std::vector<bool*> gtmat;
  std::vector<int32_t> positions;

  int32_t *p_gt = NULL;
  int32_t  n_gt = 0;
  int32_t *info_ac = NULL;
  int32_t n_ac = 0;

  int32_t query_iter = 0;

  for (int32_t k=0; odr.read(iv); ++k) { // Read haplotypes
    if ( k % verbose == 0 )
      notice("Processing %d variants at %s:%d.", k, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);
    if (bcf_get_genotypes(odr.hdr, iv, &p_gt, &n_gt) < 0) {
      error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);
    }
    if (bcf_get_info_int32(odr.hdr, iv, "AC", &info_ac, &n_ac) < 0) {continue;}
    // int32_t minor = 1;
    // if (info_ac[0] > nsamples) {minor = 0;}

    if (query_iter < (int32_t) query_pos.size() && iv->pos+1 == query_pos[query_iter]->pos1) {
      int32_t rvec[M];
      prepc.ReverseA(rvec);
      while (query_iter < (int32_t) query_pos.size() && iv->pos+1 == query_pos[query_iter]->pos1) {
        p1 = query_pos[query_iter]->pos1;
        p2 = query_pos[query_iter]->pos2;
        // Get carriers
        std::vector<int32_t> & c_id = query_pos[query_iter]->carry_id;
        std::vector<int32_t> & c_ht = query_pos[query_iter]->carry_ht;
        // Get haplotype matching
        if (c_id.size() > 1 && query_pos[query_iter]->left.size() > 0) {
          // query_id[pos] = c_id;
          // query_ht[pos] = c_ht;
          // for (uint32_t i = 0; i < c_id.size(); ++i) {
          for (uint32_t i = 0; i < query_pos[query_iter]->left.size(); ++i) {
            int32_t q_id = query_pos[query_iter]->left[i];
            std::stringstream d1, d2;
            for (uint32_t j = 0; j < c_id.size(); ++j) {
              if (c_id[j] != q_id ) {
                d1 << p1 - prepc.Dist_pref(rvec[2*q_id], rvec[2*c_id[j]+c_ht[j]]);
                d1 << ',';
                d2 << p1 - prepc.Dist_pref(rvec[2*q_id+1], rvec[2*c_id[j]+c_ht[j]]);
                d2 << ',';
              }
            }
            fprintf(wf, "%d\t%d\t%d\t%d\t%s\t%s\t%d\n", p1, p2, q_id, c_ht[i], d1.str().c_str(), d2.str().c_str(), 0 );
          }
        }
        query_iter += 1;
      }
      continue;
    }

    if (info_ac[0] < min_ac) {continue;}

    bool *y = new bool[M];
    int32_t ac = 0;
    for (int32_t i = 0; i <  nsamples; ++i) {
      int32_t g1 = p_gt[2*i];
      int32_t g2 = p_gt[2*i+1];
      if (bcf_gt_is_missing(g1)) {
        y[2*i] = 0;
      } else {
        y[2*i] = (bcf_gt_allele(g1) > 0);
      }
      if (bcf_gt_is_missing(g2)) {
        y[2*i+1] = 0;
      } else {
        y[2*i+1] = (bcf_gt_allele(g2) > 0);
      }
      ac += (y[2*i] + y[2*i+1]);
    } // Finish reading one SNP
    prepc.ForwardsAD_prefix(y, iv->pos+1);
    gtmat.push_back(y);
    positions.push_back(iv->pos+1);
  } // Finish reading all haplotypes in this chunk

  uint32_t N = positions.size();
  if (N != gtmat.size()) {
    error("Error in scanning forward: unmatched gtmat & positions");
  }

  notice("SNPs used in haplotype matrix: %d. Start right matching.", N);

  query_iter = (int32_t) query_pos.size() - 1;
  for (uint32_t k = N-1; k >= 0; --k) {
    if (positions[k] <= query_pos[query_iter]->pos1) {
      int32_t rvec[M];
      sufpc.ReverseA(rvec);
      while (positions[k] <= query_pos[query_iter]->pos1) {
        p1 = query_pos[query_iter]->pos1;
        p2 = query_pos[query_iter]->pos2;
        if (query_pos[query_iter]->right.size() > 0) {
          std::vector<int32_t> & c_id = query_pos[query_iter]->carry_id;
          std::vector<int32_t> & c_ht = query_pos[query_iter]->carry_ht;
          // Get haplotype matching
          for (uint32_t i = 0; i < query_pos[query_iter]->right.size(); ++i) {
            int32_t q_id = query_pos[query_iter]->right[i];
            std::stringstream d1, d2;
            for (uint32_t j = 0; j < c_id.size(); ++j) {
              if (c_id[j] != q_id ) {
                d1 << sufpc.Dist_suff(rvec[2*q_id], rvec[2*c_id[j]+c_ht[j]]) - p1;
                d1 << ',';
                d2 << sufpc.Dist_suff(rvec[2*q_id+1], rvec[2*c_id[j]+c_ht[j]]) - p1;
                d2 << ',';
              }
            }
            fprintf(wf, "%d\t%d\t%d\t%d\t%s\t%s\t%d\n", p1, p2, q_id, c_ht[i], d1.str().c_str(), d2.str().c_str(), 1 );
          }
        }
        query_iter -= 1;
        if (query_iter < 0) {break;}
      }
    }
    if (query_iter < 0) {break;}
    sufpc.ForwardsAD_suffix(gtmat[k], positions[k]);
  }

  fclose(wf);
  return 0;
}

















