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

  // Read the position - id file
  std::fstream rf(inPos,std::ios::in);
  std::map<int32_t, int32_t> poslist;
  int32_t pos, id;
  int32_t start, end;
  while (rf >> pos >> id) {
    poslist[pos] = id;
  }

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
    start = poslist.begin()->first;
    end   = poslist.rbegin()->first;
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
    if (ac > min_ac)
      prepc.ForwardsAD_prefix(y, iv->pos+1);
  } // Finish reading all haplotypes before the first position


  // Start finding ibd around input position-id pairs.
  int32_t st = start;
  int32_t ed = (start/pbwt_chunk + 1)*pbwt_chunk;

  // initialize ibs0 lookup
  reg = chrom + ":" + std::to_string(st) + "-" + std::to_string(ed);
  IBS0lookup ibs0finder(inVcf, reg, pgmap, delta, ibs0_chunk, min_hom_gts);

  auto iter = poslist.begin();
  bool flag = 0;
  pos = iter->first;

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
      if (ac <= min_ac)
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
    int32_t ibd = 0, ibs = 0, besthap, bestibs;
    int32_t id1, id2; // individual id
    int32_t i_s, j_s, i_p, j_p;
    int32_t dijp, dijs;
    int32_t k = 0;
    while (k <= N-1 && pos <= ed) {
      // Sorted upto and include position k
      prepc.ForwardsAD_prefix(gtmat[k], positions[k]);
      delete [] gtmat[k];
      if (positions[k] < pos) {
        k++;
        continue;
      }
      while (positions[k] > pos) {
        iter++;
        if (iter == poslist.end())
          break;
        pos = iter->first;
      }
      if (iter == poslist.end()) {
        break;
      }
std::cout << k << '\t' << pos << '\t' << positions[k] << std::endl;
      id1 = iter->second;
      switchpt pt;
      pt.pos=pos;
      pt.id =id1;
      int32_t rvec[M]; // Map haplotype id to row in suffix matrix
      for (int32_t it = 0; it < M; ++it) {
        rvec[amat[k][it]] = it;
      }
      for (int32_t h = id1*2; h <= id1*2+1; ++h) {
        i_p = prepc.FindIndex(h);
        i_s = rvec[i_p];

        // Left
        if (i_p == M-1 || prepc.d[i_p] > prepc.d[i_p+1]) {
          j_p = i_p-1;
          dijp = prepc.d[i_p];
        } else {
          j_p = i_p+1;
          dijp = prepc.d[i_p+1];
        }
        if (pt.reflist.find(prepc.a[j_p]) != pt.reflist.end() || prepc.a[j_p]/2 == id1) {
          continue;
        }
        id2 = prepc.a[j_p]/2;
        j_s = rvec[prepc.a[j_p]];
        refhap* lh = new refhap;
        lh->hapid = prepc.a[j_p];
        int32_t s_i = std::min(i_s,j_s);
        int32_t s_j = std::max(i_s,j_s);
        dijs = dmat[k][s_i+1];
        for (int32_t it = s_i+1; it <= s_j; ++it) {
          dijs = std::min(dijs,dmat[k][it]);
        }
        lh->hap = dijs - dijp;
        lh->id=id2;
        int32_t r0 = ibs0finder.FindIBS0(id1,id2,pos,0);
        int32_t l0 = ibs0finder.FindIBS0(id1,id2,pos,1);
        lh->trunc = (r0 || l0);
        r0 = (r0 < 0) ? (ibs0finder.posvec_que.back())->back() : r0;
        l0 = (l0 < 0) ? ibs0finder.start_que[0] : l0;
        lh->ibs = r0 - l0;
        if (ibd < lh->hap) {
          besthap = lh->hapid;
          ibd = lh->hap;
        }
        if (ibs < lh->ibs) {
          bestibs = lh->hapid;
          ibs = lh->ibs;
        }
        pt.reflist[lh->hapid]=lh;
        if (detailed) {
          fprintf(wf, "%d\t%d\t%d\t%d\t%d\t%d\t%d\n",id1,h/2,prepc.a[j_p], dijp,dijs,l0,r0);
        }

        // Right
        if (i_s == M-1 || dmat[k][i_s] > dmat[k][i_s+1]) {
          j_s = i_s-1;
          dijs = dmat[k][i_s];
        } else {
          j_s = i_s+1;
          dijs = dmat[k][i_s+1];
        }
        if (pt.reflist.find(amat[k][j_s]) != pt.reflist.end()) {
          continue;
        }
        id2 = amat[k][j_s]/2;
        j_p = prepc.FindIndex(amat[k][j_s]);
        refhap* rh = new refhap;
        rh->hapid = amat[k][j_s];
        int32_t p_i = std::min(i_p,j_p);
        int32_t p_j = std::max(i_p,j_p);
        dijp = prepc.d[s_i+1];
        for (int32_t it = p_i+1; it <= p_j; ++it) {
          dijp = std::max(dijp,prepc.d[it]);
        }
        rh->hap = dijs - dijp;
        rh->id=id2;
        r0 = ibs0finder.FindIBS0(id1,id2,pos,0);
        l0 = ibs0finder.FindIBS0(id1,id2,pos,1);
        rh->trunc = (r0 || l0);
        r0 = (r0 < 0) ? (ibs0finder.posvec_que.back())->back() : r0;
        l0 = (l0 < 0) ? ibs0finder.start_que[0] : l0;
        rh->ibs = r0 - l0;
        if (ibd < rh->hap) {
          besthap = rh->hapid;
          ibd = rh->hap;
        }
        if (ibs < rh->ibs) {
          bestibs = rh->hapid;
          ibs = rh->ibs;
        }
        pt.reflist[rh->hapid]=rh;
        if (detailed) {
          fprintf(wf, "%d\t%d\t%d\t%d\t%d\t%d\t%d\n",id1,h/2,prepc.a[j_p], dijp,dijs,l0,r0);
        }
      }
      fprintf(wfs, "%d\t%d\t%d\t%d\t%d\t%d\n", pos, id1,
              pt.reflist[besthap]->hap, pt.reflist[besthap]->ibs,
              pt.reflist[bestibs]->hap, pt.reflist[bestibs]->ibs);
      for (auto it : pt.reflist)
        delete it.second;
      k++;
      iter++;
      if (iter == poslist.end()) {
        flag = 1;
        break;
      }
      pos = iter->first;
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




