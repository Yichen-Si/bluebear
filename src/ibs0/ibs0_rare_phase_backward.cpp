#include "cramore.h"
#include "bcf_filter_arg.h"
#include "bcf_ordered_reader.h"
#include "compact_matrix.h"

#include "ibs0.h"
#include "pbwt_build.h"
#include "cm2bp.h"
#include "bp2cm.h"

#include <functional>
#include <iomanip>

class Candidate {
public:
  int32_t id;
  std::vector<int32_t> supporter;
  std::vector<int32_t> extension;
  bool flag; // Priority / Must in.

  Candidate(int32_t id1, int32_t id2, int32_t l, int32_t type) : id(id1) {
    supporter.push_back(id2);
    extension.push_back(l);
    flag = flag || type;
  }
};

inline bool MoreSupporters(Candidate* x, Candidate* y) {
  return (x->supporter.size() > y->supporter.size());
}

inline bool LongerExtension(Candidate* x, Candidate* y) {
  if (x->flag == y->flag) {
    std::sort(x->extension.begin(), x->extension.end());
    std::sort(y->extension.begin(), y->extension.end());
    return (x->extension.back() > y->extension.back());
  } else {
    return x->flag;
  }
}

int32_t RareIBS0PhaseBackward(int32_t argc, char** argv) {
  std::string inVcf, inMap, path_pbwt, out, reg;
  std::string outf,outf_s,outVcf;
  int32_t detailed = 0;
  std::string chrom="chr20";
  std::string mode ="b";
  int32_t verbose = 1000;
  int32_t max_flip = 3;
  int32_t min_variant = 2;
  int32_t min_hom_gts = 1;
  int32_t nsamples=0, M=0;
  int32_t st, ck, ed, edugly;
  int32_t chunk_size = 1000000;          // process suffix pbwt & ibs0 by chunk
  int32_t ibs0_chunk = 500000;
  int32_t start_pos = 0; // Starting position

  double  delta = 3.0;                   // thresholds
  double  mafcut = 0.05;
  int32_t rare_ac = 2;
  int32_t lambda = 20000, gamma = 20000, diff = 5000, nmatch = 50000;

  int32_t* p_gt = NULL;
  int32_t  n_gt = 0;

  inMap = "/net/wonderland/home/ycsi/Data/Map/plink.chr20.GRCh38.map";

  paramList pl;
  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input", NULL)
    LONG_STRING_PARAM("in-vcf",&inVcf, "Input VCF/BCF file")
    LONG_STRING_PARAM("region",&reg,"Genomic region to focus on")
    LONG_STRING_PARAM("chr",&chrom,"Chromosome to process")
    LONG_STRING_PARAM("map",&inMap, "Map file for genetic distance")
    LONG_STRING_PARAM("pbwt-path",&path_pbwt,"Pre-computed snapshot of suffix pbwt")

    LONG_PARAM_GROUP("Additional Options", NULL)
    LONG_DOUBLE_PARAM("delta",&delta, "no-ibs0 threshold (in cM) to be a proxy of IBD")
    LONG_DOUBLE_PARAM("mafcut",&mafcut, "Sample minor allele frequency for considering flipping")
    LONG_INT_PARAM("rare-ac",&rare_ac,"Rare allele sharing as evidence for IBD")
    LONG_INT_PARAM("lambda",&lambda,"Length (in bp) to the end of no-ibs0 that is not trusted")
    LONG_INT_PARAM("gamma",&gamma,"Length (in bp) of matched prefix that is requred to begin searching")
    LONG_INT_PARAM("min-diff",&diff,"The minimal difference of new matches (in bp) between flipping two individuals for a flip to be recorded")
    LONG_INT_PARAM("min-improve",&nmatch,"The minimal improve of matching (in bp) for a flip to be recorded")
    LONG_INT_PARAM("max-flip",&max_flip,"Maximal #flips per site")
    LONG_INT_PARAM("start-from",&start_pos,"Starting position. Proceed on both sides.")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output file prefix")
    LONG_STRING_PARAM("mode", &mode, "Output format of flipped bcf/vcf (b/bu/z)")
    LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")
    LONG_INT_PARAM("detailed",&detailed,"If detailed info for each pair is to be written (0/1)")

  END_LONG_PARAMS();
  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  std::ostringstream parm;
  parm << rare_ac << '-' << max_flip << '-' << (int32_t) (mafcut*100) << '-' << std::fixed << std::setprecision(1) << delta << '-' << lambda/1000 << '-' << gamma/1000 << '-' << diff/1000 << '-' << nmatch/1000;
  notice("Parameter setting: %s.", parm.str().c_str());

  // Translate between genetic & physical position.
  cm2bpMap gpmap(inMap, " ");
  bp2cmMap pgmap(inMap, " ");
  notice("Will proceed until the first position in linkage map: %d.", gpmap.minpos);

  // Output flipped bcf/vcf
  outVcf = out + "_" + parm.str() + "_st_" + std::to_string(start_pos) + "_bw.bcf";
  if (mode=="z")
    outVcf = out + "_" + parm.str() + "_st_" + std::to_string(start_pos) + "_bw.vcf.gz";
  mode = "w" + mode;
  htsFile *wbcf = hts_open(outVcf.c_str(),mode.c_str());

  // Copy header & get sample size
  htsFile *fp = hts_open(inVcf.c_str(),"r");
  bcf_hdr_t *hdr = bcf_hdr_read(fp);
  bcf_hdr_write(wbcf, hdr);
  nsamples = hdr->n[BCF_DT_SAMPLE];
  M = nsamples * 2;
  hts_close(fp);

  // If the current position is flipped
  bool flip_abs[nsamples] = {0};

  ck = start_pos/chunk_size; // Chunk in terms of chunk_size
  st = ck * chunk_size + 1;
  ed = st + chunk_size - 1;
  while (st > pgmap.centromere_st && ed < pgmap.centromere_ed ) {
    ck--;
    st = ck * chunk_size + 1;
    ed = st + chunk_size - 1;
  }
  edugly = ed;

  // Initialize ibs0 lookup blocks
  reg = chrom + ":" + std::to_string(st) + "-" + std::to_string(ed);
  RareIBS0lookup ibs0finder(inVcf, reg, pgmap, delta, rare_ac, ibs0_chunk, min_hom_gts);
  notice("Initilized ibs0 lookup (%.1f cM, %d blocks).", delta, ibs0finder.start_que.size());

  // Initialize backward pbwt
  // Read the last snp in this chunk. should be stored as a checkpoint
  reg = chrom + ":" + std::to_string(st) + "-" + std::to_string(ed);
  std::string line;
  std::string sfile = path_pbwt + "_" + reg + "_suffix.amat";

  std::ifstream sf(sfile);
  std::getline(sf, line);
  sf.close();
  std::vector<std::string> wvec;
  split(wvec, "\t", line);
  uint32_t OFFSET = wvec.size() - M;
  int32_t a[M];
  for (uint32_t i =  OFFSET; i < wvec.size(); ++i)
    a[i- OFFSET] = std::stoi(wvec[i]);

  sfile = path_pbwt + "_" + reg + "_suffix.dmat";
  sf.open(sfile);
  std::getline(sf, line);
  sf.close();
  split(wvec, "\t", line);
  int32_t d[M];
  for (uint32_t i =  OFFSET; i < wvec.size(); ++i)
    d[i- OFFSET] = std::stoi(wvec[i]);

  pbwtCursor sufpc(M, a, d);

  while (edugly > gpmap.minpos) {
    // Read and build prefix pbwt by physical chunk
    // TODO: enable efficient random access &/ customizable chunk size

    // Read genotype matrix
    reg = chrom + ":" + std::to_string(st) + "-" + std::to_string(edugly);
    ibs0finder.Update(reg);
    std::vector<GenomeInterval> intervals;
    parse_intervals(intervals, "", reg);
    BCFOrderedReader odr(inVcf, intervals);
    bcf1_t* iv = bcf_init();

    std::vector<bool*> gtmat;
    std::vector<int32_t> positions;
    std::vector<double> mafs;

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
      gtmat.push_back(y);
      positions.push_back(iv->pos+1);
      mafs.push_back(ac*1.0/M);
    } // Finish reading all haplotypes in this chunk
    int32_t N = positions.size();
    if ( N < min_variant ) {
      if (N > 0) {
        for (int32_t k = N-1; k >= 0; --k)
          sufpc.ForwardsAD_suffix(gtmat[k], positions[k]);
      }
      ck--;
      st = ck * chunk_size + 1;
      ed = st + chunk_size - 1;
      edugly = ed;
      for (auto& pt : gtmat)
        delete [] pt;
      notice("Observed only %d informative markers. Skipping the region %s", N, reg.c_str());
      continue;
    } else {
      notice("Read %d markers, start processing %s.", N, reg.c_str());
    }

    // Set up output for this block
    reg = chrom + ":" + std::to_string(st) + "-" + std::to_string(ed);
    std::ofstream wf;
    if (detailed) {
      outf = out + "_" + parm.str() + "_st_" + std::to_string(start_pos) + "_bw_" + reg + ".list";
      wf.open(outf, std::ios::trunc);
    }
    outf_s = out + "_"+ parm.str() + "_st_" + std::to_string(start_pos) + "_bw_" +reg+"_short.list";
    std::ofstream wfs(outf_s, std::ios::trunc);

    // Read the first snp in this chunk. should be stored as a checkpoint
    sfile = path_pbwt + "_" + reg + "_prefix.amat";
    sf.open(sfile);
    std::getline(sf, line);
    sf.close();
    split(wvec, "\t", line);
    OFFSET = wvec.size() - M;
    for (uint32_t i =  OFFSET; i < wvec.size(); ++i)
      a[i- OFFSET] = std::stoi(wvec[i]);

    sfile = path_pbwt + "_" + reg + "_prefix.dmat";
    sf.open(sfile);
    std::getline(sf, line);
    sf.close();
    split(wvec, "\t", line);
    for (uint32_t i =  OFFSET; i < wvec.size(); ++i)
      d[i- OFFSET] = std::stoi(wvec[i]);

    pbwtCursor prepc(M, a, d);

    // Allocate space for pbwt matricies
    std::vector<int32_t*> dmat(N,NULL);
    std::vector<int32_t*> rmat(N,NULL);
    for (int32_t i = 0; i < N; ++i) {
      dmat[i] = new int32_t[M];
      rmat[i] = new int32_t[M];
    }
    // Build prefix pbwt
    for (int32_t k = 0; k < N; ++k) {
      prepc.ForwardsAD_prefix(gtmat[k], positions[k]);
      memcpy(dmat[k], prepc.d, M*sizeof(int32_t));
      prepc.ReverseA(rmat[k]);
    }
    // Build suffix pbwt & detect switch
    int32_t h11,h12,h21,h22; // (Reflecting flips so far) index stored in suffix pbwt
    int32_t i_s, j_s, i_s_prime, j_s_prime; // row num in suffix pbwt
    int32_t dij_p, dipj_s, dijp_s; // absolute position, from pbwt divergence matrix
    for (int32_t k = N-1; k > 0; --k) {
      // Sorted after and include position k
      sufpc.ForwardsAD_suffix(gtmat[k], positions[k]);
      // gtmat[k] is no longer needed. Update to reflect flipped stat
      for (int32_t it = 0; it < nsamples; ++it) {
        if (flip_abs[it] && (gtmat[k][it*2] != gtmat[k][it*2+1]) ) {
          gtmat[k][it*2] = !gtmat[k][it*2];
          gtmat[k][it*2+1] = !gtmat[k][it*2+1];
        }
      }
      if (positions[k]>start_pos) {
        bcf_write(wbcf, odr.hdr, iv);
        continue;
        } // Start_pos & after should not be flipped
      std::vector<std::vector<int32_t> > fliprec;
      std::map<int32_t, Candidate *> flipcandy;
      for (int32_t i = 0; i < M-1; ++i) {
        h11 = sufpc.a[i]; // Original Hap ID corresponding to row i
        h12 = h11 + 1 - 2 * (h11%2);
        int32_t j = i+1;
        dij_p = sufpc.d[j];
        // Check a few hap down
        while (j < M && dij_p > positions[k] + gamma) {
          h21 = sufpc.a[j];
          h22 = h21 + 1 - 2 * (h21%2);
          if (h11/2 == h21/2) { // Probably long homozygous sites
            j++;
            if (sufpc.d[j] < dij_p)
              dij_p = sufpc.d[j];
            continue;
          }
          if (gtmat[k-1][h11] != gtmat[k-1][h21]) {
          // If a match ends at k-1
            i_s = rmat[k-1][h11]; j_s = rmat[k-1][h21];
            i_s_prime = rmat[k-1][h12]; j_s_prime = rmat[k-1][h22];
            // If flip individual 1
            int32_t lower = std::min(i_s_prime, j_s);
            int32_t upper = std::max(i_s_prime, j_s);
            dipj_s = dmat[k-1][lower + 1];
            for (int32_t it = lower + 1; it <= upper; ++it) { // Evaluate d_i'j+
              if (dmat[k-1][it] > dipj_s)
                dipj_s = dmat[k-1][it];
            }
            // If flip individual 2
            lower = std::min(i_s, j_s_prime);
            upper = std::max(i_s, j_s_prime);
            dijp_s = dmat[k-1][lower + 1];
            for (int32_t it = lower + 1; it <= upper; ++it) { // Evaluate d_ij'+
              if (dmat[k-1][it] > dijp_s)
                dijp_s = dmat[k-1][it];
            }
            int32_t thewho, theother, thelong, theshort;
            if (dijp_s < dipj_s) { // Flip individual 2
              thelong  = dijp_s; theshort = dipj_s;
              thewho   = h21/2; theother = h11/2;
            } else {              // Flip individual 1
              theshort = dijp_s; thelong  = dipj_s;
              theother = h21/2; thewho   = h11/2;
            }
            if (theshort - thelong > diff && positions[k] - thelong > nmatch) {
              int32_t flag = 0;
              int32_t nextrare = ibs0finder.FindRare(h11/2,h21/2,positions[k],0);
              int32_t prevrare = ibs0finder.FindRare(h11/2,h21/2,positions[k],1);
              int32_t previbs0 = ibs0finder.FindIBS0(h11/2,h21/2,positions[k],1);
              int32_t nextibs0 = 0;
              if (nextrare > 0 || prevrare > 0) { // Conditional on rare allele sharing
                if (prevrare > 0) {
                  if (previbs0 < 0 || previbs0 < prevrare) {
                    flag = 2;
                  } else if (mafs[k-1] < mafcut &&
                             pgmap.bp2cm(dij_p) - pgmap.bp2cm(previbs0) > delta) {
                    flag = 3;
                  }
                }
                if (flag != 2 && nextrare > 0) {
                  nextibs0 = ibs0finder.FindIBS0(h11/2,h21/2,positions[k],0);
                  if  (nextibs0 < 0 || (nextibs0 > nextrare)) {
                    flag = 2;
                  } else if (mafs[k-1] < mafcut &&
                             pgmap.bp2cm(nextibs0) - pgmap.bp2cm(previbs0) > delta) {
                    flag = 3;
                  }
                }
              }
              if (!flag && mafs[k-1] < mafcut) { // Conditional on maf & ibs0
                if (previbs0 < 0 || (pgmap.bp2cm(dij_p) - pgmap.bp2cm(previbs0) > delta)) {flag= 3;}
                if (pgmap.bp2cm(dij_p)-pgmap.bp2cm(thelong) > delta) {flag = 1;}
                if (!flag) {
                  if (nextrare == 0)
                    nextibs0 = ibs0finder.FindIBS0(h11/2,h21/2,positions[k],0);
                  if (pgmap.bp2cm(nextibs0) - pgmap.bp2cm(previbs0) > delta) {flag = 3;}
                }
              }
              if (flag) { // Flip
                fliprec.push_back(std::vector<int32_t> {positions[k+1],thewho,theother,
                                                        dij_p, theshort, thelong, flag});
                if (flipcandy.find(thewho) != flipcandy.end()) {
                  flipcandy[thewho]->supporter.push_back(theother);
                  flipcandy[thewho]->extension.push_back(positions[k]-thelong);
                } else {
                  Candidate* candy = new Candidate(thewho,theother,positions[k]-thelong,(flag == 2));
                  flipcandy[thewho] = candy;
                }
              }
            }
          }
          j++;
          if (sufpc.d[j] < dij_p)
            dij_p = sufpc.d[j];
        } // End of searching switch of pairs involving row i
      } // End of searching for this column
      // If multiple flips occur
      if (flipcandy.size() > 0) {
        int32_t adj_max_flip = max_flip;
        std::map<int32_t, std::vector<int32_t> > finalrec;
        std::vector<bool> flipped(nsamples, 0);
        if (flipcandy.size() > 1) {
          std::vector<Candidate*> toflip;
          for (auto & v : flipcandy) {
            if (v.second->flag) {adj_max_flip++;}
            toflip.push_back(v.second);
          }
          std::sort(toflip.begin(), toflip.end(), LongerExtension);
          for (auto & v : toflip) {
            if ((int32_t) finalrec.size() >= adj_max_flip) {
              break;
            }
            int32_t rm = (v->supporter).size();
            for (auto & id2 : v->supporter) {
              if (flipped[id2])
                rm--;
            }
            if (rm) {
              flipped[v->id] = 1;
              finalrec[v->id] = std::vector<int32_t>{0,v->id,0,0, 0,0,0};
            }
          }
        } else {
          flipped[flipcandy.begin()->first] = 1;
          finalrec[flipcandy.begin()->first] = std::vector<int32_t>{0,flipcandy.begin()->first,0,0, 0,0,0};
        }
        for (auto & v : fliprec) {
          if (flipped[v[1]]) {
            finalrec[v[1]][0] = v[0];
            finalrec[v[1]][2]++;      // Number of comparisons suggesting this flip
            if (v[4]-v[5] > finalrec[v[1]][3])
              finalrec[v[1]][3] = v[4]-v[5];
            finalrec[v[1]][3+v[6]]++; // Count the number of each types of flip among those comparisons
            if (detailed) {
              std::stringstream recline;
              for (auto& w : v)
                recline << w << '\t';
              recline.seekp(-1, std::ios_base::end);
              recline << '\n';
              wf << recline.str();
            }
          }
        }
        for (auto & v : finalrec) {
          std::stringstream recline;
          for (auto& w : v.second)
            recline << w << '\t';
          recline.seekp(-1, std::ios_base::end);
          recline << '\n';
          wfs << recline.str();
        }
        sufpc.SwitchIndex(flipped);
        for (int32_t it = 0; it < nsamples; ++it) {
          if (flipped[it]) {
            flip_abs[it] = !flip_abs[it];
          }
        }
      }
    } // End of processing one block

    wfs.close();
    // Output this block. Up till and NOT including start_pos
    if (edugly > start_pos-2) {
      // If it is the first block, skip sites after the starting position.
      reg = chrom + ":" + std::to_string(st) + "-" + std::to_string(start_pos-2);
      parse_intervals(intervals, "", reg);
    }
    odr.jump_to_interval(intervals[0]);
    bcf_clear(iv);
    int32_t y[M];
    odr.read(iv);
    int32_t k = 1;
    while(odr.read(iv)) {
      for (int32_t it = 0; it < nsamples; ++it) {
        y[it*2] = bcf_gt_phased(gtmat[k][it*2]);
        y[it*2+1] = bcf_gt_phased(gtmat[k][it*2+1]);
      }
      if(bcf_update_genotypes(odr.hdr, iv, y, M)) {
        error("Cannot update GT.");
      }
      if (bcf_write(wbcf, odr.hdr, iv))
        error("Error in writing bcf.");
      k++;
    }
    ck--;
    st = ck * chunk_size + 1;
    edugly = positions[0];
    ed = st + chunk_size - 1;
    // Free memory of gt matrix
    for (auto& pt : gtmat)
      delete [] pt;
    // Free memory used by the suffix pbwt
    for (auto& pt : dmat)
      delete [] pt;
    for (auto& pt : rmat)
      delete [] pt;
  } // Finish processing one chunk
  hts_close(wbcf);
  return 0;
}




