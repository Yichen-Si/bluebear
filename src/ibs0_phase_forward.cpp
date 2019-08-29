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

int32_t IBS0PhaseForward(int32_t argc, char** argv) {
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
  int32_t st, ck, ed, ibs_st, ibs_ck, ibs_ed, stugly;
  int32_t chunk_size = 1000000;          // process suffix pbwt & ibs0 by chunk
  int32_t start_pos = 0; // Starting position

  double  delta = 3.0;                   // thresholds
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
  parm << std::fixed << std::setprecision(1) << delta << '-' << lambda/1000 << '-' << gamma/1000 << '-' << diff/1000 << '-' << nmatch/1000;
  notice("Parameter setting: %s.", parm.str().c_str());

  // Translate between genetic & physical position.
  cm2bpMap gpmap(inMap, " ");
  bp2cmMap pgmap(inMap, " ");
  notice("Will proceed until the last position in linkage map: %d.", gpmap.maxpos);

  // For recording flip
  // Will first sort flip candidates by #involved pairs
  typedef std::function<bool(std::pair<int32_t, std::vector<int32_t> >, std::pair<int32_t, std::vector<int32_t> >)> Comparator;
  Comparator CompFn =
    [](std::pair<int32_t, std::vector<int32_t> > elem1 ,std::pair<int32_t, std::vector<int32_t> > elem2) {
      if ((elem1.second).size() == (elem2.second).size())
        return elem1.first > elem2.first;
      else
        return (elem1.second).size() > (elem2.second).size();};

  // Output flipped bcf/vcf
  outVcf = out + "_" + parm.str() + "_flipped_startfrom_" + std::to_string(start_pos) + "_forward.bcf";
  if (mode=="z")
    outVcf = out + "_" + parm.str() + "_flipped_startfrom_" + std::to_string(start_pos) + "_forward.vcf.gz";
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
    ck++;
    st = ck * chunk_size + 1;
    ed = st + chunk_size - 1;
  }
  stugly = st;

  // Initialize ibs0 lookup blocks. By 1Mb.
  // Number of blocks stored determined by genetic distance
  std::vector<bitmatrix*> bmatRR_que, bmatAA_que;
  std::vector<std::vector<int32_t>* > posvec_que;
  std::vector<int32_t> ibs_chunk_in_que;
  int32_t cur_ibs_ck = 0; // Currently processed position locates in which block
  ibs_ck = start_pos/chunk_size;
  ibs_st = ibs_ck * chunk_size;
  ibs_ed = (ibs_ck+1) * chunk_size-1;
  // Find the starting point for storing ibs0 look up
  while (ibs_ck > 0 && pgmap.bp2cm(start_pos) - pgmap.bp2cm(ibs_st) < delta) {
    ibs_ck--;
    ibs_ed = ibs_st - 1;
    ibs_st = std::max(ibs_ck * chunk_size, gpmap.minpos);
  }
  // Forward
  while (ibs_st < gpmap.maxpos &&  pgmap.bp2cm(ibs_st) - pgmap.bp2cm(start_pos) < delta) {
    ibs_ed = std::min(ibs_ed, gpmap.maxpos);
    reg = chrom + ":" + std::to_string(ibs_st) + "-" + std::to_string(ibs_ed);
    if (ReadIBS0Block(reg, inVcf, bmatRR_que, bmatAA_que, posvec_que, min_hom_gts, 0) > 0) {
      ibs_chunk_in_que.push_back(ibs_ck);
    }
    if (ibs_ck == start_pos/chunk_size) {
      cur_ibs_ck = (int32_t) ibs_chunk_in_que.size() - 1;
    }
    ibs_ck++;
    ibs_st = ibs_ck * chunk_size;
    ibs_ed = ibs_st + chunk_size - 1;
  } // Finish initialize ibs0 lookup blocks
  // Initialize forward pbwt
  // Read the first snp in this chunk. should be stored as a checkpoint
  reg = chrom + ":" + std::to_string(st) + "-" + std::to_string(ed);
  std::string line;
  std::string sfile = path_pbwt + "_" + reg + "_prefix.amat";

// std::cout << "initialize prefix pbwt amat " << sfile << '\n';

  std::ifstream sf(sfile);
  std::getline(sf, line);
  sf.close();
  std::vector<std::string> wvec;
  split(wvec, "\t", line);
  uint32_t OFFSET = wvec.size() - M;
  int32_t a[M];
  for (uint32_t i =  OFFSET; i < wvec.size(); ++i)
    a[i- OFFSET] = std::stoi(wvec[i]);

  sfile = path_pbwt + "_" + reg + "_prefix.dmat";

// std::cout << "initialize prefix pbwt dmat " << sfile << '\n';

  sf.open(sfile);
  std::getline(sf, line);
  sf.close();
  split(wvec, "\t", line);
  int32_t d[M];
  for (uint32_t i =  OFFSET; i < wvec.size(); ++i)
    d[i- OFFSET] = std::stoi(wvec[i]);

  pbwtCursor prepc(M, a, d);

  while (st < gpmap.maxpos) {
    // Read and build suffix pbwt by physical chunk
    // TODO: enable efficient random access &/ customizable chunk size

    // Read genotype matrix
    reg = chrom + ":" + std::to_string(stugly) + "-" + std::to_string(ed);
    std::vector<GenomeInterval> intervals;
    parse_intervals(intervals, "", reg);
    BCFOrderedReader odr(inVcf, intervals);
    bcf1_t* iv = bcf_init();

    std::vector<bool*> gtmat;
    std::vector<int32_t> positions;

    for (int32_t k=0; odr.read(iv); ++k) { // Read haplotypes
      if (bcf_get_genotypes(odr.hdr, iv, &p_gt, &n_gt) < 0) {
        error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);
      }
      bool *y = new bool[M];
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
      } // Finish reading one SNP

      if (iv->pos < start_pos) { // Start_pos & before is not flipped
        prepc.ForwardsAD_prefix(y, iv->pos+1);
        continue;
      }

      gtmat.push_back(y);
      positions.push_back(iv->pos+1);
    } // Finish reading all haplotypes in this chunk
    int32_t N = positions.size();
    if ( N < min_variant ) {
      if (N > 0) {
        for (int32_t k = 0; k < N; ++ k)
          prepc.ForwardsAD_prefix(gtmat[k], positions[k]);
      }
      ck++;
      st = ck * chunk_size + 1;
      stugly = st;
      ed = st + chunk_size - 1;
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
      outf = out + "_flip_" + parm.str() + "_startfrom_" + std::to_string(start_pos) + "_forward_" + reg + ".list";
      wf.open(outf, std::ios::trunc);
    }
    outf_s = out + "_flip_"+ parm.str() + "_startfrom_" + std::to_string(start_pos) + "_forward_" +reg+"_short.list";
    std::ofstream wfs(outf_s, std::ios::trunc);

    // Read the last snp in this chunk. should be stored as a checkpoint
    sfile = path_pbwt + "_" + reg + "_suffix.amat";
// std::cout << "Read suffix amat: " << sfile << '\n';
    sf.open(sfile);
    std::getline(sf, line);
    sf.close();
    split(wvec, "\t", line);
    OFFSET = wvec.size() - M;
    for (uint32_t i =  OFFSET; i < wvec.size(); ++i)
      a[i- OFFSET] = std::stoi(wvec[i]);

    sfile = path_pbwt + "_" + reg + "_suffix.dmat";
// std::cout << "Read suffix dmat: " << sfile << '\n';

    sf.open(sfile);
    std::getline(sf, line);
    sf.close();
    split(wvec, "\t", line);
    for (uint32_t i =  OFFSET; i < wvec.size(); ++i)
      d[i- OFFSET] = std::stoi(wvec[i]);

    pbwtCursor sufpc(M, a, d);

    // Allocate space for pbwt matricies
    std::vector<int32_t*> dmat(N,NULL);
    std::vector<int32_t*> rmat(N,NULL);
    for (int32_t i = 0; i < N; ++i) {
      dmat[i] = new int32_t[M];
      rmat[i] = new int32_t[M];
    }
    // Build suffix pbwt
    memcpy(dmat[N-1], sufpc.d, M*sizeof(int32_t));
    sufpc.ReverseA(rmat[N-1]);
    for (int32_t k = N-2; k >= 0; --k) {
      sufpc.ForwardsAD_suffix(gtmat[k], positions[k]);
      memcpy(dmat[k], sufpc.d, M*sizeof(int32_t));
      sufpc.ReverseA(rmat[k]);
    }
// std::cout << "Built suffix matrix\n";
    // Build prefix pbwt & detect switch
    int32_t h11,h12,h21,h22; // (Reflecting lips so far) index stored in prefix pbwt
    int32_t i_s, j_s, i_s_prime, j_s_prime; // row num in suffix pbwt
    int32_t dij_p, dipj_s, dijp_s; // absolute position, from pbwt divergence matrix
    int32_t bitpos = 0;
    if (stugly < start_pos) {
    // If it is the first block, skip sites before the starting position.
      reg = chrom + ":" + std::to_string(start_pos) + "-" + std::to_string(ed);
      parse_intervals(intervals, "", reg);
    }
    odr.jump_to_interval(intervals[0]);
    bcf_clear(iv);
    for (int32_t k = 0; k < N-1; ++k) {
      // Output current position
      odr.read(iv);
      if (bcf_get_genotypes(odr.hdr, iv, &p_gt, &n_gt) < 0) {
        error("Cannot find the field GT.");
      }
      int32_t y[M];
      for (int32_t it = 0; it < nsamples; ++it) {
        y[it*2] = (flip_abs[it]) ? bcf_gt_phased(gtmat[k][it*2+1]) : bcf_gt_phased(gtmat[k][it*2]);
        y[it*2+1] = (flip_abs[it]) ? bcf_gt_phased(gtmat[k][it*2]) : bcf_gt_phased(gtmat[k][it*2+1]);
      }
      if(bcf_update_genotypes(odr.hdr, iv, y, M)) {
        error("Cannot update GT.");
      }
      bcf_write(wbcf, odr.hdr, iv);

      while ((*posvec_que[cur_ibs_ck])[bitpos] < positions[k]) {bitpos++;}
      // Sorted upto and include position k
      prepc.ForwardsAD_prefix(gtmat[k], positions[k]);
      std::vector<std::vector<int32_t> > fliprec;
      std::map<int32_t, std::vector<int32_t> > flipcandy;
      for (int32_t i = 0; i < M-1; ++i) {
        h11 = prepc.a[i]; // Original Hap ID corresponding to row i
        h12 = h11 + 1 - 2 * (h11%2);
        // hap11 = (flip[h11/2]) ? (h11 + (h11%2 == 0) ? 1 : -1) : h11;
        int32_t j = i+1;
        dij_p = prepc.d[j];
        // Check a few hap down
        while (j < M && dij_p < positions[k] - gamma) {
          h21 = prepc.a[j];
          h22 = h21 + 1 - 2 * (h21%2);
          if (h11/2 == h21/2) { // Probably long homozygous sites
            j++;
            if (prepc.d[j] > dij_p)
              dij_p = prepc.d[j];
            continue;
          }
          if (gtmat[k+1][h11] != gtmat[k+1][h21]) { // If a match ends at k+1
            i_s = rmat[k+1][h11]; j_s = rmat[k+1][h21];
            i_s_prime = rmat[k+1][h12]; j_s_prime = rmat[k+1][h22];
            // If flip individual 1
            int32_t lower = std::min(i_s_prime, j_s);
            int32_t upper = std::max(i_s_prime, j_s);
            dipj_s = dmat[k+1][lower + 1];
            for (int32_t it = lower + 1; it <= upper; ++it) { // Evaluate d_i'j+
              if (dmat[k+1][it] < dipj_s)
                dipj_s = dmat[k+1][it];
            }
            // If flip individual 2
            lower = std::min(i_s, j_s_prime);
            upper = std::max(i_s, j_s_prime);
            dijp_s = dmat[k+1][lower + 1];
            for (int32_t it = lower + 1; it <= upper; ++it) { // Evaluate d_ij'+
              if (dmat[k+1][it] < dijp_s)
                dijp_s = dmat[k+1][it];
            }
            int32_t flag = 0;
            if (dijp_s - dipj_s > diff && dijp_s - positions[k] > nmatch) { // Consider flip individual 2
              if (pgmap.bp2cm(dijp_s) - pgmap.bp2cm(dij_p) > delta) {
                flag = 1;
              } else {        // Need to evaluate no-ibs0
                int32_t bpos = bitpos / 8;  // This is not exact
                int32_t ibs_ck_to_look = cur_ibs_ck;
                int32_t nextibs0 = IBS0inOneBlock(bmatRR_que[ibs_ck_to_look],
                                                  bmatAA_que[ibs_ck_to_look],
                                                  h11/2, h21/2,0,bpos);
                while (nextibs0 == -1 &&
                       ibs_ck_to_look < ((int32_t) bmatRR_que.size())-1) {
                  ibs_ck_to_look++;
                  nextibs0 = IBS0inOneBlock(bmatRR_que[ibs_ck_to_look],
                                            bmatAA_que[ibs_ck_to_look],
                                            h11/2, h21/2,0);
                }
                nextibs0 = (nextibs0 > 0) ? (*posvec_que[ibs_ck_to_look])[nextibs0] : posvec_que[ibs_ck_to_look]->back();
                if (nextibs0 - positions[k] > lambda) {
                  // Not too close to ibs0
                  if (pgmap.bp2cm(nextibs0) - pgmap.bp2cm(dij_p) > delta) {
// std::cout << dij_p << '\t' << nextibs0 << '\t' << pgmap.bp2cm(dij_p) << '\t' << pgmap.bp2cm(nextibs0) << '\n';
                    flag = 2;
                  } else { // Need to check the previous ibs0
                    ibs_ck_to_look = cur_ibs_ck;
                    int32_t previbs0 = IBS0inOneBlock(bmatRR_que[ibs_ck_to_look],
                                                      bmatAA_que[ibs_ck_to_look],
                                                      h11/2, h21/2,1,bpos);
                    while (previbs0 == -1 && ibs_ck_to_look > 0) {
                      ibs_ck_to_look--;
                      previbs0 = IBS0inOneBlock(bmatRR_que[ibs_ck_to_look],
                                                bmatAA_que[ibs_ck_to_look],
                                                h11/2, h21/2,1);
                    }
                    if (previbs0 > 0)
                      previbs0 = (*posvec_que[ibs_ck_to_look])[previbs0];
                    if (pgmap.bp2cm(nextibs0) - pgmap.bp2cm(previbs0) > delta) {
                      flag = 3;
                    }
                  }
                }
              }
              if (flag) { // Flip
                fliprec.push_back(std::vector<int32_t> {positions[k+1],h21/2,h11/2,
                                                       dij_p, dipj_s, dijp_s, flag});
                if (flipcandy.find(h21/2) != flipcandy.end()) {
                  flipcandy[h21/2].push_back(h11/2);
                } else {
                  flipcandy[h21/2] = std::vector<int32_t> {h11/2};
                }
              }
            } else if (dipj_s - dijp_s > diff && dipj_s - positions[k] > nmatch) { // Consider flip individual 1
              if (pgmap.bp2cm(dipj_s) - pgmap.bp2cm(dij_p) > delta) {
                flag = 1;
              } else {        // Need to evaluate no-ibs0
                int32_t bpos = bitpos / 8;  // This is not exact
                int32_t ibs_ck_to_look = cur_ibs_ck;
                int32_t nextibs0 = IBS0inOneBlock(bmatRR_que[ibs_ck_to_look],
                                                  bmatAA_que[ibs_ck_to_look],
                                                  h11/2, h21/2,0,bpos);
                while (nextibs0 == -1 &&
                       ibs_ck_to_look < ((int32_t) bmatRR_que.size())-1) {
                  ibs_ck_to_look++;
                  nextibs0 = IBS0inOneBlock(bmatRR_que[ibs_ck_to_look],
                                            bmatAA_que[ibs_ck_to_look],
                                            h11/2, h21/2,0);
                }
                nextibs0 = (nextibs0 > 0) ? (*posvec_que[ibs_ck_to_look])[nextibs0] : posvec_que[ibs_ck_to_look]->back();
                if (nextibs0 - positions[k] > lambda) {
                  // Not too close to ibs0
                  if (pgmap.bp2cm(nextibs0) - pgmap.bp2cm(dij_p) > delta) {
// std::cout << dij_p << '\t' << nextibs0 << '\t' << pgmap.bp2cm(dij_p) << '\t' << pgmap.bp2cm(nextibs0) << '\n';
                    flag = 2;
                  } else { // Need to check the previous ibs0
                    ibs_ck_to_look = cur_ibs_ck;
                    int32_t previbs0 = IBS0inOneBlock(bmatRR_que[ibs_ck_to_look],
                                                      bmatAA_que[ibs_ck_to_look],
                                                      h11/2, h21/2,1,bpos);
                    while (previbs0 == -1 && ibs_ck_to_look > 0) {
                      ibs_ck_to_look--;
                      previbs0 = IBS0inOneBlock(bmatRR_que[ibs_ck_to_look],
                                                bmatAA_que[ibs_ck_to_look],
                                                h11/2, h21/2,1);
                    }
                    if (previbs0 > 0)
                      previbs0 = (*posvec_que[ibs_ck_to_look])[previbs0];
                    if (pgmap.bp2cm(nextibs0) - pgmap.bp2cm(previbs0) > delta) {
                      flag = 3;
                    }
                  }
                }
              }
              if (flag) { // Flip
                fliprec.push_back(std::vector<int32_t> {positions[k+1],h11/2,h21/2,
                                                        dij_p, dijp_s, dipj_s, flag});
                if (flipcandy.find(h11/2) != flipcandy.end()) {
                  flipcandy[h11/2].push_back(h21/2);
                } else {
                  flipcandy[h11/2] = std::vector<int32_t> {h21/2};
                }
              }
            }
          }
          j++;
          if (prepc.d[j] > dij_p)
            dij_p = prepc.d[j];
        } // End of searching switch of pairs involving row i
      } // End of searching for this column
      // If multiple flips may occur, sort involved id by #occurence
      if (flipcandy.size() > 0) {
        std::map<int32_t, std::vector<int32_t> > finalrec;
        std::vector<bool> flipped(nsamples, 0);
        if (flipcandy.size() > 1) {
          std::set<std::pair<int32_t, std::vector<int32_t> >, Comparator> IDToFlip(flipcandy.begin(), flipcandy.end(), CompFn);
          for (auto & id : IDToFlip) {
            if ((int32_t) finalrec.size() >= max_flip) {
              break;
            }
            int32_t rm = (id.second).size();
            for (auto & id2 : id.second) {
              if (flipped[id2])
                rm--;
            }
            if (rm) {
              flipped[id.first] = 1;
              finalrec[id.first] = std::vector<int32_t>{0,id.first,0,0, 0,0,0};
            }
          }
        } else {
          flipped[flipcandy.begin()->first] = 1;
          finalrec[flipcandy.begin()->first] = std::vector<int32_t>{0,flipcandy.begin()->first,0,0, 0,0,0};
        }
        for (auto & v : fliprec) {
          if (flipped[v[1]]) {
            finalrec[v[1]][0] = v[0];
            finalrec[v[1]][2]++;
            if (v[5]-v[4] > finalrec[v[1]][3])
              finalrec[v[1]][3] = v[5]-v[4];
            finalrec[v[1]][3+v.back()]++;
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
        prepc.SwitchIndex(flipped);
        for (int32_t it = 0; it < nsamples; ++it) {
          if (flipped[it])
            flip_abs[it] = !flip_abs[it];
        }
      }
    } // End of processing one block
    wfs.close();
    ck++;
    st = ck * chunk_size + 1;
    stugly = positions.back();
    ed = st + chunk_size - 1;
    // Free memory of gt matrix
    for (auto& pt : gtmat)
      delete [] pt;
    // Free memory used by the suffix pbwt
    for (auto& pt : dmat)
      delete [] pt;
    for (auto& pt : rmat)
      delete [] pt;
    // Slide the window of ibs0 lookup
    cur_ibs_ck++;
// std::cout << st << '\t' << ed << '\t' << ibs_chunk_in_que.size() << '\n';

    while (ibs_chunk_in_que.size() > 1 && pgmap.bp2cm(stugly) - pgmap.bp2cm(posvec_que[0]->back()) > delta) {
// std::cout << pgmap.bp2cm(st) << '\t' << pgmap.bp2cm((*posvec_que[0]).back()) << '\n';
      delete posvec_que[0]; posvec_que.erase(posvec_que.begin());
      delete bmatRR_que[0]; bmatRR_que.erase(bmatRR_que.begin());
      delete bmatAA_que[0]; bmatAA_que.erase(bmatAA_que.begin());
      ibs_chunk_in_que.erase(ibs_chunk_in_que.begin());
      cur_ibs_ck--;
// std::cout << "After\t" << ibs_chunk_in_que.size() << '\t' <<  cur_ibs_ck << '\n';
    }
    ibs_ck = ibs_chunk_in_que.back() + 1;
    ibs_st = ibs_ck * chunk_size + 1;
    ibs_ed = (ibs_ck+1) * chunk_size;
    while (pgmap.bp2cm(posvec_que.back()->back()) - pgmap.bp2cm(ed) < delta) {
// std::cout << "Check\t" << ibs_st << '\t' <<  ibs_ed << '\n';
      if (ibs_st > pgmap.centromere_st && ibs_ed < pgmap.centromere_ed) {
        ibs_ck++;
        ibs_st = ibs_ck * chunk_size + 1;
        ibs_ed = (ibs_ck+1) * chunk_size;
        continue;
      }
      if (ibs_st > gpmap.maxpos)
        break;
      reg = chrom + ":" + std::to_string(ibs_st) + "-" + std::to_string(ibs_ed);
      if (ReadIBS0Block(reg, inVcf, bmatRR_que, bmatAA_que, posvec_que, min_hom_gts) > 0){
        ibs_chunk_in_que.push_back(ibs_st/chunk_size);
      }
      ibs_ck++;
      ibs_st = ibs_ck * chunk_size + 1;
      ibs_ed = (ibs_ck+1) * chunk_size;
    } // Finish adding new ibs0 lookup blocks

// if (pgmap.bp2cm(st) - pgmap.bp2cm((*posvec_que[0])[0]) < delta || pgmap.bp2cm(posvec_que.back()->back()) - pgmap.bp2cm(ed) < delta) {
// std::cout << "Finish adding new ibs0 lookup blocks\t" << ibs_chunk_in_que.size() << '\t' <<  cur_ibs_ck << '\n';
// std::cout << "First pos: " << (*posvec_que[0])[0] << ' ' << pgmap.bp2cm(st) - pgmap.bp2cm((*posvec_que[0])[0]) << "\tLast: " << posvec_que.back()->back() << ' ' << pgmap.bp2cm(posvec_que.back()->back()) - pgmap.bp2cm(ed) << '\n';
// }
// if (cur_ibs_ck >= 0)
//   std::cout << (*posvec_que[cur_ibs_ck])[0] << '\t' << (*posvec_que[cur_ibs_ck]).back() << '\n';

  } // Finish processing one chunk
  hts_close(wbcf);
  return 0;
}




