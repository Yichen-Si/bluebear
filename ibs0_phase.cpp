#include "cramore.h"
#include "bcf_filter_arg.h"
#include "bcf_ordered_reader.h"
#include "compact_matrix.h"

#include "pbwt_build.h"
#include "cm2bp.h"
#include "bp2cm.h"

#include <functional>

int32_t IBS0inOneBlock(bitmatrix* bmatRR, bitmatrix* bmatAA,
                       int32_t i, int32_t j,
                       bool reverse = 0, int32_t start = 0) {
  uint8_t* iRR = bmatRR->get_row_bits(i);
  uint8_t* iAA = bmatAA->get_row_bits(i);
  uint8_t* jRR = bmatRR->get_row_bits(j);
  uint8_t* jAA = bmatAA->get_row_bits(j);
  int32_t k;
  if (reverse) {
    if (!start)
      start = bmatRR->nbytes_col-1;
    for(k=start; k >= 0; --k) {
      if ( ( iRR[k] ^ jRR[k] ) & ( iAA[k] ^ jAA[k] ) ) { // IBS0 exists
        break;
      }
    }
    return k * 8 + 7;
  } else {
    for(k=start; k < bmatRR->nbytes_col; ++k) {
      if ( ( iRR[k] ^ jRR[k] ) & ( iAA[k] ^ jAA[k] ) ) { // IBS0 exists
        break;
      }
    }
    if (k == bmatRR->nbytes_col)
      return -1;
    return k * 8;
  }
}

int32_t IBS0Phase(int32_t argc, char** argv) {
  std::string inVcf, inMap, path_suffix, out, reg;
  std::string outf,outf_s;
  std::string chrom="chr20";
  int32_t verbose = 1000;
  int32_t min_variant = 2;
  int32_t min_hom_gts = 1;
  int32_t nsamples=0, M=0;
  int32_t st, ck, ed, ibs_st, ibs_ck, ibs_ed, stugly;
  int32_t chunk_size = 1000000;          // process suffix pbwt & ibs0 by chunk

  double  delta = 3.0;                   // thresholds
  int32_t lambda = 20000, gamma = 20000, diff = 5000;

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
    LONG_STRING_PARAM("pbwt-suffix",&path_suffix,"Pre-computed snapshot of suffix pbwt")

    LONG_PARAM_GROUP("Additional Options", NULL)
    LONG_DOUBLE_PARAM("delta",&delta, "no-ibs0 threshold (in cM) to be a proxy of IBD")
    LONG_INT_PARAM("lambda",&lambda,"Length (in bp) to the end of no-ibs0 that is not trusted")
    LONG_INT_PARAM("gamma",&gamma,"Length (in bp) of matched prefix that is requred to begin searching")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output file prefix")
    LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")

  END_LONG_PARAMS();
  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  // Translate between genetic & physical position.
  cm2bpMap gpmap(inMap, " ");
  bp2cmMap pgmap(inMap, " ");

  htsFile *fp = hts_open(inVcf.c_str(),"r");
  bcf_hdr_t *hdr = bcf_hdr_read(fp);
  nsamples = hdr->n[BCF_DT_SAMPLE];
  M = nsamples * 2;
  hts_close(fp);

  // To record flip
  std::vector<bool> flip(nsamples, 0);
  // Will sort flip candidates by #involved pairs
  typedef std::function<bool(std::pair<int32_t, std::vector<int32_t> >, std::pair<int32_t, std::vector<int32_t> >)> Comparator;
  Comparator CompFn =
    [](std::pair<int32_t, std::vector<int32_t> > elem1 ,std::pair<int32_t, std::vector<int32_t> > elem2) {
      if ((elem1.second).size() == (elem2.second).size())
        return elem1.first > elem2.first;
      else
        return (elem1.second).size() > (elem2.second).size();};

  ck = 0;
  st = ck * chunk_size + 1;
  stugly = st;
  ed = st + chunk_size - 1;

  // Initialize ibs0 lookup blocks. By 1Mb.
  // Number of blocks stored determined by genetic distance
  std::vector<bitmatrix*> bmatRR_que, bmatAA_que;
  std::vector<std::vector<int32_t>* > posvec_que;
  std::vector<int32_t> ibs_chunk_in_que;
  int32_t cur_ibs_ck = 0; // Currently processed position locates in the first block
  ibs_ck = 0;

  while (ibs_ck * chunk_size < gpmap.cm2bp(delta+pgmap.bp2cm(ed))) {
    ibs_st = ibs_ck * chunk_size + 1;
    ibs_ck++;
    ibs_ed = ibs_ck * chunk_size;
    reg = chrom + ":" + std::to_string(ibs_st) + "-" + std::to_string(ibs_ed);
    std::vector<GenomeInterval> ibs0interval;
    parse_intervals(ibs0interval, "", reg);
    BCFOrderedReader ibs0odr(inVcf, ibs0interval);
    bcf1_t* ibs0iv = bcf_init();
    bitmatrix *bmatRR = new bitmatrix(nsamples);
    bitmatrix *bmatAA = new bitmatrix(nsamples);
    std::vector<int32_t> *posvec = new std::vector<int32_t>;
    uint8_t* gtRR = (uint8_t*)calloc(nsamples, sizeof(uint8_t));
    uint8_t* gtAA = (uint8_t*)calloc(nsamples, sizeof(uint8_t));
    for(int32_t k=0; ibs0odr.read(ibs0iv); ++k) {  // read markers
      if ( bcf_get_genotypes(ibs0odr.hdr, ibs0iv, &p_gt, &n_gt) < 0 ) {
        error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(ibs0odr.hdr, ibs0iv->rid), ibs0iv->pos+1);
      }
      memset(gtRR, 0, nsamples);
      memset(gtAA, 0, nsamples);
      int32_t gcs[3] = {0,0,0};
      for(int32_t i=0; i < nsamples; ++i) {
        int32_t g1 = p_gt[2*i];
        int32_t g2 = p_gt[2*i+1];
        int32_t geno;
        if ( bcf_gt_is_missing(g1) || bcf_gt_is_missing(g2) ) {
          //geno = 0;
        } else {
          geno = ((bcf_gt_allele(g1) > 0) ? 1 : 0) + ((bcf_gt_allele(g2) > 0) ? 1 : 0);
          if ( geno == 0 )    { gtRR[i] = 1; }
          else if ( geno == 2 ) { gtAA[i] = 1; }
          ++gcs[geno];
        }
      }
      if (( gcs[0] >= min_hom_gts ) && ( gcs[2] >= min_hom_gts )) {
        bmatRR->add_row_bytes(gtRR);
        bmatAA->add_row_bytes(gtAA);
        posvec->push_back(ibs0iv->pos+1);
      }
    }
    free(gtRR);
    free(gtAA);
    bmatRR->transpose();
    bmatAA->transpose();
    bmatRR_que.push_back (bmatRR);
    bmatAA_que.push_back (bmatAA);
    posvec_que.push_back (posvec);
    ibs_chunk_in_que.push_back(ibs_ck-1);
  } // Finish initialize ibs0 lookup blocks

  // Initialize pbwt Cursor. Build prefix pbwt from begining.
  // TODO: start from a close snapshot
  pbwtCursor prepc(M, gpmap.minpos-1);

  while (ed < gpmap.maxpos) {
// if (st > (int32_t) 3e6)
//   break;
    // Read and build suffix pbwt by physical chunk
    // TODO: enable efficient random access &/ customizable chunk size
    reg = chrom + ":" + std::to_string(st) + "-" + std::to_string(ed);
    outf = out + "_flip_"+reg+".list";
    std::ofstream wf(outf, std::ios::trunc);
    outf_s = out + "_flip_"+reg+"_short.list";
    std::ofstream wfs(outf_s, std::ios::trunc);
    std::string line;
    std::string sfile = path_suffix + "_" + reg + "_suffix.amat";
    // Read the last snp in this chunk. should be stored as a checkpoint
    std::ifstream sf(sfile);
    std::getline(sf, line); // std::string line = GetLastLine(sfile);
    sf.close();
    std::vector<std::string> wvec;
    split(wvec, "\t", line);
    uint32_t OFFSET = wvec.size() - M;
    int32_t a[M];
    for (uint32_t i =  OFFSET; i < wvec.size(); ++i)
      a[i- OFFSET] = std::stoi(wvec[i]);

    sfile = path_suffix + "_" + reg + "_suffix.dmat";
    sf.open(sfile);
    std::getline(sf, line); // std::string line = GetLastLine(sfile);
    sf.close();
    split(wvec, "\t", line);
    int32_t d[M];
    for (uint32_t i =  OFFSET; i < wvec.size(); ++i)
      d[i- OFFSET] = std::stoi(wvec[i]);

    pbwtCursor sufpc(M, a, d);

    // Read genotype matrix
    reg = chrom + ":" + std::to_string(stugly) + "-" + std::to_string(ed);
    notice("Processing %s.", reg.c_str());
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
      notice("Read %d markers, start extending pbwt for %s.", N, reg.c_str());
    }

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
    // Build prefix pbwt & detect switch
    int32_t h11,h12,h21,h22; // (Reflecting lips so far) index stored in prefix pbwt
    int32_t i_s, j_s, i_s_prime, j_s_prime; // row num in suffix pbwt
    int32_t dij_p, dipj_s, dijp_s; // absolute position, from pbwt divergence matrix
    for (int32_t k = 0; k < N-1; ++k) {
      // Sorted upto and include position k
      prepc.ForwardsAD_prefix(gtmat[k], positions[k]);
      std::vector<std::vector<int32_t> > fliprec;
      // TODO: resolve multiple flips at one position
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
            if (dij_p < prepc.d[j])
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
            if (dijp_s - dipj_s > diff) { // Consider flip individual 2
              if (pgmap.bp2cm(dijp_s) - pgmap.bp2cm(dij_p) > delta) {
                flag = 1;
              } else {        // Need to evaluate no-ibs0
                int32_t bpos = k / 8;  // This is not exact
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
                nextibs0 = (*posvec_que[ibs_ck_to_look])[nextibs0];
                if (nextibs0 == -1 || nextibs0 - positions[k] > lambda) {
                  // Not too close to ibs0
                  if (pgmap.bp2cm(nextibs0) - pgmap.bp2cm(dij_p) > delta) {
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
                    if (pgmap.bp2cm(nextibs0) - pgmap.bp2cm(previbs0)) {
                      flag = 3;
                    }
                  }
                }
              }
              if (flag) { // Flip
                fliprec.push_back(std::vector<int32_t> {positions[k],h21/2,h11/2,
                                                       dij_p, dipj_s, dijp_s, flag});
                if (flipcandy.find(h21/2) != flipcandy.end()) {
                  flipcandy[h21/2].push_back(h11/2);
                } else {
                  flipcandy[h21/2] = std::vector<int32_t> {h11/2};
                }
              }
            } else if (dipj_s - dijp_s > diff) { // Consider flip individual 1
              if (pgmap.bp2cm(dipj_s) - pgmap.bp2cm(dij_p) > delta) {
                flag = 1;
              } else {        // Need to evaluate no-ibs0
                int32_t bpos = k / 8;  // This is not exact
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
                nextibs0 = (*posvec_que[ibs_ck_to_look])[nextibs0];
                if (nextibs0 == -1 || nextibs0 - positions[k] > lambda) {
                  // Not too close to ibs0
                  if (pgmap.bp2cm(nextibs0) - pgmap.bp2cm(dij_p) > delta) {
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
                    if (pgmap.bp2cm(nextibs0) - pgmap.bp2cm(previbs0)) {
                      flag = 3;
                    }
                  }
                }
              }
              if (flag) { // Flip
                fliprec.push_back(std::vector<int32_t> {positions[k],h11/2,h21/2,
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
          if (dij_p < prepc.d[j])
            dij_p = prepc.d[j];
        } // End of searching switch of pairs involving row i
      } // End of searching for this column
      // If multiple flips may occur, sort involved id by #occurence
      if (flipcandy.size() > 0) {
        std::set<int32_t> flipped;
        if (flipcandy.size() > 1) {
          std::set<std::pair<int32_t, std::vector<int32_t> >, Comparator> IDToFlip(flipcandy.begin(), flipcandy.end(), CompFn);
          for (auto & id : IDToFlip) {
            int32_t rm = (id.second).size();
            for (auto & id2 : id.second) {
              if (flipped.find(id2) != flipped.end())
                rm--;
            }
            if (rm) {
              flipped.insert(id.first);
            }
          }
        } else {
          flipped.insert(flipcandy.begin()->first);
        }
        std::map<int32_t, std::vector<int32_t> > finalrec;
        for (auto & v : flipped) {
          finalrec[v] = std::vector<int32_t>{0,v,0,0};
        }
        for (auto & v : fliprec) {
          if (flipped.find(v[1]) != flipped.end()) {
            finalrec[v[1]][0] = v[0];
            finalrec[v[1]][2]++;
            if (v[5]-v[4] > finalrec[v[1]][3])
              finalrec[v[1]][3] = v[5]-v[4];
            std::stringstream recline;
            for (auto& w : v)
              recline << w << '\t';
            recline.seekp(-1, std::ios_base::end);
            recline << '\n';
            wf << recline.str();
          }
        }
        for (auto & v : flipped) {
          std::stringstream recline;
          for (auto& w : finalrec[v])
            recline << w << '\t';
          recline.seekp(-1, std::ios_base::end);
          recline << '\n';
          wfs << recline.str();
        }
      }
    } // End of processing one block
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
    while (pgmap.bp2cm(st) - pgmap.bp2cm((*posvec_que[0]).back()) > delta) {
      delete posvec_que[0]; posvec_que.erase(posvec_que.begin());
      delete bmatRR_que[0]; bmatRR_que.erase(bmatRR_que.begin());
      delete bmatAA_que[0]; bmatAA_que.erase(bmatAA_que.begin());
      ibs_chunk_in_que.erase(ibs_chunk_in_que.begin());
      cur_ibs_ck--;
    }
    ibs_ck = ibs_chunk_in_que.back() + 1;

    while (pgmap.bp2cm(posvec_que.back()->back()) - pgmap.bp2cm(ed) < delta) {
      ibs_st = ibs_ck * chunk_size + 1;
      ibs_ck++;
      ibs_ed = ibs_ck * chunk_size;
      if (ibs_st > pgmap.centromere_st && ibs_ed < pgmap.centromere_ed)
        continue;
      if (ibs_st > gpmap.maxpos)
        break;
      reg = chrom + ":" + std::to_string(ibs_st) + "-" + std::to_string(ibs_ed);
      std::vector<GenomeInterval> ibs0interval;
      parse_intervals(ibs0interval, "", reg);
      BCFOrderedReader ibs0odr(inVcf, ibs0interval);
      bcf1_t* ibs0iv = bcf_init();
      bitmatrix *bmatRR = new bitmatrix(nsamples);
      bitmatrix *bmatAA = new bitmatrix(nsamples);
      std::vector<int32_t> *posvec = new std::vector<int32_t>;
      uint8_t* gtRR = (uint8_t*)calloc(nsamples, sizeof(uint8_t));
      uint8_t* gtAA = (uint8_t*)calloc(nsamples, sizeof(uint8_t));
      for(int32_t k=0; ibs0odr.read(ibs0iv); ++k) {  // read marker
        // extract genotype and apply genotype level filter
        if ( bcf_get_genotypes(ibs0odr.hdr, ibs0iv, &p_gt, &n_gt) < 0 ) {
          error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(ibs0odr.hdr, ibs0iv->rid), ibs0iv->pos+1);
        }
        memset(gtRR, 0, nsamples);
        memset(gtAA, 0, nsamples);
        int32_t gcs[3] = {0,0,0};
        for(int32_t i=0; i < nsamples; ++i) {
          int32_t g1 = p_gt[2*i];
          int32_t g2 = p_gt[2*i+1];
          int32_t geno;
          if ( bcf_gt_is_missing(g1) || bcf_gt_is_missing(g2) ) {
            //geno = 0;
          } else {
            geno = ((bcf_gt_allele(g1) > 0) ? 1 : 0) + ((bcf_gt_allele(g2) > 0) ? 1 : 0);
            if ( geno == 0 )    { gtRR[i] = 1; }
            else if ( geno == 2 ) { gtAA[i] = 1; }
            ++gcs[geno];
          }
        }
        if (( gcs[0] >= min_hom_gts ) && ( gcs[2] >= min_hom_gts )) {
          bmatRR->add_row_bytes(gtRR);
          bmatAA->add_row_bytes(gtAA);
          posvec->push_back(ibs0iv->pos+1);
        }
      }
      free(gtRR);
      free(gtAA);
      bmatRR->transpose();
      bmatAA->transpose();
      bmatRR_que.push_back (bmatRR);
      bmatAA_que.push_back (bmatAA);
      posvec_que.push_back (posvec);
      ibs_chunk_in_que.push_back(ibs_ck-1);
    } // Finish adding new ibs0 lookup blocks
    wf.close();
    wfs.close();
  } // Finish processing one chunk
  return 0;
}














