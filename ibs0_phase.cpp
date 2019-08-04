#include "cramore.h"
#include "bcf_filter_arg.h"
#include "bcf_ordered_reader.h"
#include "compact_matrix.h"

#include "pbwt_build.h"
#include "cm2bp.h"
#include "bp2cm.h"

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
  std::string outf;
  std::string chrom="chr20";
  int32_t verbose = 1000;
  int32_t min_variant = 2;
  int32_t min_hom_gts = 1;
  int32_t nsamples=0, M=0;
  int32_t st, ck, ed, ibs_st, ibs_ck, ibs_ed, stugly;
  int32_t chunk_size = 1000000;        // process suffix pbwt & ibs0 by chunk

  double  delta = 2.0;                   // thresholds
  int32_t lambda = 10000, gamma = 10000;

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
  // Output flip position
  std::vector<std::string> outlines;
  outf = out + "_flip.list";

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
    ibs_chunk_in_que.push_back(ibs_ck);
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
  } // Finish initialize ibs0 lookup blocks

  // Initialize pbwt Cursor. Build prefix pbwt from begining.
  // TODO: start from a close snapshot
  pbwtCursor prepc(M, gpmap.minpos-1);

  std::ofstream wf(outf);
  while (ed < gpmap.maxpos) {
if (st > (int32_t) 1e6)
  break;
    // Read and build suffix pbwt by physical chunk
    // TODO: enable efficient random access &/ customizable chunk size
    reg = chrom + ":" + std::to_string(st) + "-" + std::to_string(ed);
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
// int32_t avec[M];
    // Build suffix pbwt
    memcpy(dmat[N-1], sufpc.d, M*sizeof(int32_t));
    sufpc.ReverseA(rmat[N-1]);
    for (int32_t k = N-2; k >= 0; --k) {
      sufpc.ForwardsAD_suffix(gtmat[k], positions[k]);
      memcpy(dmat[k], sufpc.d, M*sizeof(int32_t));
      // if (k == 4)
        // memcpy(avec, sufpc.a, M*sizeof(int32_t));
      sufpc.ReverseA(rmat[k]);
    }

// // Check suffix matrix:
// notice("Checking suffix: \n");
// for (int32_t it=0; it<M; it++) {
//   for (int32_t jt=4; jt<14; jt++) {
//     std::cout << gtmat[jt][avec[it]] << ' ';
//   }
//   std::cout<<'\n';
// }
// std::cout<<'\n';

    // Build prefix pbwt & detect switch
    int32_t h11,h12,h21,h22; // (Reflecting current flips) index stored in prefix pbwt matrix
    int32_t i_s, j_s, i_s_prime, j_s_prime; // index in suffix pbwt
    int32_t dij_p, dipj_s, dijp_s; // absolute position, from pbwt divergence matrix
    for (int32_t k = 0; k < N-1; ++k) {
      // Sorted upto and include position k
      prepc.ForwardsAD_prefix(gtmat[k], positions[k]);
if (positions[k] == 146289||positions[k] == 146330) {
  std::cout << positions[k] << '\n';
  std::cout << "a: ";
  for (int32_t it=0; it<M; it++)
    std::cout << prepc.a[it] << ' ';
  std::cout << "\n";
  std::cout << "g: ";
  for (int32_t it=0; it<M; it++)
    std::cout << gtmat[k-1][prepc.a[it]] << ' ';
  std::cout << "\n";
  std::cout << "g: ";
  for (int32_t it=0; it<M; it++)
    std::cout << gtmat[k][prepc.a[it]] << ' ';
  std::cout << "\n";
  std::cout << "d: ";
  for (int32_t it=0; it<M; it++)
    std::cout << prepc.a[it] << ":" <<  prepc.d[it] << ' ';
  std::cout << "\n";
}

      std::set<int32_t> toswitch;

// int32_t CHECKPT=20;
// if (k == CHECKPT) {
//   notice("Checking prefix: \n");
//   int32_t avec[M];
//   memcpy(avec, prepc.a, M*sizeof(int32_t));
// std::cout<<'\n';
// for (int32_t it=0; it<M; it++) {
//   for (int32_t jt=CHECKPT-10; jt<=CHECKPT; jt++) {
//     std::cout << gtmat[jt][avec[it]] << ' ';
//   }
//   std::cout<<'\n';
// }
// std::cout<<'\n';
// for (int32_t it=0; it<nsamples; ++it) {
//   if (flip[it]) {
//     std::cout << "Flipped " << it << ";";
//     for (int32_t jt = 0; jt < M; ++jt) {
//       if (avec[jt] == it*2) {
//         std::cout << jt<< ";";
//         avec[jt] = it*2+1;}
//       else if (avec[jt] == it*2+1) {
//         avec[jt] = it*2;
//         std::cout << jt<< ";";}
//     }
//     std::cout<<'\t';
//   }
// }
// std::cout<<'\n';
// for (int32_t it=0; it<M; it++) {
//   for (int32_t jt=CHECKPT-10; jt<=CHECKPT; jt++) {
//     std::cout << gtmat[jt][avec[it]] << ' ';
//   }
//   std::cout<<'\n';
// }
// std::cout<<'\n';
// wf.close();
// exit(EXIT_FAILURE);
// }

      for (int32_t i = 0; i < M-1; ++i) {
        h11 = prepc.a[i]; // Haplotype ID. In original input data.
        h12 = h11 + 1 - 2 * (h11%2);
        // h11 = (flip[h11/2]) ? (h11 + (h11%2 == 0) ? 1 : -1) : h11;
        // h12 = h11 + 1 - 2 * (h11%2);
        int32_t j = i+1;
        dij_p = prepc.d[j];
        // Check a few hap down
        while (j < M && dij_p < positions[k] - gamma) {
          h21 = prepc.a[j];
          h22 = h21 + 1 - 2 * (h21%2);
          // h21 = (flip[h21/2]) ? (h21 + (h21%2 == 0) ? 1 : -1) : h21;
          // h22 = h21 + 1 - 2 * (h21%2);
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
            if (dijp_s > dipj_s) { // Consider flip individual 2
              if (pgmap.bp2cm(dijp_s) - pgmap.bp2cm(dij_p) > delta) {
                flag = 1;
              } else {        // Need to evaluate no-ibs0
                int32_t bpos = k / 8;  // TODO this is not accurate
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
                if (nextibs0 - positions[k] > lambda) { // Not too close to ibs0
                  if (pgmap.bp2cm(nextibs0) - pgmap.bp2cm(dij_p) > delta) {
                    flag = 2;
                  } else { // Need to check the previous ibs0
                    ibs_ck_to_look = cur_ibs_ck;
                    int32_t previbs0 = IBS0inOneBlock(bmatRR_que[ibs_ck_to_look],
                                                      bmatAA_que[ibs_ck_to_look],
                                                      h11/2, h21/2,1,bpos);
                    while (previbs0 == -1 && ibs_ck_to_look > 0) {
                      ibs_ck_to_look--;
                      nextibs0 = IBS0inOneBlock(bmatRR_que[ibs_ck_to_look],
                                                bmatAA_que[ibs_ck_to_look],
                                                h11/2, h21/2,1);
                    }
                    previbs0 = (*posvec_que[ibs_ck_to_look])[previbs0];
                    if (pgmap.bp2cm(nextibs0) - pgmap.bp2cm(previbs0)) {
                      flag = 3;
                    }
                  }
                }
              }
              if (flag) { // Flip
                flip[h21/2] = !flip[h21/2];
                // POS, ID to flip, ID as ref,
                // hap match before flipping & hap match after fipping
                std::vector<int32_t> rec{positions[k],i,j,h21/2, h11/2,
                  dij_p, dipj_s, dijp_s};
                std::stringstream recline;
                for (auto& v : rec)
                  recline << v << '\t';
                recline<< "Flip2 " << flag << '\t' << gtmat[k][h11]<<gtmat[k][h21] <<'\t'<<gtmat[k+1][h11]<<gtmat[k+1][h21]<<gtmat[k+1][h22] << '\n';
                // recline.seekp(-1, std::ios::end);
                // recline << '\n';
                wf << recline.str();
                toswitch.insert(h21/2);
                // prepc.SwitchHapIndex(h21, h22);
              }
            } else if (dijp_s < dipj_s) { // Consider flip individual 1
              if (pgmap.bp2cm(dipj_s) - pgmap.bp2cm(dij_p) > delta) {
                flag = 1;
              } else {        // Need to evaluate no-ibs0
                int32_t bpos = k / 8;  // TODO this is not accurate
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
                if (nextibs0 - positions[k] > lambda) { // Not too close to ibs0
                  if (pgmap.bp2cm(nextibs0) - pgmap.bp2cm(dij_p) > delta) {
                    flag = 2;
                  } else { // Need to check the previous ibs0
                    ibs_ck_to_look = cur_ibs_ck;
                    int32_t previbs0 = IBS0inOneBlock(bmatRR_que[ibs_ck_to_look],
                                                      bmatAA_que[ibs_ck_to_look],
                                                      h11/2, h21/2,1,bpos);
                    while (previbs0 == -1 && ibs_ck_to_look > 0) {
                      ibs_ck_to_look--;
                      nextibs0 = IBS0inOneBlock(bmatRR_que[ibs_ck_to_look],
                                                bmatAA_que[ibs_ck_to_look],
                                                h11/2, h21/2,1);
                    }
                    previbs0 = (*posvec_que[ibs_ck_to_look])[previbs0];
                    if (pgmap.bp2cm(nextibs0) - pgmap.bp2cm(previbs0)) {
                      flag = 3;
                    }
                  }
                }
              }
              if (flag) { // Flip
                flip[h11/2] = !flip[h11/2];
                // KeepGoDown = 0;
                std::vector<int32_t> rec{positions[k], i, j, h11/2, h21/2,
                  dij_p, dijp_s, dipj_s};
                std::stringstream recline;
                for (auto& v : rec)
                  recline << v << '\t';
                recline<< "Flip1 " << flag << '\t' << gtmat[k][h11]<<gtmat[k][h21] <<'\t'<<gtmat[k+1][h11]<<gtmat[k+1][h21]<<gtmat[k+1][h12] << '\n';
                // recline.seekp(-1, std::ios::end);
                // recline << '\n';
                wf << recline.str();
                // prepc.SwitchHapIndex(h11, h12);
                toswitch.insert(h11/2);
                break;
              }
            }
          }
          j++;
          if (dij_p < prepc.d[j])
            dij_p = prepc.d[j];
        } // End of searching switch of pairs involving row i
      } // End of searching for this column
      if (toswitch.size() > 0) {
        wf << "Flip happening: ";
        for (auto it = toswitch.begin(); it != toswitch.end(); ++it) {
          wf << *it << '\t';
          prepc.SwitchHapIndex((*it)*2, (*it)*2+1);
        }
      wf <<'\n';
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
    while (pgmap.bp2cm(st) - pgmap.bp2cm((*posvec_que[0]).back()) > delta) {
      delete posvec_que[0]; posvec_que.erase(posvec_que.begin());
      delete bmatRR_que[0]; bmatRR_que.erase(bmatRR_que.begin());
      delete bmatAA_que[0]; bmatAA_que.erase(bmatAA_que.begin());
      ibs_chunk_in_que.erase(ibs_chunk_in_que.begin());
      cur_ibs_ck++;
    }
    ibs_ck = ibs_chunk_in_que.back() + 1;

    while (pgmap.bp2cm(posvec_que.back()->back()) - pgmap.bp2cm(ed) < delta) {
      ibs_chunk_in_que.push_back(ibs_ck);
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
    } // Finish adding new ibs0 lookup blocks
  }
  return 0;
}














