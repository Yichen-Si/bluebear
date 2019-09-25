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

int32_t pbwt_dist_suff(int32_t* d, int32_t i, int32_t j) {
  if (i > j) {
    int32_t tmp = j; j = i; i = tmp;
  }
  int32_t dist = d[i+1];
  for (int32_t it = i+1; it <= j; ++it) {
    if (d[it] < dist) {
      dist = d[it];
    }
  }
  return dist;
}

// Use haplotype matching length around rare alleles to detect blip

int32_t RareBlipOnlyForward(int32_t argc, char** argv) {
  std::string inVcf, inMap, path_pbwt, out, reg;
  std::string outf, outf_s, outVcf;
  int32_t detailed = 1;
  std::string chrom="chr20";
  // std::string mode ="b";
  int32_t verbose = 10000;
  int32_t min_variant = 2;
  int32_t nsamples=0, M=0;
  int32_t st, ck, ed, stugly;
  int32_t pbwt_chunk = 1000000;          // process suffix pbwt & ibs0 by chunk
  int32_t start_pos = 0; // Starting position

  int32_t rare_ac = 0;
  int32_t diff = 5000;
  double  mafcut = 0.0;

  int32_t* p_gt = NULL;
  int32_t  n_gt = 0;

  paramList pl;
  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input", NULL)
    LONG_STRING_PARAM("in-vcf",&inVcf, "Input VCF/BCF file")
    LONG_STRING_PARAM("region",&reg,"Genomic region to focus on")
    LONG_STRING_PARAM("chr",&chrom,"Chromosome to process")
    LONG_STRING_PARAM("map",&inMap, "Map file for genetic distance")
    LONG_STRING_PARAM("pbwt-path",&path_pbwt,"Pre-computed snapshot of suffix pbwt")

    LONG_PARAM_GROUP("Additional Options", NULL)
    LONG_DOUBLE_PARAM("rare-maf",&mafcut, "Sample minor allele frequency for considering blipping")
    LONG_INT_PARAM("rare-ac",&rare_ac,"Rare allele sharing as evidence for IBD")
    LONG_INT_PARAM("min-diff",&diff,"The minimal difference of new matches (in bp) between flipping two individuals for a flip to be recorded")
    LONG_INT_PARAM("start-from",&start_pos,"Starting position. Proceed on both sides.")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output file prefix")
    // LONG_STRING_PARAM("mode", &mode, "Output format of flipped bcf/vcf (b/bu/z)")
    LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")
    LONG_INT_PARAM("detailed",&detailed,"If detailed info for each pair is to be written (0/1)")

  END_LONG_PARAMS();
  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  // sanity check of input arguments
  if ( inVcf.empty() || out.empty() || inMap.empty() ) {
    error("[E:%s:%d %s] --in-vcf, --out, --map are required parameters",__FILE__,__LINE__,__FUNCTION__);
  }

  // Translate between genetic & physical position.
  bp2cmMap pgmap(inMap, " ");
  notice("Will proceed until the last position in linkage map: %d.", pgmap.maxpos);

  // Copy header & get sample size
  htsFile *fp = hts_open(inVcf.c_str(),"r");
  bcf_hdr_t *hdr = bcf_hdr_read(fp);
  // bcf_hdr_write(wbcf, hdr);
  nsamples = hdr->n[BCF_DT_SAMPLE];
  M = nsamples * 2;
  hts_close(fp);

  // Process rare allele count threshold
  if (rare_ac >  0) {mafcut = rare_ac*1.0/M;}
  if (rare_ac == 0) {rare_ac = (int32_t) (mafcut*M);}
  if (rare_ac == 0) {
    error("One of --rare-maf and --rare-ac is required.");
  }
  std::ostringstream parm;
  parm << std::fixed << std::setprecision(1) << mafcut*1000 << '-' << diff/1000;
  notice("Parameter setting: %s.", parm.str().c_str());

  // // Output flipped bcf/vcf
  // outVcf = out + "_" + parm.str() + "_st_" + std::to_string(start_pos) + "_fw.bcf";
  // if (mode=="z")
  //   outVcf = out + "_" + parm.str() + "_st_" + std::to_string(start_pos) + "_fw.vcf.gz";
  // mode = "w" + mode;
  // htsFile *wbcf = hts_open(outVcf.c_str(),mode.c_str());

  ck = start_pos/pbwt_chunk; // Chunk in terms of pbwt_chunk
  st = ck * pbwt_chunk + 1;
  ed = st + pbwt_chunk - 1;
  while (st > pgmap.centromere_st && ed < pgmap.centromere_ed ) {
    ck++;
    st = ck * pbwt_chunk + 1;
    ed = st + pbwt_chunk - 1;
  }
  stugly = st;

  // // Initialize ibs0 lookup blocks
  // reg = chrom + ":" + std::to_string(st) + "-" + std::to_string(ed);
  // RareIBS0lookup ibs0finder(inVcf, reg, pgmap, delta, rare_ac, ibs0_chunk, min_hom_gts);
  // notice("Initilized ibs0 lookup (%.1f cM, %d blocks).", delta, ibs0finder.start_que.size());

  // Initialize forward pbwt
  // Read the first snp in this chunk. should be stored as a checkpoint
  reg = chrom + ":" + std::to_string(st) + "-" + std::to_string(ed);
  std::string line;
  std::string sfile = path_pbwt + "_" + reg + "_prefix.amat";
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
  sf.open(sfile);
  std::getline(sf, line);
  sf.close();
  split(wvec, "\t", line);
  int32_t d[M];
  for (uint32_t i =  OFFSET; i < wvec.size(); ++i)
    d[i- OFFSET] = std::stoi(wvec[i]);

  pbwtCursor prepc(M, a, d);

  while (st < pgmap.maxpos) {
    // Read and build suffix pbwt by physical chunk
    // TODO: enable efficient random access &/ customizable chunk size

    // Read genotype matrix
    reg = chrom + ":" + std::to_string(stugly) + "-" + std::to_string(ed);
    // ibs0finder.Update(reg);
    std::vector<GenomeInterval> intervals;
    parse_intervals(intervals, "", reg);
    BCFOrderedReader odr(inVcf, intervals);
    bcf1_t* iv = bcf_init();

    std::vector<bool*> gtmat;
    std::vector<int32_t> positions;
    std::vector<double> mafs;
    std::vector<bool> oppo;

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

      if (iv->pos < start_pos) { // Start_pos & before is not flipped
        prepc.ForwardsAD_prefix(y, iv->pos+1);
        // bcf_write(wbcf, odr.hdr, iv);
        continue;
      }
      if (ac > M/2) {
        oppo.push_back(1);
        ac = M - ac;
      } else {
        oppo.push_back(0);
      }
      gtmat.push_back(y);
      positions.push_back(iv->pos+1);
      mafs.push_back(ac*1.0/M);
    } // Finish reading all haplotypes in this chunk
    int32_t N = positions.size();
    if ( N < min_variant ) {
      if (N > 0) {
        for (int32_t k = 0; k < N; ++ k)
          prepc.ForwardsAD_prefix(gtmat[k], positions[k]);
      }
      ck++;
      st = ck * pbwt_chunk + 1;
      stugly = st;
      ed = st + pbwt_chunk - 1;
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
    outf = out + "_" + parm.str() + "_st_" + std::to_string(start_pos) + "_fw_" + reg + ".blip.list";
    wf.open(outf, std::ios::trunc);
    std::ofstream wfs;
    outf_s = out + "_" + parm.str() + "_st_" + std::to_string(start_pos) + "_fw_" + reg + ".blip.list.short";
    wfs.open(outf_s, std::ios::trunc);

    // Read the last snp in this chunk. should be stored as a checkpoint
    sfile = path_pbwt + "_" + reg + "_suffix.amat";
    sf.open(sfile);
    std::getline(sf, line);
    sf.close();
    split(wvec, "\t", line);
    OFFSET = wvec.size() - M;
    for (uint32_t i =  OFFSET; i < wvec.size(); ++i)
      a[i- OFFSET] = std::stoi(wvec[i]);

    sfile = path_pbwt + "_" + reg + "_suffix.dmat";
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
    for (int32_t k = N-1; k >= 0; --k) {
      sufpc.ForwardsAD_suffix(gtmat[k], positions[k]);
      memcpy(dmat[k], sufpc.d, M*sizeof(int32_t));
      sufpc.ReverseA(rmat[k]);
    }
    // Build prefix pbwt & detect blip
    // if (stugly < start_pos) {
    // // If it is the first block, skip sites before the starting position.
    //   reg = chrom + ":" + std::to_string(start_pos) + "-" + std::to_string(ed);
    //   parse_intervals(intervals, "", reg);
    // }
    // odr.jump_to_interval(intervals[0]);
    // bcf_clear(iv);
    for (int32_t k = 0; k < N-2; ++k) {
      // Output current position
      // odr.read(iv);
      // if (bcf_get_genotypes(odr.hdr, iv, &p_gt, &n_gt) < 0) {
      //   error("Cannot find the field GT.");
      // }
      // int32_t y[M];
      // for (int32_t it = 0; it < nsamples; ++it) {
      //   y[it*2] = (flip_abs[it]) ? bcf_gt_phased(gtmat[k][it*2+1]) : bcf_gt_phased(gtmat[k][it*2]);
      //   y[it*2+1] = (flip_abs[it]) ? bcf_gt_phased(gtmat[k][it*2]) : bcf_gt_phased(gtmat[k][it*2+1]);
      // }
      // if(bcf_update_genotypes(odr.hdr, iv, y, M)) {
      //   error("Cannot update GT.");
      // }
      // bcf_write(wbcf, odr.hdr, iv);

      // Sorted upto and include position k
      prepc.ForwardsAD_prefix(gtmat[k], positions[k]);
      if (mafs[k+1] > mafcut) {continue;}
      bool minor = 1;
      if (oppo[k+1]) { // 0 labels the minor allele
        minor = 0;
      }
      // Find rare allele carriers
      std::vector<int32_t> carriers;
      std::vector<int32_t> prefindx;
      std::vector<int32_t> suffindx;
      for (int32_t i = 0; i < nsamples; ++i) {
        if ((gtmat[k+1][2*i]==minor) || (gtmat[k+1][2*i+1] == minor))
          carriers.push_back(i);
      }
      int32_t m = (int32_t) carriers.size();
      // Find index for correponding haps in pbwt matrices
      int32_t rvec[M];
      prepc.ReverseA(rvec);
      for (auto &v : carriers) {
        prefindx.push_back(rvec[2*v]);
        prefindx.push_back(rvec[2*v+1]);
        suffindx.push_back(rmat[k+2][2*v]);
        suffindx.push_back(rmat[k+2][2*v+1]);
      }
      // Find max. pairwise hap match
      // Could be more efficient.
      int32_t hapmat[prefindx.size()][prefindx.size()];
      int32_t max_dpref = positions[k+1], max_dsuff = positions[k+1];
      int32_t max_prefindx = 0, max_suffindx = 0;
      for (uint32_t r = 0; r < prefindx.size()-1; ++r) {
        for (uint32_t c = r+1; c < prefindx.size(); ++c) {
          // upper triangular is pref. distance (mis-match position)
          hapmat[r][c] = prepc.Dist_pref(prefindx[r], prefindx[c]);
          if (hapmat[r][c] < max_dpref) {
            max_dpref = hapmat[r][c];
            // max_prefindx = prefindx[r];
            max_prefindx = r;
          }
          // lower triangular is suff. distance (mis-match position)
          hapmat[c][r] = pbwt_dist_suff(dmat[k+2], suffindx[r], suffindx[c]);
          if (hapmat[c][r] > max_dsuff) {
            max_dsuff = hapmat[c][r];
            // max_suffindx = suffindx[r];
            max_suffindx = r;
          }
        }
      }
// std::cout << max_dpref << '\t' << max_dsuff << '\t' << positions[k+1] << "\n";

      // Use the longest hap. to detect blip
      int32_t left[2], right[2], hap[2];
      int32_t ip, lindx, ct = 0;
      for (int32_t it = 0; it < m; ++it) {
        for (int32_t jt = 0; jt < 2; ++jt) {
          ip = it*2+jt;
          if (ip == max_prefindx) {
            left[jt] = max_dpref;
          } else {
            left[jt] = hapmat[std::min(ip,max_prefindx)][std::max(ip,max_prefindx)];
          }
          if (ip == max_suffindx) {
            right[jt] = max_dsuff;
          } else {
            right[jt] = hapmat[std::max(ip,max_suffindx)][std::min(ip,max_suffindx)];
          }
          hap[jt] = right[jt]-left[jt];
// std::cout << left[jt] << '\t' << right[jt] << "\n";
        }

        lindx = (hap[0] > hap[1]) ? 0 : 1;
        if (hap[lindx] - hap[1-lindx] > diff && gtmat[k+1][carriers[it]*2+lindx] != 1) {
          ct++;
          gtmat[k+1][carriers[it]*2+1-lindx] = 0;
          gtmat[k+1][carriers[it]*2+lindx] = 1;
          std::stringstream recline;
          recline << positions[k+1] << '\t' << m << '\t' << carriers[it] << '\t' << hap[lindx] << '\t' << hap[1-lindx] << '\n';
          wf << recline.str();
// std::cout << positions[k+1] << '\t' << m << '\t' << carriers[it] << '\t' << hap[lindx] << '\t' << hap[1-lindx] << '\n';
        }
      }
      if (ct > 0) {
        std::stringstream recline;
        recline << positions[k+1] << '\t' << m << '\t' << ct << '\n';
        wfs << recline.str();
      }
    } // End of processing one block
    wf.close();
    wfs.close();
    ck++;
    st = ck * pbwt_chunk + 1;
    stugly = positions.back();
    ed = st + pbwt_chunk - 1;
    // Free memory of gt matrix
    for (auto& pt : gtmat)
      delete [] pt;
    // Free memory used by the suffix pbwt
    for (auto& pt : dmat)
      delete [] pt;
    for (auto& pt : rmat)
      delete [] pt;
  } // Finish processing one chunk
  // hts_close(wbcf);
  return 0;
}
