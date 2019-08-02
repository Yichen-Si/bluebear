#include "cramore.h"
#include "bcf_filter_arg.h"
#include "bcf_ordered_reader.h"

#include "pbwt_build.h"

void pbwtSuffix(pbwtCursor& pc, std::vector<bool*>& gtmat, std::vector<int32_t>& posvec, std::vector<int32_t*>& dmat, std::vector<int32_t*>& rmat) {
  int32_t N = (int32_t) posvec.size();
  for (int32_t k = N-1; k >= 0; --k) {
    pc.ForwardsAD(gtmat[k], posvec[k]);
    memcpy(dmat[k], pc.d, pc.M*sizeof(int32_t));
    pc.ReverseA(rmat[k]);
  }
}

// Input bcf file
// Output suffix pbwt by chunk
int32_t pbwtBuildSuffix(int32_t argc, char** argv) {

  std::string inVcf, reg, out;
  int32_t chunksize = 100000;   // Read genotype by chunk
  int32_t store_interval = 500; // Store snapshot of pbwt
  int32_t nsamples=0, M=0;
  int32_t min_variant = 1;

  paramList pl;
  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Sites", NULL)
    LONG_STRING_PARAM("in-vcf",&inVcf, "Input VCF/BCF file")
    LONG_STRING_PARAM("region",&reg,"Genomic region to focus on")

    LONG_PARAM_GROUP("Additional Options", NULL)
    LONG_INT_PARAM("chunk-size",&chunksize, "Read, process and write by chunk")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output file prefix")
    LONG_INT_PARAM("store-per-markers",&store_interval, "Snapshot per x markers")

  END_LONG_PARAMS();
  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  // std::string inVcf = "/net/wonderland/home/ycsi/IBD/ibs0/Data/test/test.10.bcf";
  // std::string reg = "chr20:3100000-3200000";

  int32_t start = 0, end = 300000000;
  std::vector<std::string> v;
  std::string chrom = "chr20";
  if (!reg.empty()) {
    split(v, ":-", reg);
    if (!str2int32(v[1], start) || !str2int32(v[2], end)) {
      error("Invalid region.");
    }
    chrom = v[0];
  }

  int32_t nchunk = (end-start)/chunksize;
  int32_t ck=nchunk;
  int32_t st, ed;
  for (ck = nchunk; ck >= 0; --ck) {
    st = start + ck*chunksize;
    ed = st + chunksize;
    reg = chrom + ":" + std::to_string(st) + "-" + std::to_string(ed);
    std::vector<GenomeInterval> intervals;
    parse_intervals(intervals, "", reg);
    BCFOrderedReader odr(inVcf, intervals);
    bcf1_t* iv = bcf_init();
    if (odr.read(iv)) {
      nsamples = bcf_hdr_nsamples(odr.hdr);
      M = nsamples * 2;
      break;
    }
  }
  if (M == 0) {
    error("Invalid --region or unreadable --inVcf.");
  }

  pbwtCursor pc(M);
  nchunk = ck;

  notice("Started Reading VCF, identifying %d samples. Will procees in %d chunks and store snapshot every %d markers.", nsamples, nchunk+1, store_interval);

  for (ck = nchunk; ck >= 0; --ck) {

    st = start + ck*chunksize;
    ed = st + chunksize;
    if (st == end) continue;
    if (ed > end) ed = end;
    reg = chrom + ":" + std::to_string(st) + "-" + std::to_string(ed);
    notice("Processing %s.", reg.c_str());
    std::vector<GenomeInterval> intervals;
    parse_intervals(intervals, "", reg);
    BCFOrderedReader odr(inVcf, intervals);
    bcf1_t* iv = bcf_init();

    std::vector<bool*> gtmat;
    std::vector<int32_t> positions;

    int32_t *p_gt = NULL;
    int32_t  n_gt = 0;

    // read marker
    for (int32_t k=0; odr.read(iv); ++k) {
      // if (k > 10) {break;}
      if (bcf_get_genotypes(odr.hdr, iv, &p_gt, &n_gt) < 0) {
        error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);
      }

      bool *y;
      y = new bool[M];
      for (int32_t i = 0; i <  nsamples; ++i) {
        int32_t g1 = p_gt[2*i];
        int32_t g2 = p_gt[2*i+1];
        if (bcf_gt_is_missing(g1)) {
          y[2*i] = 0;
        } else {
          y[2*i] = ((bcf_gt_allele(g1) > 0) ? 1 : 0);
        }
        if (bcf_gt_is_missing(g2)) {
          y[2*i+1] = 0;
        } else {
          y[2*i+1] = ((bcf_gt_allele(g2) > 0) ? 1 : 0);
        }
      }
      gtmat.push_back(y);
      positions.push_back(iv->pos+1);
    }

    int32_t N = positions.size();
    if ( N < min_variant ) {
      notice("Observed only %d informative markers. Skipping this chunk...", N);
      continue;
    } else {
      notice("Read %d markers, start extending pbwt for %s.", N, reg.c_str());
    }

    // // Allocate space for pbwt matricies
    // std::vector<int32_t*> dmat(N,NULL);
    // std::vector<int32_t*> amat(N,NULL);
    // // std::vector<int32_t*> rmat(N,NULL);
    // for (int32_t i = 0; i < N; ++i) {
    //   dmat[i] = new int32_t[M];
    //   amat[i] = new int32_t[M];
    // }

    // Build pbwt (suffix) for this block

    std::string dout = "", aout = "";
    for (int32_t k = N-1; k >= 0; --k) {
      pc.ForwardsAD(gtmat[k], positions[k]);

      if (k % store_interval == 0) {
        // memcpy(dmat[k], pc.d, pc.M*sizeof(int32_t));
        // memcpy(amat[k], pc.a, pc.M*sizeof(int32_t));
        dout += (chrom + '\t' + std::to_string(positions[k]));
        for (int32_t j = 0; j < M; ++j)
          dout += ('\t' + std::to_string(pc.d[j]));
        dout += '\n';
        aout += (chrom + '\t' + std::to_string(positions[k]));
        for (int32_t j = 0; j < M; ++j)
          aout += ('\t' + std::to_string(pc.a[j]));
        aout += '\n';
      }
    }

    // pbwtSuffix(pc, gtmat, positions, dmat, rmat);

    // Write to file
    std::string outf = out + "_" + reg + ".dmat";
    htsFile* wf = hts_open(outf.c_str(), "w");
    hprintf(wf, "%s", dout.c_str());
    hts_close(wf);

    outf = out + "_" + reg + ".amat";
    wf = hts_open(outf.c_str(), "w");
    hprintf(wf, "%s", aout.c_str());
    hts_close(wf);

    // for (int32_t i = 0; i < N; ++i) {
    //   delete [] gtmat[i];
    //   delete [] dmat[i];
    //   delete [] rmat[i];
    // }

  }

  // // Print out test results. Sorted from the first position
  // for (uint32_t i = 0; i < positions.size(); ++i) {
  //   for (int32_t j = 0; j < M; ++j) {
  //     std::cout << gtmat[i][pc.a[j]] << ' ';
  //   }
  //   std::cout << std::endl;
  // }
  // std::cout << std::endl;
  // for (int32_t j = 0; j < M; ++j) {
  //   std::cout << pc.d[j]-positions[0] << ' ';
  // }


  // // Print out test results. Sorted from the first position
  // // Same as above, check the reverse index
  // for (int32_t i = 0; i < 10; ++i) {
  //   for (int32_t j = 0; j < M; ++j) {
  //     for (int32_t k = 0; k < M; ++k) {
  //       if (rmat[0][k] == j) {
  //         std::cout << gtmat[i][k] << ' ';
  //       }
  //     }
  //   }
  //   std::cout << std::endl;
  // }

  return 0;
}
















