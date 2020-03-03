#include "cramore.h"
#include "utils.h"
#include "bcf_filter_arg.h"
#include "bcf_ordered_reader.h"
#include "bcf_ordered_writer.h"
#include "pbwt_build.h"
#include "bp2cm.h"
#include "rare_variant_ibs0.h"

void UpdateInfo_IBD(bcf_hdr_t *hdr, bcf1_t *nv, RareVariant *rare) {
  std::stringstream ss;
  ss << rare->subset1.size() << "," << rare->subset2.size();
  bcf_update_info_string(hdr, nv, "BiCluster", ss.str().c_str());
  bcf_update_info_int32(hdr, nv, "AvgDist_bp", &(rare->AvgDist), 1);
  bcf_update_info_float(hdr, nv, "AvgDist_cM", &(rare->AvgDist_cm), 1);
  bcf_update_info_float(hdr, nv, "MedDist_cM", &(rare->MedDist_cm), 1);
}

// goal -- get all pairwise ibd regions covering shared
//         rare variants from the given region (for simulated data)
int32_t AnnotateIBDAroundRare(int32_t argc, char** argv) {
  std::string inVcf, inMap, chrom, reg, oreg;
  std::string out;
  int32_t verbose = 10000;
  int32_t max_rare_ac = 10;
  int32_t start, end;
  int32_t min_pos = -1, max_pos = -1;
  int32_t chunksize = 1000000;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Sites", NULL)
    LONG_STRING_PARAM("in-vcf",&inVcf, "Input VCF/BCF file")
    LONG_STRING_PARAM("region",&reg,"Genomic region to focus on")
    LONG_STRING_PARAM("map",&inMap, "Map file for genetic distance")

    LONG_PARAM_GROUP("Additional Options", NULL)
    LONG_INT_PARAM("min-pos",&min_pos, "The start point to search for ibd (bp)")
    LONG_INT_PARAM("max-pos",&max_pos, "The end point to search for ibd (bp)")

    LONG_INT_PARAM("max-rare",&max_rare_ac, "Maximal minor allele count to be considered as anchor for IBD")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output VCF file name")
    LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")

  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  // sanity check of input arguments
  if ( inVcf.empty() || reg.empty() || inMap.empty() || out.empty() ) {
    error("[E:%s:%d %s] --in-vcf, --region, --inMap, --out are required parameters",__FILE__,__LINE__,__FUNCTION__);
  }

  // Genetic map
  bp2cmMap pgmap(inMap, " ", "", -1, -1);
  notice("Read map. min %d; max %d; cst %d; ced %d.", pgmap.minpos, pgmap.maxpos, pgmap.centromere_st, pgmap.centromere_ed);
  if (min_pos < 0)
    min_pos = pgmap.minpos;
  if (max_pos < 0)
    max_pos = pgmap.maxpos;

  // Region to process
  oreg = reg;
  std::vector<std::string> v;
  split(v, ":-", reg);
  chrom = v[0];
  if (!str2int32(v[1], start) || !str2int32(v[2], end)) {
    error("Invalid region.");
  }

  // bcf reader
  // TODO: currently assume a small enough region.
  std::vector<GenomeInterval> intervals;
  parse_intervals(intervals, "", reg);
  BCFOrderedReader odr(inVcf, intervals);
  bcf1_t* iv = bcf_init();

  // bcf writer (site only)
  BCFOrderedWriter odw(out.c_str(),0);
  bcf_hdr_t* hnull = bcf_hdr_subset(odr.hdr, 0, 0, 0);
  bcf_hdr_remove(hnull, BCF_HL_FMT, NULL);
  odw.set_hdr(hnull);
  char buffer[65536];
  if ( bcf_hdr_id2int(odr.hdr, BCF_DT_ID, "BiCluster" ) < 0 ) {
    sprintf(buffer,"##INFO=<ID=BiCluster,Number=1,Type=String,Description=\"Sizes of the two subtree\">\n");
    bcf_hdr_append(odw.hdr, buffer);
  }
  if ( bcf_hdr_id2int(odr.hdr, BCF_DT_ID, "AvgDist_bp" ) < 0 ) {
    sprintf(buffer,"##INFO=<ID=AvgDist_bp,Number=1,Type=Integer,Description=\"Average no-IBS0 length between two clusters, in bp\">\n");
    bcf_hdr_append(odw.hdr, buffer);
  }
  if ( bcf_hdr_id2int(odr.hdr, BCF_DT_ID, "AvgDist_cM" ) < 0 ) {
    sprintf(buffer,"##INFO=<ID=AvgDist_cM,Number=1,Type=Float,Description=\"Average no-IBS0 length between two clusters, in cM\">\n");
    bcf_hdr_append(odw.hdr, buffer);
  }
  if ( bcf_hdr_id2int(odr.hdr, BCF_DT_ID, "MedDist_cM" ) < 0 ) {
    sprintf(buffer,"##INFO=<ID=MedDist_cM,Number=1,Type=Float,Description=\"Median no-IBS0 length between two clusters, in cM\">\n");
    bcf_hdr_append(odw.hdr, buffer);
  }
  odw.write_hdr();

  int32_t nVariant = 0, nFinished = 0;
  int32_t nsamples = bcf_hdr_nsamples(odr.hdr);
  int32_t M = nsamples * 2;

  notice("Started Reading site information from VCF file, identifying %d samples", nsamples);
  if ( nsamples == 0 )
    error("FATAL ERROR: The VCF does not have any samples with genotypes");

  std::map< std::pair<int32_t, int32_t>, std::vector<RareVariant*> > idpair_l;
  std::map< std::pair<int32_t, int32_t>, std::vector<RareVariant*> > idpair_r;
  std::vector<RareVariant*> snplist;
  int32_t* p_gt = NULL;
  int32_t  n_gt = 0;
  int32_t *info_ac = NULL;
  int32_t  n_ac = 0, pos = 0;
  int32_t ac = 0;

  // First pass.
  // Find indiv. carrying rare alleles.
  for(int32_t k=0; odr.read(iv); ++k) {
    if (iv->pos == pos) {
      if (iv->pos == snplist.back()->iv->pos) {
        delete snplist.back();
        snplist.pop_back();
        nVariant--;
      }
      continue; // Ignore tri-allelic sites
    }
    // periodic message to user
    if ( k % verbose == 0 )
      notice("Processing %d markers at %s:%d. Recorded %d rare variants", k, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1, nVariant);

    if ( bcf_get_genotypes(odr.hdr, iv, &p_gt, &n_gt) < 0 ) {
      error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);
    }
    if (bcf_get_info_int32(odr.hdr, iv, "AC", &info_ac, &n_ac) < 0) {continue;}
    if (info_ac[0] < 2) {continue;}
    ac = info_ac[0];
    if (ac < 2 || ac > max_rare_ac) {continue;}

    std::vector<int32_t> carry;
    for(int32_t i=0; i < nsamples; ++i) {
      int32_t g1 = p_gt[2*i];
      int32_t g2 = p_gt[2*i+1];
      int32_t geno;
      if ( bcf_gt_is_missing(g1) || bcf_gt_is_missing(g2) ) {
        continue;
      } else {
        geno = ((bcf_gt_allele(g1) > 0) ? 1 : 0) + ((bcf_gt_allele(g2) > 0) ? 1 : 0);
        if (geno == 1) {
          carry.push_back(i);
        }
      }
    } // Found carriers (have to be het)
    // TODO: here we ignore rare varaiants with any hom carrier
    if ((int32_t) carry.size() != ac) {continue;}

    // Creat a new RareVariant record
    pos = iv->pos;
    bcf1_t* nv = bcf_dup(iv);
    bcf_unpack(nv, BCF_UN_INFO);
    bcf_subset(odw.hdr, nv, 0, 0);
    RareVariant* rare = new RareVariant(nv, ac);
    rare->Add_id(carry);
    snplist.push_back(rare);
    nVariant++;
  }

  notice("Finished first pass. Processed %d rare variants across %d samples; finished %d.", nVariant, nsamples, nFinished);
  if (nVariant < 1) {
    notice("Not enough rare variants in the given region");
    odw.close();
    return 0;
  }

  // Build pbwt from beginning.

  reg = chrom + ":" + std::to_string(min_pos) + "-" + std::to_string(end+1);
  parse_intervals(intervals, "", reg);
  odr.jump_to_interval(intervals[0]);
  bcf_clear(iv);

  // Initialize pbwt Cursor
  pbwtCursor pc(M, min_pos);

  RareVariant *rv;

  // Build forward pbwt
  int32_t pt = 0; // index of the next rara variant in snplist
  for (int32_t k=0; odr.read(iv); ++k) {
    if ( k % verbose == 0 )
      notice("Forward pbwt. Processing %d markers at %s:%d.", k, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);
    if (bcf_get_genotypes(odr.hdr, iv, &p_gt, &n_gt) < 0) {
      error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);
    }
    if (bcf_get_info_int32(odr.hdr, iv, "AC", &info_ac, &n_ac) < 0) {continue;}
    ac = info_ac[0];

    if ( iv->pos == (snplist[pt]->iv)->pos ) {
      // Check haplotype matching before this position
      if ( ((int32_t) snplist[pt]->id_list.size() ) == ac ) {
        rv = snplist[pt];
        std::vector<int32_t> hapcarry; // which haplotype contains the rare allele
        for ( auto & v : rv->id_list ) {
          int32_t g1 = p_gt[2*v];
          int32_t g2 = p_gt[2*v+1];
          if (bcf_gt_allele(g1) > 0 || bcf_gt_allele(g2) > 0) {
            hapcarry.push_back((bcf_gt_allele(g1) > 0)? 0 : 1);
          } else { // Never happen
            error("Error in identifying rare variant carriers");
          }
        }
        int32_t rvec[M];
        pc.ReverseA(rvec);
        int32_t h1,h2,fin;
        for (int32_t i = 0; i < ac-1; i++) {
          h1 = rvec[ rv->id_list[i] * 2 + hapcarry[i]];
          for (int32_t j = i+1; j < ac; j++) {
            h2 = rvec[ rv->id_list[j] * 2 + hapcarry[j]];
// std::cout << iv->pos << "\tLeft: " << pc.Dist_pref(h1,h2) << '\n';
            fin = rv->Add_half( rv->id_list[i], rv->id_list[j], pc.Dist_pref(h1,h2), 1);
          }
        }
      }
      pt++;
      if (pt >= (int32_t) snplist.size()) {break;}
    }
    // Update pbwt cursor
    if (ac < max_rare_ac) {continue;}
    bool y[M];
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
    pc.ForwardsAD_prefix(y, iv->pos+1);

  }

  // Build backward pbwt
  pt = snplist.size() - 1;

  int32_t nchunk = (max_pos-start)/chunksize;
  int32_t ck=nchunk;
  int32_t st, ed;

  for (ck = nchunk; ck >= 0; --ck) {
    st = start + ck*chunksize + 1;
    ed = start + ck*chunksize + chunksize;
    reg = chrom + ":" + std::to_string(st) + "-" + std::to_string(ed);
    parse_intervals(intervals, "", reg);
    odr.jump_to_interval(intervals[0]);
    bcf_clear(iv);
    if (odr.read(iv)) {
      break;
    }
  }
  nchunk = ck + 1;

  // Initialize pbwt Cursor
  pc.Reset(M, max_pos);

  // Have to read by chunk
  for (ck = 0; ck <= nchunk; ++ck) {

    ed = max_pos - ck*chunksize;
    st = ed  - chunksize + 1;
    if (ed <= start) {break;}
    if (st < start) st = start;
    reg = chrom + ":" + std::to_string(st) + "-" + std::to_string(ed);
    notice("Processing %s.", reg.c_str());
    parse_intervals(intervals, "", reg);
    odr.jump_to_interval(intervals[0]);
    bcf_clear(iv);

    std::vector<bool*> gtmat;
    std::vector<int32_t> positions;
    std::vector<int32_t> mac;

    int32_t *p_gt = NULL;
    int32_t  n_gt = 0;
    int32_t *info_ac = NULL;
    int32_t n_ac = 0;

    // read marker
    for (int32_t k=0; odr.read(iv); ++k) {
      if (bcf_get_genotypes(odr.hdr, iv, &p_gt, &n_gt) < 0) {
        error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);
      }
      if (bcf_get_info_int32(odr.hdr, iv, "AC", &info_ac, &n_ac) < 0) {continue;}
      ac = info_ac[0];
      if (ac < 2) {continue;}
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
      positions.push_back(iv->pos);
      mac.push_back(ac);
    }

    int32_t N = positions.size();
    if ( N < 1 ) {
      notice("Observed 0 informative markers. Skipping this chunk...");
      continue;
    } else {
      notice("Read %d markers, start backward pbwt for %s.", N, reg.c_str());
    }

    for (int32_t k = N-1; k >= 0; --k) {
      ac = mac[k];
      if (positions[k] == snplist[pt]->iv->pos) {
        // Check haplotype matching before this position
        if ( ((int32_t) snplist[pt]->id_list.size() ) == ac ) {
          rv = snplist[pt];
          std::vector<int32_t> hapcarry; // which haplotype contains the rare allele
          for ( auto & v : rv->id_list ) {
            int32_t g1 = gtmat[k][2*v];
            int32_t g2 = gtmat[k][2*v+1];
            if (g1|| g2) {
              hapcarry.push_back((g1)? 0 : 1);
            } else { // Never happen
              error("Error in identifying rare variant carriers");
            }
          }
          int32_t rvec[M];
          pc.ReverseA(rvec);
          int32_t h1,h2,fin;
          for (int32_t i = 0; i < ac-1; i++) {
            h1 = rvec[ rv->id_list[i] * 2 + hapcarry[i]];
            for (int32_t j = i+1; j < ac; j++) {
              h2 = rvec[ rv->id_list[j] * 2 + hapcarry[j]];
              fin = rv->Add_half( rv->id_list[i], rv->id_list[j], pc.Dist_suff(h1,h2), 2);
            }
          }
        }
        pt--;
        if (pt < 0) {break;}
      }
      if (mac[k] < max_rare_ac) {continue;}
      pc.ForwardsAD_suffix(gtmat[k], positions[k]);
    } // Finish building pbwt for this block
    for (int32_t k = 0; k < N; ++k) {
      delete [] gtmat[k];
    }
  }

  notice("Finish finding IBD. Start output");
  // Process & Output
  for (auto & ptr : snplist) {
    if (ptr->IfDone()) {
      // Build tree
      ptr->Organize(pgmap);
      // Output this variant
      bcf1_t* nv = bcf_dup(ptr->iv);
      bcf_unpack(nv, BCF_UN_INFO);
      UpdateInfo_IBD(odw.hdr, nv, ptr);
      odw.write(nv);
      bcf_destroy(nv);
    } else {
      odw.write(ptr->iv);
    }
    // Delete the record
    delete ptr;
    ptr = NULL;
  }

  odw.close();

  return 0;
}










