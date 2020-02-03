#include "cramore.h"
#include "utils.h"
#include "bcf_filter_arg.h"
#include "bcf_ordered_reader.h"
#include "bcf_ordered_writer.h"

// goal -- (For simulated data) merge nearby rare variants to create artificial parallel mutation
int32_t VCFInsertPM(int32_t argc, char** argv) {
  std::string inVcf, outVcf, reg;
  // std::string out, outf;
  int32_t verbose = 5000;
  int32_t maxac = 10;
  int32_t winsize = 200;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Files", NULL)
    LONG_STRING_PARAM("in-vcf",&inVcf, "Input VCF/BCF file")
    LONG_STRING_PARAM("out-vcf",&outVcf, "Output VCF/BCF file")
    LONG_STRING_PARAM("region",&reg,"Genomic region to focus on")

    LONG_PARAM_GROUP("Additional Options", NULL)
    LONG_INT_PARAM("max-ac",&maxac, "Maximal allele count to consider as focal rare variants")
    LONG_INT_PARAM("win-size",&winsize, "Maximal window size to merge as one locus")

    LONG_PARAM_GROUP("Output Options", NULL)
    // LONG_STRING_PARAM("out", &out, "Output file prefix")
    LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")

  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  // sanity check of input arguments
  if ( inVcf.empty() || outVcf.empty() ) {
    error("[E:%s:%d %s] --in-vcf, --outVcf are required parameters",__FILE__,__LINE__,__FUNCTION__);
  }

  // bcf reader
  std::vector<GenomeInterval> intervals;
  parse_intervals(intervals, "", reg);
  BCFOrderedReader odr(inVcf, intervals);
  bcf1_t* iv = bcf_init();

  // Output modified bcf file
  BCFOrderedWriter odw(outVcf.c_str(), 0);
  odw.set_hdr(odr.hdr);

  int32_t nsamples = bcf_hdr_nsamples(odr.hdr);
  notice("Identifying %d samples from Vcf file", nsamples);
  if ( nsamples == 0 )
    error("FATAL ERROR: The VCF does not have any samples with genotypes");

   // Add info field indicating pm
  char buffer[65536];
  if ( bcf_hdr_id2int(odw.hdr, BCF_DT_ID, "PMconfig" ) < 0 ) {
    sprintf(buffer,"##INFO=<ID=PMconfig,Number=.,Type=Integer,Description=\"AC for each of the original variants merged to the new record\">\n");
    bcf_hdr_append(odw.hdr, buffer);
  }
  if ( bcf_hdr_id2int(odw.hdr, BCF_DT_ID, "PMid" ) < 0 ) {
    sprintf(buffer,"##INFO=<ID=PMid,Number=.,Type=Integer,Description=\"Carrier id for original variants merged to the new record\">\n");
    bcf_hdr_append(odw.hdr, buffer);
  }
  bcf_hdr_sync(odw.hdr);
  odw.write_hdr();

  // Scan upto window size
  int32_t *p_gt = NULL;
  int32_t  n_gt = 0;
  int32_t *info_ac = NULL;
  int32_t  n_ac = 0;
  int32_t prepos = 0, pos = 0;
  std::vector<bcf1_t* > ivstack;
  std::vector<int32_t>  acstack;
  int32_t acsum = 0;

  while (odr.read(iv)) {
    if (bcf_get_info_int32(odr.hdr, iv, "AC", &info_ac, &n_ac) < 0) {
      continue;
    }
    break;
  }
  odw.write(iv);
  pos = iv->pos;
  prepos = pos;

  int32_t npm = 0;
  for (int32_t k = 0; odr.read(iv); ++k) {
    if (bcf_get_info_int32(odr.hdr, iv, "AC", &info_ac, &n_ac) < 0) {continue;}
    if (info_ac[0] < 1) {continue;}

    if ( k % verbose == 0 )
      notice("Processing %d variants at %s:%d. Created %d PM variants", k, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1, npm);
    if (info_ac[0] >= maxac) {
      odw.write(iv);
      continue;
    }
      if (ivstack.size() > 1 && acsum < maxac && pos - prepos > winsize) { // Merge
        bcf1_t *nv = bcf_dup(ivstack.back());
        // bcf1_t *nv = bcf_dup(ivstack[0]);
        std::vector<int32_t> carry;
        std::vector<int32_t> config;
        int32_t newac = acstack.back();
        if ( bcf_get_genotypes(odr.hdr, nv, &p_gt, &n_gt) < 0 ) {
          error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr.hdr, nv->rid), nv->pos+1);
        }
        int32_t m_gt[2*nsamples];

        for (int32_t it = 0; it < nsamples; ++it) {
          if (bcf_gt_allele(p_gt[2*it]) > 0 || bcf_gt_allele(p_gt[2*it+1]) > 0) {
            carry.push_back(it);
          }
          m_gt[2*it] = p_gt[2*it]; m_gt[2*it+1] = p_gt[2*it+1];
        }
        config.push_back(acstack.back());

        for (uint32_t j = 0; j < ivstack.size()-1; ++j) {
          if ( bcf_get_genotypes(odr.hdr, ivstack[j], &p_gt, &n_gt) < 0 ) {
            error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr.hdr, ivstack[j]->rid), ivstack[j]->pos+1);
          }
          bool flag = 0;
          std::vector<int32_t> tmpcarry;
          for (int32_t it = 0; it < nsamples; ++it) {
            if (bcf_gt_allele(p_gt[2*it]) > 0 || bcf_gt_allele(p_gt[2*it+1]) > 0) {
              if (std::find(carry.begin(),carry.end(),it) != carry.end()) {
                flag = 1; break;
              }
              tmpcarry.push_back(it);
              m_gt[2*it] = p_gt[2*it]; m_gt[2*it+1] = p_gt[2*it+1];
            }
          }
          if (flag) { // Have overlapped carriers
            for (auto it : tmpcarry) {
              m_gt[2*it] = bcf_gt_phased(0); m_gt[2*it+1] = bcf_gt_phased(0);
            }
          } else {
            for (auto it : tmpcarry) {
              carry.push_back(it);
            }
            newac += acstack[j];
            config.push_back(acstack[j]);
          }
        }
        if (config.size() > 1) { // Output a new record
          nv->pos = iv->pos-1;
          float af = newac*0.5/nsamples;
          bcf_update_info_int32(odw.hdr, nv, "AC", &newac, 1);
          bcf_update_info_float(odw.hdr, nv, "AF", &af, 1);
          bcf_update_info_float(odw.hdr, nv, "MAF", &af, 1);
          bcf_update_genotypes(odw.hdr, nv, m_gt, nsamples*2);
          bcf_update_info_int32(odw.hdr, nv, "PMid", &carry[0], carry.size());
          bcf_update_info_int32(odw.hdr, nv, "PMconfig", &config[0], config.size());
          odw.write(nv);
          npm++;
          // New window
          prepos = pos;
          acstack.clear();
          for (auto & v : ivstack) {
            bcf_destroy(v);
            v = NULL;
          }
          ivstack.clear();
          acsum = 0;
        }
        bcf_destroy(nv);
      }
      else if (acsum > maxac || pos - prepos >= winsize) {
      prepos = pos;
      acstack.clear();
      for (auto & v : ivstack) {
        bcf_destroy(v);
        v = NULL;
      }
      ivstack.clear();
      acsum = 0;
      }
    ivstack.push_back(bcf_dup(iv));
    acstack.push_back(info_ac[0]);
    pos = iv->pos;
    acsum += info_ac[0];
    if (info_ac[0] > 1) {
      odw.write(iv);
    }
  }
  notice("Finished");
  odw.close();

  return 0;
}





