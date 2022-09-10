#include "parsimony_tree.h"

int32_t parsimonyTree::add_variable_site_vcf(const bcf_hdr_t* odr.hdr, const bcf1_t* iv) {

	int32_t *p_gt, *info_ac, *info_an;
	p_gt = info_ac = info_an = NULL;
	int32_t n_gt = 0, n_ac = 0, n_an = 0;
	int32_t n_samples = bcf_hdr_nsamples(odr.hdr);
    // Allele info
    if (bcf_get_info_int32(odr.hdr, iv, "AC", &info_ac, &n_ac) < 0) {return 0;}
    if (bcf_get_info_int32(odr.hdr, iv, "AN", &info_an, &n_an) < 0) {return 0;}
    std::string ref(iv->d.allele[0]);
    std::string alt(iv->d.allele[1]);
    // extract genotype
    if ( bcf_get_genotypes(odr.hdr, iv, &p_gt, &n_gt) < 0 ) {
        error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);
    }
    variableSites.push_back(allele(iv->pos+1, ref, alt, info_ac[0], info_ac[0]));
    siteID.push_back(iv->pos+1);
    // TODO: now assume haploid and use only the first allele
    // TODO: need to handle missingness
    gt = new bool[n_samples];
    for (int32_t i = 0; i < n_samples; ++i) {
        int32_t g1 = p_gt[2*i];
    	gt[i] = bcf_gt_is_missing(g1);
    	if (! gt[i])  {
    		gt[i] = bcf_gt_allele(g1) > 0;
    	}
    }
    gtmat.push_back(gt);
}
