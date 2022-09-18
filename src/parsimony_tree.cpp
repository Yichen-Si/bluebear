#include "parsimony_tree.h"

int32_t parsimonyTree::add_allele(int32_t _s, std::string _a, std::string _b, int32_t _x, int32_t _y) {
    if (_s == last_site_id) {
        variableSites.back().add_allele(_a, _b, _x, _y);
        return 0;
    } else {
        variableSites.push_back(variableSite(_S, _a, _b, _x, _y));
        return 1;
    }
}

int32_t parsimonyTree::add_variable_site_vcf_hyploid(const bcf_hdr_t* odr.hdr, const bcf1_t* iv) {

	int32_t *p_gt, *info_ac, *info_an;
	p_gt = info_ac = info_an = NULL;
	int32_t n_gt = 0, n_ac = 0, n_an = 0;
    // Allele info
    if (bcf_get_info_int32(odr.hdr, iv, "AC", &info_ac, &n_ac) < 0) {return 0;}
    if (bcf_get_info_int32(odr.hdr, iv, "AN", &info_an, &n_an) < 0) {return 0;}
    std::string ref(iv->d.allele[0]);
    std::string alt(iv->d.allele[1]);
    // extract genotype
    if ( bcf_get_genotypes(odr.hdr, iv, &p_gt, &n_gt) < 0 ) {
        error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);
    }
    if (last_site_id == iv->pos+1) {

    }
    add_allele(iv->pos+1, ref, alt, info_ac[0], info_ac[0]);
    // TODO: now assume haploid and use only the first allele
    // TODO: need to handle missingness
    gt = new bool[nleaves];
    for (int32_t i = 0; i < nleaves; ++i) {
        int32_t g1 = p_gt[2*i];
        if (bcf_gt_is_missing(g1)) {
            gt[i] = 0;
        } else {
    		gt[i] = ((bcf_gt_allele(g1) > 0) ? 1 : 0);
    	}
    }
    gtmat.push_back(gt);
}
