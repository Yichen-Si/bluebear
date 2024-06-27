#include "cramore.h"
#include "bcf_filter_arg.h"
#include "bcf_ordered_reader.h"

#include "pbwt_build.h"

void pbwtSuffix(pbwtCursor& pc, std::vector<bool*>& gtmat, std::vector<int32_t>& posvec, std::vector<int32_t*>& dmat, std::vector<int32_t*>& rmat) {
	int32_t N = (int32_t) posvec.size();
	for (int32_t k = N-1; k >= 0; --k) {
		pc.ForwardsAD_suffix(gtmat[k], posvec[k]);
		memcpy(dmat[k], pc.d, pc.M*sizeof(int32_t));
		pc.ReverseA(rmat[k]);
	}
}

// Input bcf file
// Output suffix pbwt by chunk
int32_t pbwtBuildSuffix(int32_t argc, char** argv) {

	std::string inVcf, chrom, reg, out;
	int32_t verbose = 10000;
	int32_t chunksize = 100000;   // Read genotype by chunk
	int32_t max_store_interval = 100000; // bp
	int32_t store_interval = 10000;      // markers
	int32_t min_variant = 1;
	int32_t min_ac = 10;
	double max_missing = 0.05;
	bool haploid = false;
	bool snp_only = false;

	paramList pl;
	BEGIN_LONG_PARAMS(longParameters)
		LONG_PARAM_GROUP("Input Sites", NULL)
		LONG_STRING_PARAM("in-vcf",&inVcf, "Input VCF/BCF file")
		LONG_STRING_PARAM("region",&reg,"Genomic region to focus on")

		LONG_PARAM_GROUP("Additional Options", NULL)
		LONG_INT_PARAM("chunk-size",&chunksize, "Read, process and write by chunk")
		LONG_INT_PARAM("min-ac",&min_ac, "Minimum allele count to include")
		LONG_PARAM("haploid",&haploid, "Assume input is haploid coded in diploid vcf form and use only the first allele e.g. for chrY")
		LONG_PARAM("snp-only",&snp_only, "Only consider SNPs")

		LONG_PARAM_GROUP("Output Options", NULL)
		LONG_STRING_PARAM("out", &out, "Output file prefix")
		LONG_INT_PARAM("store-per-markers",&store_interval, "Snapshot per X markers (after filtering)")
		LONG_INT_PARAM("max-store-interval",&max_store_interval, "Snapshot at least per X bp")
		LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")

	END_LONG_PARAMS();
	pl.Add(new longParams("Available Options", longParameters));
	pl.Read(argc, argv);
	pl.Status();

	if ( inVcf.empty() || reg.empty() || out.empty() ) {
		error("[E:%s:%d %s] --in-vcf, --region, --out are required parameters",__FILE__,__LINE__,__FUNCTION__);
	}

	int32_t start = 0, end = 300000000;
	std::vector<std::string> v;
	split(v, ":-", reg);
	if (!str2int32(v[1], start) || !str2int32(v[2], end)) {
		error("Invalid region.");
	}
	chrom = v[0];

	std::vector<GenomeInterval> intervals;
	BCFOrderedReader odr(inVcf, intervals);
	bcf1_t* iv = bcf_init();
	int32_t nsamples = bcf_hdr_nsamples(odr.hdr);
	int32_t M = nsamples;
	if (!haploid) {
		M *= 2;
	}
	if (nsamples == 0) {
		error("Did not find any samples in the VCF file.");
	}
	notice("Identifying %d samples. Will store snapshot every %d markers or after %d bp", nsamples, store_interval, max_store_interval);

	int32_t nchunk = (end-start)/chunksize;
	int32_t st, ed;

	pbwtCursor pc(M, start + (nchunk+1) * chunksize + 1);

	// Output
	std::string d_outf = out + ".suffix.dmat";
	std::string a_outf = out + ".suffix.amat";
	std::ofstream d_wf, a_wf;
	d_wf.open(d_outf);
	a_wf.open(a_outf);

	int32_t k = 0, n_rec = 0, n_out = 0;
	int32_t *p_gt = NULL;
	int32_t  n_gt = 0;
	int32_t *info_ac = NULL, *info_an = NULL;
	int32_t  n_ac = 0, n_an = 0;
	int32_t last_stored_pos = -1;

	for (int32_t ck = nchunk; ck >= 0; --ck) {

		st = start + ck*chunksize + 1; // 1-based inclusive
		ed = start + ck*chunksize + chunksize; // 1-based inclusive
		if (st > end) continue;
		if (ed > end) ed = end;
		reg = chrom + ":" + std::to_string(st) + "-" + std::to_string(ed);

		parse_intervals(intervals, "", reg);
		if (!odr.jump_to_interval(intervals[0])) {
			notice("Failed to jump to interval %s.", reg.c_str());
			continue;
		}
		notice("Processing %s.", reg.c_str());

		std::vector<bool*> gtmat;
		std::vector<int32_t> positions;
		// read marker
		for (k=0; odr.read(iv); ++k) {
			if ( k % verbose == 0 ) {
				notice("Processed %d markers %s:%d; kept %d.", k, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1, n_rec);
			}
			if (snp_only && !bcf_is_snp(iv)) {continue;}
			if (bcf_get_info_int32(odr.hdr, iv, "AC", &info_ac, &n_ac) < 0) {continue;}
			if (bcf_get_info_int32(odr.hdr, iv, "AN", &info_an, &n_an) < 0) {continue;}
			if (info_ac[0] < min_ac || info_ac[0] > info_an[0] - min_ac) {continue;}
			if (bcf_get_genotypes(odr.hdr, iv, &p_gt, &n_gt) < 0) {
				error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);
			}

			int32_t major = (info_ac[0] > info_an[0]/2) ? 1 : 0;
			int32_t n_miss = 0;
			bool *y = new bool[M];
			for (int32_t i = 0; i < M; ++i) {
				int32_t gt = p_gt[i];
				if (haploid) {
					gt = p_gt[i*2];
				}
				if (bcf_gt_is_missing(gt)) {
					n_miss++;
					y[i] = major;
				} else {
					y[i] = ((bcf_gt_allele(gt) > 0) ? 1 : 0);
				}
			}
			if (n_miss > max_missing * M) {
				delete [] y;
				continue;
			}
			gtmat.push_back(y);
			positions.push_back(iv->pos+1);
		}

		int32_t N = positions.size();
		if ( N == 0 ) {
			notice("Observed 0 informative markers. Skipped %s.", reg.c_str());
			continue;
		} else {
			notice("Read %d markers, start extending pbwt for %s.", N, reg.c_str());
		}

		// Build pbwt (suffix) for this block
		for (k = N-1; k >= 0; --k) {
			if (n_rec == 0) {
				pc.Reset(M, positions[k]);
			}
			pc.ForwardsAD_suffix(gtmat[k], positions[k]);
			if (n_rec % store_interval == 0 || last_stored_pos - positions[k] > max_store_interval) {
				d_wf << chrom << '\t' << positions[k];
				for (int32_t j = 0; j < M; ++j)
					d_wf << '\t' << pc.d[j];
				d_wf << '\n';
				a_wf << chrom << '\t' << positions[k];
				for (int32_t j = 0; j < M; ++j)
					a_wf << '\t' << pc.a[j];
				a_wf << '\n';
				last_stored_pos = positions[k];
				n_out++;
			}
			n_rec++;
		}

		for (int32_t k = 0; k < N; ++k) {
			delete [] gtmat[k];
		}
		notice("Chunk %s finished. Kept %d markers, snapshot at %d positions so far.", reg.c_str(), n_rec, n_out);
	}

	delete [] p_gt;
	bcf_destroy(iv);
	odr.close();
	d_wf.close();
	a_wf.close();

	return 0;
}
