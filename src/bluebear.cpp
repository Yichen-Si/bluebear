#include "cramore.h"
#include "commands.h"
#include "utils.h"

KHASH_MAP_INIT_STR(vdict, bcf_idinfo_t)
typedef khash_t(vdict) vdict_t;

int32_t KmerCount(int32_t argc, char** argv);
int32_t cmdVcfAddContexte(int32_t argc, char** argv);
int32_t KmerSFS(int32_t argc, char** argv);
int32_t MultiAllelicSites(int32_t argc, char** argv);

int32_t AnnotateAge(int32_t argc, char** argv);
int32_t AnnotatePM(int32_t argc, char** argv);
int32_t AnnotateRefAC(int32_t argc, char** argv);
int32_t cmdVcfSummerRefAC(int32_t argc, char** argv);
int32_t AnnotateRareWithTotalAC(int32_t argc, char** argv);
int32_t cmdVcfSingletonRefAC(int32_t argc, char** argv);

int32_t cmdVcfSampleSummary(int32_t argc, char** argv);
int32_t cmdVcfCountSingleton(int32_t argc, char** argv);

// int32_t tmpUpdateCMInfo(int32_t argc, char** argv); // tmp
int32_t cmdVcfIBS0Pairwise(int32_t argc, char** argv); // zz
int32_t cmdVcfIBS0Baseline(int32_t argc, char** argv); // zz
int32_t cmdVcfIBS0Unconditional(int32_t argc, char** argv); // zz
int32_t cmdVcfIBS0View(int32_t argc, char** argv);
int32_t cmdVcfIBS0Flank(int32_t argc, char** argv); // zz
int32_t cmdVcfIBS0full(int32_t argc, char** argv); // zz
int32_t RandomPairIBS0(int32_t argc, char** argv);
int32_t VCFInsertPM(int32_t argc, char** argv);
int32_t CDFInfoInt(int32_t argc, char** argv);
int32_t CDFInfoFloatByKmer(int32_t argc, char** argv);
int32_t CDFInfoFloatByAnn(int32_t argc, char** argv); // tmp
int32_t CDFInfoFloatSimu(int32_t argc, char** argv); // tmp

int32_t IBS0PairwiseScan(int32_t argc, char** argv);
int32_t AnnotateIBS0AroundRare_Small(int32_t argc, char** argv);
int32_t AnnotateIBDAroundRare(int32_t argc, char** argv);

int32_t cmdIBS0Bridge(int32_t argc, char** argv);
int32_t mapPt2Interval(int32_t argc, char** argv);

int32_t cmdVcfRareShare(int32_t argc, char** argv);
int32_t cmdVcfSFS(int32_t argc, char** argv);

int32_t pbwtBuildSuffix(int32_t argc, char** argv);
int32_t pbwtBuildPrefix(int32_t argc, char** argv);
int32_t hapIBDpbwtLeft(int32_t argc, char** argv);
int32_t hapIBDpbwtRight(int32_t argc, char** argv);

int32_t HapibdVSnoibs0(int32_t argc, char** argv);

// int32_t IBS0Phase(int32_t argc, char** argv);
int32_t trioSwitchDetect(int32_t argc, char** argv);
int32_t trioSwitchDetect_onepass(int32_t argc, char** argv);
int32_t IBS0PhaseForward(int32_t argc, char** argv);
int32_t IBS0PhaseBackward(int32_t argc, char** argv);
int32_t RareIBS0PhaseForward(int32_t argc, char** argv);
int32_t RareIBS0PhaseBackward(int32_t argc, char** argv);
int32_t RareOnlyIBS0PhaseForward(int32_t argc, char** argv);
int32_t RareOnlyIBS0PhaseBackward(int32_t argc, char** argv);

int32_t RareBlipOnlyForward(int32_t argc, char** argv);
int32_t RareBlipOnlyStats(int32_t argc, char** argv);

int32_t IBS0AddPMDistr(int32_t argc, char** argv);
int32_t IBS0AddNoncarrierControl(int32_t argc, char** argv);
int32_t IBS0AddOutgroupControl(int32_t argc, char** argv);
int32_t IBS0Triallelic(int32_t argc, char** argv);

int32_t cmdDgeShuffle(int32_t argc, char** argv);

int32_t test(int32_t argc, char** argv);

int32_t main(int32_t argc, char** argv) {
  commandList cl;

  BEGIN_LONG_COMMANDS(longCommandlines)
    LONG_COMMAND("kmer-count",&KmerCount, "Count total occurence of kmers from fasta file (by chr)")
    LONG_COMMAND("kmer-anno",&cmdVcfAddContexte, "Add kmer context and CpG info to VCF/BCF")

    LONG_COMMAND("folded-sfs",&cmdVcfSFS, "Folded SFS from BCF/VCF")
    LONG_COMMAND("ctx-sfs",&KmerSFS, "Kmer context specific SFS from BCF/VCF")
    LONG_COMMAND("triallelic",&MultiAllelicSites, "Find tri-allelic sites from BCF/VCF")

    LONG_COMMAND("anno-age",&AnnotateAge, "Annotate variant age")
    LONG_COMMAND("anno-pm",&AnnotatePM, "Annotate posterior probability of being a parallel mutation")
    LONG_COMMAND("anno-total-carrier",&AnnotateRareWithTotalAC, "Annotate a query dataset with the ultra-rare variant AC and carriers in a reference dataset")
    LONG_COMMAND("anno-ref-ac",&AnnotateRefAC, "Annotate a query dataset with the AC in a reference dataset (optional: annotate carrier ancestry info)")
    LONG_COMMAND("singleton-ref-ac",&cmdVcfSingletonRefAC, "Summarize the AC in a large ref. pop of singletons in a smaller sample")
    LONG_COMMAND("rare-ref-ac",&cmdVcfSummerRefAC, "Summarize the AC in a large ref. pop of (rare) variants in a smaller sample")

    LONG_COMMAND("vcf-sample-summary",&cmdVcfSampleSummary, "Sample-level summary from BCF/VCF")
    LONG_COMMAND("vcf-count-singleton",&cmdVcfCountSingleton, "Count & list singletons from BCF/VCF")

    // Temperary ibs0 tests
    // LONG_COMMAND("vcf-ibs0-correct-cm",&tmpUpdateCMInfo, "Update genetic distance bug")
    LONG_COMMAND("vcf-ibs0-pairwise",&cmdVcfIBS0Pairwise, "Pairwise IBS0 and rare allele sharing from BCF/VCF")
    LONG_COMMAND("vcf-ibs0-baseline",&cmdVcfIBS0Baseline, "Pairwise IBS0 without rare alleles from BCF/VCF")
    LONG_COMMAND("vcf-ibs0-unconditional",&cmdVcfIBS0Unconditional, "Pairwise unconditional IBS0 from BCF/VCF")
    LONG_COMMAND("vcf-ibs0-view",&cmdVcfIBS0View, "Dichotomized IBS0 length distribution from BCF/VCF")
    LONG_COMMAND("vcf-ibs0-flank",&cmdVcfIBS0Flank, "Flanking IBS0 length (null distribution) from BCF/VCF")
    LONG_COMMAND("vcf-ibs0",&cmdVcfIBS0full, "Unconditional no-IBS0 length, full output from BCF/VCF")
    LONG_COMMAND("ibs0-scan",&IBS0PairwiseScan, "Pairwise IBS0 and rare allele sharing from BCF/VCF")
    LONG_COMMAND("vcf-ibs0-pos",&RandomPairIBS0, "Given a list of positions, randomly sample a pair of individuals for each, and get flanking ibs0 loci from BCF/VCF")
    LONG_COMMAND("insert-pm",&VCFInsertPM, "Merge nearby rare variants to create artificial parallel mutation from & to BCF/VCF")
    LONG_COMMAND("info-int-cdf",&CDFInfoInt, "Summarize empirical CDF from VCF INFO")
    LONG_COMMAND("info-flt-cdf",&CDFInfoFloatByKmer, "Summarize empirical CDF from VCF INFO (by kmer)")
    LONG_COMMAND("info-flt-cdf-ann",&CDFInfoFloatByAnn, "Summarize empirical CDF from VCF INFO (by functional annotations)")
    LONG_COMMAND("info-flt-cdf-simu",&CDFInfoFloatSimu, "Summarize empirical CDF from VCF INFO")


    LONG_COMMAND("anno-ibd",&AnnotateIBDAroundRare, "Annotate pairwise IBD among rare allele carriers from input BCF/VCF")
    LONG_COMMAND("anno-ibs0-small",&AnnotateIBS0AroundRare_Small, "Annotate pairwise IBS0 among rare allele carriers from input BCF/VCF")

    LONG_COMMAND("ibs0-bridge",&cmdIBS0Bridge, "Connecting unconditional boundary output for long-range no-IBS0")

    LONG_COMMAND("pt-to-interval",&mapPt2Interval, "Map points to intervals, both indexed by the carriers.")


    LONG_COMMAND("vcf-rare-share",&cmdVcfRareShare, "Rare variants sharing matrix from BCF/VCF")

    LONG_COMMAND("pbwt-suffix",&pbwtBuildSuffix, "Build pbwt backwards from BCF/VCF")
    LONG_COMMAND("pbwt-prefix",&pbwtBuildPrefix, "Build pbwt forwards from BCF/VCF")

    LONG_COMMAND("hap-ibd-left",&hapIBDpbwtLeft, "Haplotype matching length preceeding a shared rare allele from BCF/VCF")
    LONG_COMMAND("hap-ibd-right",&hapIBDpbwtRight, "Haplotype matching length succeeding a shared rare allele from BCF/VCF")

    LONG_COMMAND("ibd-ibs",&HapibdVSnoibs0, "Haplotype & ibs0 based IBD around given positions from BCF/VCF")

    // LONG_COMMAND("ibs0-phase",&IBS0Phase, "Find switch error in phased BCF/VCF")
    LONG_COMMAND("ibs0-phase-forward",&IBS0PhaseForward, "Starting from an arbitrary position, proceed froward to find switch error in phased BCF/VCF")
    LONG_COMMAND("ibs0-phase-backward",&IBS0PhaseBackward, "Starting from an arbitrary position, proceed backward to find switch error in phased BCF/VCF")

    LONG_COMMAND("rare-ibs0-phase-forward",&RareIBS0PhaseForward, "Starting from an arbitrary position, proceed froward to find switch error in phased BCF/VCF")
    LONG_COMMAND("rare-ibs0-phase-backward",&RareIBS0PhaseBackward, "Starting from an arbitrary position, proceed backward to find switch error in phased BCF/VCF")

    LONG_COMMAND("rare-phase-forward",&RareOnlyIBS0PhaseForward, "Based on shared rare allele, proceed froward to find switch error in phased BCF/VCF")
    LONG_COMMAND("rare-phase-backward",&RareOnlyIBS0PhaseBackward, "Based on shared rare allele, proceed backward to find switch error in phased BCF/VCF")

    LONG_COMMAND("rare-blip",&RareBlipOnlyForward, "Based on shared rare allele and haplotype matching, proceed froward to detect blip in phased BCF/VCF")
    LONG_COMMAND("rare-blip-stats",&RareBlipOnlyStats, "Statistics of ibd around rare variants for detecting blip in phased BCF/VCF")

    LONG_COMMAND("trio-switch",&trioSwitchDetect, "Detect switch errors in trio-child")
    LONG_COMMAND("trio-switch-full",&trioSwitchDetect_onepass, "Detect switch errors & blips in trio-child")

LONG_COMMAND("btw-carrier-sets",&IBS0AddPMDistr, "Reference distribution for parallel mutations")
LONG_COMMAND("add-null-by-site",&IBS0AddNoncarrierControl, "Add no-IBS0 length between non-carriers. Currently only for doubletons")
LONG_COMMAND("sample-outgroup-by-site",&IBS0AddOutgroupControl, "Add average no-IBS0 length between a randomly sampled non-carriers and the carriers. Currently only tested for doubletons")
LONG_COMMAND("ibs0-tri-allelic",&IBS0Triallelic, "No-IBS0 between carriers of different alleles at tri-allelic sites")

LONG_COMMAND("dge-shuffle",&cmdDgeShuffle, "Shuffle digital expression matrix")

    LONG_COMMAND("test",&test, "Test")

  END_LONG_COMMANDS();

  cl.Add(new longCommands("Available Commands", longCommandlines));

  if ( argc < 2 ) {
    // fprintf(stderr, " Copyright (c) 2009-2017 by Hyun Min Kang and Adrian Tan\n");
    fprintf(stderr, " Licensed under the Apache License v2.0 http://www.apache.org/licenses/\n\n");
    fprintf(stderr, "To run a specific command      : %s [command] [options]\n",argv[0]);
    fprintf(stderr, "For detailed instructions, run : %s --help\n",argv[0]);
    cl.Status();
    return 1;
  }
  else {
    if ( strcmp(argv[1],"--help") == 0 ) {
      cl.HelpMessage();
    }
    else {
      return cl.Read(argc, argv);
    }
  }
  return 0;
}
