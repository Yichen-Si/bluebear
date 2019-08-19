#include "cramore.h"
#include "commands.h"
#include "utils.h"

KHASH_MAP_INIT_STR(vdict, bcf_idinfo_t)
typedef khash_t(vdict) vdict_t;

int32_t cmdVcfIBS0Pairwise(int32_t argc, char** argv);
int32_t cmdVcfIBS0Baseline(int32_t argc, char** argv);
int32_t cmdVcfIBS0Unconditional(int32_t argc, char** argv);
int32_t cmdVcfIBS0View(int32_t argc, char** argv);
int32_t cmdVcfIBS0Flank(int32_t argc, char** argv);
int32_t cmdVcfRareShare(int32_t argc, char** argv);

int32_t pbwtBuildSuffix(int32_t argc, char** argv);
int32_t pbwtBuildPrefix(int32_t argc, char** argv);

int32_t IBS0Phase(int32_t argc, char** argv);

int32_t test(int32_t argc, char** argv);

int32_t main(int32_t argc, char** argv) {
  commandList cl;

  BEGIN_LONG_COMMANDS(longCommandlines)
    LONG_COMMAND("vcf-ibs0-pairwise",&cmdVcfIBS0Pairwise, "Pairwise IBS0 and rare allele sharing from BCF/VCF")
    LONG_COMMAND("vcf-ibs0-baseline",&cmdVcfIBS0Baseline, "Pairwise IBS0 without rare alleles from BCF/VCF")
    LONG_COMMAND("vcf-ibs0-unconditional",&cmdVcfIBS0Unconditional, "Pairwise unconditional IBS0 from BCF/VCF")
    LONG_COMMAND("vcf-ibs0-view",&cmdVcfIBS0View, "Dichotomized IBS0 length distribution from BCF/VCF")
    LONG_COMMAND("vcf-ibs0-flank",&cmdVcfIBS0Flank, "Flanking IBS0 length (null distribution) from BCF/VCF")
    LONG_COMMAND("vcf-rare-share",&cmdVcfRareShare, "Rare variants sharing matrix from BCF/VCF")

    LONG_COMMAND("pbwt-suffix",&pbwtBuildSuffix, "Build pbwt backwards from BCF/VCF")
    LONG_COMMAND("pbwt-prefix",&pbwtBuildPrefix, "Build pbwt forwards from BCF/VCF")

    LONG_COMMAND("ibs0-phase",&IBS0Phase, "Find switch error in phased BCF/VCF")

    LONG_COMMAND("test",&test, "Test")

  END_LONG_COMMANDS();

  cl.Add(new longCommands("Available Commands", longCommandlines));

  if ( argc < 2 ) {
    fprintf(stderr, " Copyright (c) 2009-2017 by Hyun Min Kang and Adrian Tan\n");
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
