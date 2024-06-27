#include "utils.h"
#include "cramore.h"
#include "hts_utils.h"
#include <chrono>
#include <random>
#include <numeric>
#include "bcf_ordered_reader.h"

int32_t temp(int32_t argc, char** argv) {

    // std::string inBin, out, outf;
    // int32_t verbose = 10000;
    // int32_t debug = 0;
    // int32_t acmin = 2;

    // paramList pl;
    // BEGIN_LONG_PARAMS(longParameters)
    //     LONG_PARAM_GROUP("Input Sites", NULL)
    //     LONG_STRING_PARAM("in-bin",&inBin, "Input binary file")
    //     LONG_PARAM_GROUP("Output Options", NULL)
    //     LONG_STRING_PARAM("out", &out, "Output file prefix")
    //     LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")
    //     LONG_INT_PARAM("debug",&debug,"Debug")
    //     LONG_INT_PARAM("min-ac",&acmin,"Minimum allele count")
    // END_LONG_PARAMS();
    // pl.Add(new longParams("Available Options", longParameters));
    // pl.Read(argc, argv);
    // pl.Status();

    // std::vector<BinaryVariant> mixStorage;
    // parseBinaryFile(inBin, mixStorage);
    // int32_t nValues = mixStorage.size();
    // int32_t nalleles = mixStorage[0].nBits;
    // notice("Loaded %d variants", nValues);

    // outf = out + ".dprime.pairs.tsv.gz";
    // htsFile* wf = hts_open(outf.c_str(), "wz");
    // if (wf == NULL) {
    //     error("Could not open file: %s", outf.c_str());
    // }

    // int32_t npairs = 0, nrec = 0;
    // for (int32_t i = 1; i < nValues; ++i) {
    //     if (mixStorage[i].mOnes < acmin) {
    //         npairs += i;
    //         continue;
    //     }
    //     int32_t n1x = mixStorage[i].mOnes;
    //     for (int32_t j = 0; j < i; ++j) {
    //         npairs++;
    //         if (mixStorage[j].mOnes < acmin) {
    //             continue;
    //         }
    //         int32_t nx1 = mixStorage[j].mOnes;
    //         int32_t n11 = binary_intersect(mixStorage[i], mixStorage[j]);
    //         double dp = dprime(nalleles, n1x, nx1, n11);
    //         hprintf(wf, "%d\t%d\t%d\t%d\t%d\t%.6f\n", i, j, n1x, nx1, n11, dp);
    //         nrec++;
    //         if (nrec % verbose == 0) {
    //             notice("Processed %d pairs, recorded %d", npairs, nrec);
    //         }
    //     }
    //     if (debug && nrec > debug) {
    //         break;
    //     }
    // }

    // notice("Processed %d pairs, recorded %d", npairs, nrec);
    // hts_close(wf);

    return 0;
}
