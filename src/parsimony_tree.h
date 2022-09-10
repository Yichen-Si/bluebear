#ifndef PARSIMONY_TREE_H
#define PARSIMONY_TREE_H

#include "bcf_ordered_reader.h"
#include "pbwt_build.h"

struct parsimonyNode
{
    int32_t size, seqID, birthIter, nChild;
    parsimonyNode* parent;
    std::vector<parsimonyNode*> children;
};

struct allele
{
    int32_t siteID;
    std::string ref, alt;
    int32_t sampleAC, currentAC;
    // double score;
    allele(int32_t _s, std::string _a, std::string _b, int32_t _x, int32_t _y) : siteID(_s), ref(_a), alt(_b), sampleAC(_x), currentAC(_y) {};
};

class parsimonyTree
{
    private:

    public:

    std::vector<parsimonyNode*> activeNodes;
    std::vector<parsimonyNode*> pastNodes;
    std::vector<int8_t*> gtmat; // temporary? M x N
    std::vector<allele> variableSites;
    std::map<int32_t, std::vector<allele> > activeAlleles;

    parsimonyTree() {};
    ~parsimonyTree() {};

    int32_t add_variable_site_vcf(bcf1_t* iv) {};
    int32_t read_sequences() {};
    int32_t build_pbwt() {};
    int32_t run_greedy_parsimony() {};
    int32_t clearn() {};
    int32_t score_recurrent_mutation() {};
    int32_t remove_one_recurrent_mutation() {};
    int32_t clearn_rebuild() {};

};

#endif
