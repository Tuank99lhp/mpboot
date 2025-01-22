//
//  phylosupertreeunlinked.h
//  tree
//
//  Created by Minh Bui on 2/5/18.
//

#ifndef phylosupertreeunlinked_h
#define phylosupertreeunlinked_h

#include "phylosupertree.h"

/**
    Super-tree with separate partition tree topologies
 */
class PhyloSuperTreeUnlinked : public PhyloSuperTree {
public:
    /**
     constructors
     */
    PhyloSuperTreeUnlinked(Params &params);

    bool isSuperTreeUnlinked() {
        return true;
    }

    void readTree(istream &in, bool &is_rooted);

    void setAlignment(Alignment *alignment);

    void setParams(Params& params);

    void deleteAllPartialLh();

    void initializePLL(Params &params);

    void mapTrees();

    void initializeAllPartialPars();

    void clearAllPartialLH();

    int computeParsimony();

    double getBestScore();

    void printBestScores(int numBestScore);

    void computeParsimonyTree(const char *out_prefix, Alignment *alignment);

    string getTreeString();

    void readTreeString(const string &tree_string);

    void setRootNode(char *my_root);

    void computeBranchLengths();

    void printTree(ostream & out, int brtype = WT_BR_LEN);

    void printResultTree(string suffix = "");

    string doNNISearch(int &nniCount, int &nniSteps);

    double doTreeSearch();

    void summarizeBootstrap(Params &params);

    void writeUFBootTrees(Params &params, StrVector &removed_seqs, StrVector &twin_seqs);

    void optimizeAlignments(Params &params);

    void computeInitialPLLTrees(Params &params);

    void initCandidateUnlinkedTreeSet(Params &params, int numInitTrees);

    void reportUnlinkedPhyloAnalysis(Params &params, string &original_model, vector<ModelInfo> &model_info, StrVector &removed_seqs, StrVector &twin_seqs);

};

#endif /* phylosupertreeunlinked_h */
