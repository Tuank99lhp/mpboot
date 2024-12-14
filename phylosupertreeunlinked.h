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

    void mapTrees();

    void computeParsimonyTree(const char *out_prefix, Alignment *alignment);

    bool isBifurcating(Node *node = NULL, Node *dad = NULL);

    string getTreeString();

    void readTreeString(const string &tree_string);

    void setRootNode(char *my_root);

    void computeBranchLengths();

    void printTree(ostream & out, int brtype = WT_BR_LEN);

    void printResultTree(string suffix = "");

    double treeLength(Node *node = NULL, Node *dad = NULL);

    double treeLengthInternal(double epsilon, Node *node = NULL, Node *dad = NULL);

    string doNNISearch(int &nniCount, int &nniSteps);

    double doTreeSearch();

    void summarizeBootstrap(Params &params);

    void writeUFBootTrees(Params &params, StrVector &removed_seqs, StrVector &twin_seqs);

};

#endif /* phylosupertreeunlinked_h */
