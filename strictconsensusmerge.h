#ifndef STRICTCONSENSUSMERGE_H
#define STRICTCONSENSUSMERGE_H

#include "genetree.h"

class StrictConsensusMerge {
public:
    StrictConsensusMerge(const StrVector &sourceTrees, const map<string, int> &seqNameToIndex);

    ~StrictConsensusMerge();

    GeneTree *getSCMTree();

private:
    void rerootOnLowestCommonIndexPath(GeneTree *tree, Split commonMask);

    void collapsePathNotFound(
        map<Split, vector<pair<Split, GeneNode*>>> &treeRelevantSplits, 
        map<Split, vector<pair<Split, GeneNode*>>> &otherTreeRelevantSplits, 
        GeneTree *tree
    );

    void getNextPair(int &tree1Index, int &tree2Index);

    void pairwiseMerger(GeneTree *tree1, GeneTree *tree2);

    int getOverlap(GeneTree *tree1, GeneTree *tree2);

    vector<GeneTree*> trees;
};

#endif