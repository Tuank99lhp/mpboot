#ifndef GENETREE_H
#define GENETREE_H

#include "iqtree.h"
#include "genenode.h"

class GeneTree: public IQTree { 
public:
    GeneTree();

    GeneTree(Alignment *aln);

    GeneTree(const string &treeString);

    ~GeneTree();

    GeneNode *newNode(int node_id = -1, const char* node_name = NULL) override;

    GeneNode *newNode(int node_id, int node_name) override;

    void setOriginalNodeIndex(const map<string, int> &seqNameToIndex);

    void addSplitOfNode(GeneNode *node);

    void removeSplitOfNode(GeneNode *node);

    Split* getLeavesSplit();

    void dfsEncodeSplits(GeneNode *node, GeneNode *parent);

    void encodeSplits();

    void resetAfterMergeSCM();

    void makeFirstChildOfRoot(GeneNode *node);

    void reroot(GeneNode *node);

    void deroot();

    void collapseEdge(GeneNode *head, GeneNode *tail);

    void collapseSubtree(GeneNode *node, bool do_collapse);
    
    void readTreeString(const string &treeString);

    void printResultTree(string fileName, bool isAppend);

    Params treeParams;

    int totalLeafNum = 0;

    map<Split, GeneNode*> splitToNode;

    // pseudo root for unrooted tree
    GeneNode *seedNode = NULL;
};

#endif