#ifndef GENENODE_H
#define GENENODE_H

#include "phylonode.h"
#include "split.h"

class GeneNode : public PhyloNode {
    friend class GeneTree;

public:
    GeneNode();

    GeneNode(int id);

    GeneNode(int id, int name);

    GeneNode(int id, const char *name);

    ~GeneNode();

    void setOriginalIndex(int index);

    int getOriginalIndex() const;

    void addNeighbor(Node *node, double length = 0, int id = -1) override;

    bool removeNeighbor(GeneNode *node);

    void addChild(GeneNode *child);

    void removeChild(GeneNode *child);

    vector<GeneNode*> getChildren();

    int getNumChildren();

    int originalIndex;

    GeneNode *parent = NULL;

    Split *split = NULL;
};

class GeneNeighbor : public PhyloNeighbor {
    friend class GeneNode;
    friend class GeneTree;

public:
    GeneNeighbor(GeneNode *tail) : PhyloNeighbor(tail, -1) {}

    ~GeneNeighbor() {}
};

#endif