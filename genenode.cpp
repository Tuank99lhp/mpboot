#include "genenode.h"

GeneNode::GeneNode() : PhyloNode() {
    // Constructor implementation
}

GeneNode::GeneNode(int id) : PhyloNode(id) {
    // Constructor implementation
}

GeneNode::GeneNode(int id, int name) : PhyloNode(id, name) {
    // Constructor implementation
}

GeneNode::GeneNode(int id, const char *name) : PhyloNode(id, name) {
    // Constructor implementation
}

GeneNode::~GeneNode() {
    // Destructor implementation
}

void GeneNode::setOriginalIndex(int index) {
    originalIndex = index;
}

int GeneNode::getOriginalIndex() const {
    return originalIndex;
}

void GeneNode::addNeighbor(Node *node, double length, int id) {
    neighbors.push_back(new GeneNeighbor((GeneNode*)node));
}

bool GeneNode::removeNeighbor(GeneNode *node) {
    for (auto it = neighbors.begin(); it != neighbors.end(); ++it) {
        if ((*it)->node == node) {
            delete *it;
            neighbors.erase(it);
            return true;
        }
    }
    assert(0);
}

void GeneNode::addChild(GeneNode *child) {
    addNeighbor(child);
    child->addNeighbor(this);

    assert(child->parent == NULL);
    child->parent = this;
}

void GeneNode::removeChild(GeneNode *child) {
    removeNeighbor(child);
    child->removeNeighbor(this);
    
    assert(child->parent == this);
    child->parent = NULL;
}

vector<GeneNode*> GeneNode::getChildren() {
    vector<GeneNode*> children;
    FOR_NEIGHBOR_IT(this, parent, it) {
        children.push_back((GeneNode *)(*it)->node);
    }
    return children;
}

int GeneNode::getNumChildren() {
    return degree() - (parent ? 1 : 0);
}