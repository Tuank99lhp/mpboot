#include "genetree.h"

GeneTree::GeneTree() : IQTree() {
    // Constructor implementation
}

GeneTree::GeneTree(Alignment *aln) : IQTree(aln) {
    // Constructor implementation
}

GeneTree::GeneTree(const string &treeString) : IQTree() {
    readTreeString(treeString);
    
    seedNode = (GeneNode*) root->neighbors[0]->node;
    assert(!seedNode->isLeaf());
}

GeneTree::~GeneTree() {
    // Destructor implementation
}

GeneNode *GeneTree::newNode(int node_id, const char* node_name) {
    return new GeneNode(node_id, node_name);
}
GeneNode *GeneTree::newNode(int node_id, int node_name) {
    return new GeneNode(node_id, node_name);
}

void GeneTree::setOriginalNodeIndex(const map<string, int> &seqNameToIndex) {
    totalLeafNum = seqNameToIndex.size();
    
    NodeVector taxa;
    getTaxa(taxa);

    for (auto taxon : taxa) {
        assert(taxon->isLeaf() && seqNameToIndex.find(taxon->name) != seqNameToIndex.end());
        ((GeneNode*) (taxon))->setOriginalIndex(seqNameToIndex.at(taxon->name));
    }
}

void GeneTree::addSplitOfNode(GeneNode *node) {
    assert(node->split);
    assert(splitToNode.find(*node->split) == splitToNode.end());
    splitToNode[*node->split] = node;
}

void GeneTree::removeSplitOfNode(GeneNode *node) {
    assert(node->split);
    assert(splitToNode[*node->split] == node);
    assert(splitToNode.erase(*node->split) == 1);
}

Split* GeneTree::getLeavesSplit() {
    assert(seedNode != NULL);
    if (splitToNode.empty() || seedNode->split == NULL) {
        assert(splitToNode.empty() && seedNode->split == NULL);
        encodeSplits();
    }
    return seedNode->split;
}

void GeneTree::dfsEncodeSplits(GeneNode *node, GeneNode *parent) {
    node->parent = parent;
    assert(node->getNumChildren() != 1);
    
    if (node->split) {
        delete node->split;
    }
    node->split = new Split(totalLeafNum, 0);

    if (node->isLeaf()) {
        node->split->addTaxon(((GeneNode*) node)->getOriginalIndex());
    }

    for (auto child: node->getChildren()) {
        dfsEncodeSplits(child, node);
        *node->split += *child->split;
    }

    addSplitOfNode(node);
}

void GeneTree::encodeSplits() {
    assert(seedNode != NULL && totalLeafNum > 0);

    if (seedNode->split != NULL) {
        return;
    }
    
    splitToNode.clear();

    if (seedNode->getNumChildren() == 2) {
        deroot();
    }

    dfsEncodeSplits((GeneNode*) seedNode, NULL);
}

void GeneTree::resetAfterMergeSCM() {
    delete seedNode->split;
    seedNode->split = NULL;
    splitToNode.clear();
}

void GeneTree::makeFirstChildOfRoot(GeneNode *node) {
    assert(node->parent != NULL);
    reroot(node->parent);
    assert(seedNode == node->parent);

    for (int i = 0 ; i < seedNode->neighbors.size(); ++i) {
        if (seedNode->neighbors[i]->node == node) {
            swap(seedNode->neighbors[i], seedNode->neighbors[0]);
            return;
        }
    }
    assert(0);
}

void GeneTree::reroot(GeneNode *node) {
    GeneNode *oldPar = node->parent;
    if (oldPar == NULL) {
        assert(seedNode == node);
        return;
    }

    if (oldPar != seedNode) {
        reroot(oldPar);
    }

    assert(oldPar == seedNode && oldPar->getNumChildren() > 2);

    removeSplitOfNode(oldPar);
    removeSplitOfNode(node);

    oldPar->removeChild(node);
    *oldPar->split -= *node->split;

    node->addChild(oldPar);
    *node->split += *oldPar->split;
        
    addSplitOfNode(oldPar);
    addSplitOfNode(node);

    seedNode = node;
}

void GeneTree::deroot() {
    if (!seedNode || seedNode->degree() != 2) {
        return;
    }

    cout << "Derooting tree" << endl;

    assert(rooted);

    GeneNode *child0 = (GeneNode*) seedNode->neighbors[0]->node;
    GeneNode *child1 = (GeneNode*) seedNode->neighbors[1]->node;

    if (child0->degree() < 3 && child1->degree() < 3) {
        cout << "Error tree's topology, both children have degree < 3" << endl;
        assert(0);
    }
    if (child0->degree() < 3) {
        swap(child0, child1);
    }

    collapseEdge(seedNode, child0);
}

void GeneTree::collapseEdge(GeneNode *head, GeneNode *tail) {
    assert(tail->isLeaf() == false);
    assert(head->parent == tail || tail->parent == head);
    
    bool down = (tail->parent == head);
    
    // Collapse should be from head = parent to tail = child in this version
    assert(down);

    if (down) {
        head->removeChild(tail);
    
        for (auto neighbor: tail->getChildren()) {
            tail->removeChild(neighbor);
            head->addChild(neighbor);
        }
    } else {
        tail->removeChild(head);
        GeneNode *parent = tail->parent;
        if (parent) {
            parent->removeChild(tail);
            parent->addChild(head);
        }
    }

    delete tail;
}

void GeneTree::collapseSubtree(GeneNode *node, bool do_collapse) {
    if (node->isLeaf()) {
        return;
    }

    for (auto child: node->getChildren()) {
        collapseSubtree(child, true);
    }

    if (do_collapse) {
        assert(node->parent != NULL);
        collapseEdge(node->parent, node);
    }
}

void GeneTree::readTreeString(const string &treeString) {
    stringstream str;
	str << treeString;
	str.seekg(0, ios::beg);
	freeNode();
	readTree(str, rooted);
}

void GeneTree::printResultTree(string fileName, bool isAppend) {
    printTree(fileName.c_str(), WT_SORT_TAXA | WT_NEWLINE | (isAppend ? WT_APPEND : 0));
}