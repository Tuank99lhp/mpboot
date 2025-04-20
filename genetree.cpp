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

    // Collapse all node has degree 2, except the seedNode
    dfsFixTree(seedNode, NULL);
}

void GeneTree::dfsFixTree(GeneNode *node, GeneNode *parent) {
    node->parent = parent;
    
    if (node->isLeaf()) {
        return;
    }

    for (auto child: node->getChildren()) {
        dfsFixTree(child, node);
    }

    if (node->parent != NULL && node->degree() == 2) {
        collapseEdge(node->parent, node);
    }
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

    // cout << "Derooting tree" << endl;

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

void GeneTree::reInitializeTree(GeneNode *node, GeneNode* parent) {
    if (!node) {
        setRootLeaf(NULL);
        node = (GeneNode*) root;
        nodeNum = getNumTaxa();
        leafNum = 0;
        branchNum = 0;
    }

    node->id = (node->isLeaf() ? leafNum++ : nodeNum++);

    FOR_NEIGHBOR_IT(node, parent, it) {
        (*it)->id = branchNum;
        (*it)->node->findNeighbor(node)->id = branchNum;
        branchNum++;
        reInitializeTree((GeneNode*)(*it)->node, node);
    }
}

GeneNode* GeneTree::getAnyOtherLeaf(GeneNode *node, GeneNode *parent) {
    if (!node) {
        assert(root);
        node = (GeneNode*) root;
    } else if (node->isLeaf()) {
        return (node->name != root->name) ? node : NULL;
    }

    FOR_NEIGHBOR_IT(node, parent, it) {
        GeneNode *child = (GeneNode*) (*it)->node;
        GeneNode *leaf = getAnyOtherLeaf(child, node);
        if (leaf) {
            return leaf;
        }
    }
    return NULL;
}

void GeneTree::setRootLeaf(char *my_root) {
    string root_name;
    if (my_root) {
        root_name = my_root;
    } else {
        if (root && root->isLeaf()) {
            return;
        }
        if (aln) {
            root_name = aln->getSeqName(0);
        } else {
            root = getAnyOtherLeaf();
            assert(root && root->isLeaf());
            return;
        }
    }
    root = findLeafName(root_name);
    assert(root);
}

void GeneTree::getPolytomies(vector<GeneNode*> &polytomies, GeneNode *node, GeneNode *parent) {
    if (!node) {
        node = (GeneNode*) root;
    }

    if (node->degree() > 3) {
        polytomies.push_back(node);
    }

    FOR_NEIGHBOR_IT(node, parent, it) {
        GeneNode *child = (GeneNode*) (*it)->node;
        getPolytomies(polytomies, child, node);
    }
}

void GeneTree::getRelabelMap(map<string, string> &relabel, map<string, GeneNode*> &delabel, int &label, GeneNode *node, GeneNode *parent) {
    node->parent = parent;

    if (node->isLeaf()) {
        relabel[node->name] = to_string(label);
        assert(parent != NULL);
    }

    for (auto child: node->getChildren()) {
        if (parent == NULL) {
            delabel[to_string(++label)] = child;
        }
        getRelabelMap(relabel, delabel, label, child, node);
    }
}

void GeneTree::relabelAndCollapse(const map<string, string> &relabel) {
    setRootLeaf(NULL);
    relabelTree(relabel);
    deleteDuplicateNode();
    GeneNode *newRoot = getAnyOtherLeaf();
    if (newRoot != NULL) {
        assert(root->name != newRoot->name);
        root = newRoot;
        deleteDuplicateNode();
    } else {
        assert(getNumTaxa() <= 2);
    }
    reInitializeTree();
    if (leafNum > 2) {
        map<string, int> checkUniqueLabel;
        doubleCheckUniqueName(checkUniqueLabel);
    }
}

void GeneTree::doubleCheckUniqueName(map<string, int> &checkUniqueLabel, GeneNode *node, GeneNode *parent) {
    if (!node) {
        node = (GeneNode*) root;
    }

    if (node->isLeaf()) {
        assert(checkUniqueLabel.find(node->name) == checkUniqueLabel.end());
        checkUniqueLabel[node->name] = 1;
    }

    FOR_NEIGHBOR_IT(node, parent, it) {
        GeneNode *child = (GeneNode*) (*it)->node;
        doubleCheckUniqueName(checkUniqueLabel, child, node);
    }
}

void GeneTree::relabelTree(const map<string, string> &relabel, GeneNode *node, GeneNode *parent) {
    if (!node) {
        node = (GeneNode*) root;
    }

    if (node->isLeaf()) {
        assert(relabel.find(node->name) != relabel.end());
        node->name = relabel.at(node->name);
    }

    FOR_NEIGHBOR_IT(node, parent, it) {
        GeneNode *child = (GeneNode*) (*it)->node;
        relabelTree(relabel, child, node);
    }
}

void GeneTree::delabelTree(map<string, GeneNode*> &delabel, GeneNode *node, GeneNode *parent) {
    if (!node) {
        assert(root && root->isLeaf());
        node = (GeneNode*) root->neighbors[0]->node;
    }

    node->parent = parent;

    if (node->isLeaf()) {
        assert(delabel.find(node->name) != delabel.end());
        GeneNode *subtree = delabel[node->name];
        delabel.erase(node->name);

        assert(subtree->parent != NULL && parent != NULL);
        
        subtree->parent->removeChild(subtree);
        parent->removeChild(node);
        parent->addChild(subtree);
        
        delete node;
        return;
    }
    
    for (auto child: node->getChildren()) {
        delabelTree(delabel, child, node);
    }

    if (parent == NULL && delabel.empty() == false) {
        for (auto [_, subtree]: delabel) {
            assert(subtree->parent != NULL);
            subtree->parent->removeChild(subtree);
            node->addChild(subtree);
        }
        delabel.clear();
    }
}

void GeneTree::deleteDuplicateNode(GeneNode *node, GeneNode *parent) {
    if (!node) {
        node = (GeneNode*) root;
    } else if (node->isLeaf()) {
        node->parent = parent;
        return;
    }
    
    node->parent = parent;

    for (auto child: node->getChildren()) {
        deleteDuplicateNode(child, node);
    }

    bool doDelete = true;

    GeneNode *prevChild = NULL;

    for (auto child: node->getChildren()) {
        if (child->isLeaf() == false) {
            doDelete = false;
            break;
        }

        if (prevChild == NULL) {
            prevChild = child;
        } else if (child->name != prevChild->name) {
            doDelete = false;
            break;
        }
    }

    if (doDelete && parent != NULL) {
        parent->removeChild(node);
        node->removeChild(prevChild);
        parent->addChild(prevChild);

        for (auto child: node->getChildren()) {
            delete child;
        }
        delete node;
    }
}