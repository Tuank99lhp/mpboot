#include "phylosupertreeunlinked.h"

PhyloSuperTreeUnlinked::PhyloSuperTreeUnlinked(Params &params): PhyloSuperTree(params, true) {
    this->params = &params;

    for (auto it = begin(); it != end(); it++) {
        GeneTree* tree = (GeneTree*)(*it);
        tree->treeParams = *(this->params);
    }

    if (params.aln_file) {
        conAln = new Alignment(params.aln_file, params.sequence_type, params.intype);
        conAln->checkGappySeq();
    } else {
        // Temporarily not supporting
        assert(false);
    }
}

PhyloSuperTreeUnlinked::~PhyloSuperTreeUnlinked() {
    if (conAln) {
        delete conAln;
    }
    if (mrpAln) {
        delete mrpAln;
    }
    if (mrpTree) {
        delete mrpTree;
    }
    if (scmTree) {
        delete scmTree;
    }
}

StrVector PhyloSuperTreeUnlinked::getAllSeqNames() {
    if (allSeqNames.size() > 0) {
        return allSeqNames;
    }

    if (conAln) {
        for (int i = 0; i < conAln->getNSeq(); ++i) {
            string seqName = conAln->getSeqName(i);
            allSeqNames.push_back(seqName);
            seqNameToIndex[seqName] = i;
        }
    } else {
        for (auto it = begin(); it != end(); it++) {
            GeneTree* tree = (GeneTree*)(*it);
            for (int i = 0; i < tree->aln->getNSeq(); ++i) {
                string seqName = tree->aln->getSeqName(i);
                if (seqNameToIndex.find(seqName) == seqNameToIndex.end()) {
                    allSeqNames.push_back(seqName);
                    seqNameToIndex[seqName] = allSeqNames.size() - 1;
                }
            }
        }
    }

    return allSeqNames;
}

void PhyloSuperTreeUnlinked::runGeneTreesReconstruction() {
    for (auto it = begin(); it != end(); it++) {
        GeneTree* tree = (GeneTree*)(*it);
        runOptimizeAndReconstruction(tree->treeParams, tree);
        // TODO: get best tree or greedy consensus tree or random tree
    }
}

void PhyloSuperTreeUnlinked::printGeneTrees() {
    string treeFile(this->params->out_prefix);
    treeFile += ".genetreesfile";
    // open treeFile and remove all
    ofstream outFile(treeFile.c_str());
    outFile.close();
    
    for (auto it = begin(); it != end(); it++) {
        GeneTree* tree = (GeneTree*)(*it);
        tree->setRootLeaf(NULL);
        tree->printResultTree(treeFile, true);
    }
}

void PhyloSuperTreeUnlinked::dfsMRP(Node* u, Node* pa, int &time, vector<pair<int, int>> &eulerInternalBranch, map<string, int> &leafIndex) {
    int in = ++time;
    FOR_NEIGHBOR_IT(u, pa, it){
        Node* v = (*it)->node;
        dfsMRP(v, u, time, eulerInternalBranch, leafIndex);
    }
    int out = time;

    if (u->isLeaf()) {
        assert(in == out);
        leafIndex[u->name] = in;
    } else if (!pa->isLeaf()) {
        eulerInternalBranch.push_back(make_pair(in, out));
    }
}

void PhyloSuperTreeUnlinked::buildMRPMatrix() {
    StrVector seqNames = getAllSeqNames();
    StrVector sequences(seqNames.size());
    
    for (auto it = begin(); it != end(); it++) {
        GeneTree* tree = (GeneTree*)(*it);
        int time = 0;
        vector<pair<int, int>> eulerInternalBranch;
        map<string, int> leafIndex;

        tree->setRootLeaf(NULL);
        assert(tree->root->isLeaf());

        leafIndex[tree->root->name] = 0;
        dfsMRP(tree->root->neighbors[0]->node, tree->root, time, eulerInternalBranch, leafIndex);

        for (int i = 0; i < seqNames.size(); ++i) {
            if (leafIndex.find(seqNames[i]) != leafIndex.end()) {

                int leafInd = leafIndex[seqNames[i]];
                for (auto [in, out] : eulerInternalBranch) {
                    sequences[i] += (leafInd >= in && leafInd <= out) ? '1' : '0';
                }

            } else {
                sequences[i] += string(eulerInternalBranch.size(), '?');
            }
        }
    }

    mrpAln = new Alignment(seqNames, sequences, params->sequence_type);
    mrpAln->printPhylip(cout);
}

void PhyloSuperTreeUnlinked::doMRP() {
    buildMRPMatrix();
    mrpTree = new GeneTree(mrpAln);
    mrpTree->treeParams = *(this->params);
    runOptimizeAndReconstruction(mrpTree->treeParams, mrpTree);
    // TODO: get best tree or greedy consensus tree or random tree
    switch (params->mrp_type) {
        case MRPType::MRP_GREEDY:
            break;
        case MRPType::MRP_RANDOM:
            break;
        case MRPType::MRP_BEST:
            break;
        default:
            break;
    }
}

void PhyloSuperTreeUnlinked::printResultWithMRPTree() {
    assert(conAln);
    assert(mrpTree);

    IQTree *tmpTree = new IQTree(conAln);
    tmpTree->copyTree(mrpTree);

    Params newParams = *(this->params);
    tmpTree->setParams(newParams);

    cout << "\nSCORE OF MRP TREE: " << tmpTree->computeParsimony() << endl;
    tmpTree->printResultTree();
    
    delete tmpTree;
}

void PhyloSuperTreeUnlinked::doSCM() {
    StrVector sourcesTree;
    for (auto it = begin(); it != end(); it++) {
        GeneTree* tree = (GeneTree*)(*it);
        sourcesTree.push_back(tree->getTreeString());
    }
    getAllSeqNames();
    StrictConsensusMerge scm(sourcesTree, seqNameToIndex);
    scm.run();
}

void PhyloSuperTreeUnlinked::printSCMTree() {
}