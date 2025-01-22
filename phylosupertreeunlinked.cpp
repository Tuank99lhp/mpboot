//
//  phylosupertreeunlinked.cpp
//  tree
//
//  Created by Minh Bui on 2/5/18.
//

#include "phylosupertreeunlinked.h"
#include "iqtree.h"
#include "timeutil.h"
#include "sprparsimony.h"

extern Params *globalParam;

PhyloSuperTreeUnlinked::PhyloSuperTreeUnlinked(Params &params): PhyloSuperTree(params, true) {

}

void PhyloSuperTreeUnlinked::readTree(istream &in, bool &is_rooted) {
    for (iterator it = begin(); it != end(); it++) {
        (*it)->rooted = globalParam->is_rooted;
        (*it)->readTree(in, (*it)->rooted);
        is_rooted |= (*it)->rooted;
    }
}

void PhyloSuperTreeUnlinked::setAlignment(Alignment *alignment) {
    assert(alignment->isSuperAlignment());
    SuperAlignment *saln = (SuperAlignment*)alignment;
    assert(saln->partitions.size() == size());
    for (int i = 0; i < size(); i++) {
        at(i)->setAlignment(saln->partitions[i]);
    }
}

void PhyloSuperTreeUnlinked::setParams(Params& params) {
    PhyloSuperTree::setParams(params);
    for (auto it = begin(); it != end(); it++) {
        ((IQTree*)(*it))->setParams(params);
        ((IQTree*)(*it))->curScore = -DBL_MAX;
    }
}

void PhyloSuperTreeUnlinked::deleteAllPartialLh() {
    for (auto it = begin(); it != end(); it++)
        (*it)->deleteAllPartialLh();
}

void PhyloSuperTreeUnlinked::initializePLL(Params &params) {
    for (auto it = begin(); it != end(); it++)
        ((IQTree*)(*it))->initializePLL(params);
}

void PhyloSuperTreeUnlinked::mapTrees() {
    // do nothing here as partition trees are unlinked
}

void PhyloSuperTreeUnlinked::initializeAllPartialPars() {
    for (auto it = begin(); it != end(); it++)
        (*it)->initializeAllPartialPars();
}

void PhyloSuperTreeUnlinked::clearAllPartialLH() {
    for (auto it = begin(); it != end(); it++)
        (*it)->clearAllPartialLH();
}

int PhyloSuperTreeUnlinked::computeParsimony() {
    int score = 0;
    for (auto it = begin(); it != end(); it++) {
        (*it)->curScore = -(*it)->computeParsimony();
        score += -(*it)->curScore;
    }
    return score;
}

double PhyloSuperTreeUnlinked::getBestScore() {
    double score = 0.0;
    for (auto it = begin(); it != end(); it++) {
        IQTree *iqtree = (IQTree*)(*it);
        score += iqtree->getBestScore();
    }
    return score;
}

void PhyloSuperTreeUnlinked::printBestScores(int numBestScore) {
    for (auto it = begin(); it != end(); it++) {
        IQTree *iqtree = (IQTree*)(*it);
        cout << "Tree " << it - begin() << ": ";
        iqtree->printBestScores(iqtree->candidateTrees.popSize);
    }
}

void PhyloSuperTreeUnlinked::computeParsimonyTree(const char *out_prefix, Alignment *alignment) {
    SuperAlignment *saln = (SuperAlignment*)alignment;
    assert(saln->partitions.size() == size());
    for (int i = 0; i < size(); i++) {
        at(i)->computeParsimonyTree(NULL, saln->partitions[i]);
    }
    if (out_prefix) {
        string file_name = out_prefix;
        file_name += ".parstree";
        try {
            ofstream out;
            out.open(file_name.c_str());
            for (int i = 0; i < size(); i++) {
                at(i)->printTree(out, WT_NEWLINE);
            }
            out.close();
        } catch (...) {
            outError("Cannot write to file ", file_name);
        }
    }
}

string PhyloSuperTreeUnlinked::getTreeString() {
    stringstream tree_stream;
    for (iterator it = begin(); it != end(); it++)
        (*it)->printTree(tree_stream, WT_TAXON_ID + WT_BR_LEN + WT_SORT_TAXA);
    return tree_stream.str();
}

void PhyloSuperTreeUnlinked::readTreeString(const string &tree_string) {
    stringstream str;
    str << tree_string;
    str.seekg(0, ios::beg);
    for (iterator it = begin(); it != end(); it++) {
        (*it)->freeNode();
        (*it)->readTree(str, (*it)->rooted);
        (*it)->assignLeafNames();
        (*it)->clearAllPartialLH();
    }
}

void PhyloSuperTreeUnlinked::setRootNode(char *my_root) {
    // DOES NOTHING
}

void PhyloSuperTreeUnlinked::computeBranchLengths() {
    // DOES NOTHING
}

void PhyloSuperTreeUnlinked::printTree(ostream &out, int brtype) {
    for (iterator tree = begin(); tree != end(); tree++)
        (*tree)->printTree(out, brtype);
}

void PhyloSuperTreeUnlinked::printResultTree(string suffix) {
    string tree_file_name = params->out_prefix;
    tree_file_name += ".treefile";
    if (suffix.compare("") != 0) {
        tree_file_name += "." + suffix;
    }
    ofstream out;
    out.open(tree_file_name.c_str());
    for (iterator tree = begin(); tree != end(); tree++)
        (*tree)->printTree(out, WT_BR_LEN | WT_BR_LEN_FIXED_WIDTH | WT_SORT_TAXA | WT_NEWLINE);
    out.close();
    if (verbose_mode >= VB_MED)
        cout << "Best tree printed to " << tree_file_name << endl;
}

string PhyloSuperTreeUnlinked::doNNISearch(int &nniCount, int &nniSteps) {
    assert(0);
}

double PhyloSuperTreeUnlinked::doTreeSearch() {
    double score = 0.0;

    VerboseMode saved_mode = verbose_mode;
    // verbose_mode = VB_QUIET;
    bool saved_print_ufboot_trees = params->print_ufboot_trees;
    params->print_ufboot_trees = false;
    
    for (iterator tree = begin(); tree != end(); tree++) {
        IQTree *iqtree = (IQTree*)(*tree);
        iqtree->initializePLLSankoff(*iqtree->params);
        score += iqtree->doTreeSearch();
        cout << "Tree " << tree - begin() << " score: " << -iqtree->bestScore << endl;
    }

    verbose_mode = saved_mode;
    params->print_ufboot_trees = saved_print_ufboot_trees;

    if (score > bestScore) {
        bestScore = score;
        bestTreeString = getTreeString();
    }
    curScore = score;
    
    // candidateTrees.update(getTreeString(), curScore);
    printResultTree();
    return bestScore;
}

void PhyloSuperTreeUnlinked::summarizeBootstrap(Params &params) {
    for (auto tree = begin(); tree != end(); tree++) {
        ((IQTree*)*tree)->summarizeBootstrap(params);
    }
}

void PhyloSuperTreeUnlinked::writeUFBootTrees(Params &params, StrVector &removed_seqs, StrVector &twin_seqs) {
    assert(0);
    // //    IntVector tree_weights;
    // int i, j;
    // string filename = params.out_prefix;
    // filename += ".ufboot";
    // ofstream out(filename.c_str());
    
    // for (auto tree = begin(); tree != end(); tree++) {
    //     MTreeSet trees;

    //     trees.init(((IQTree*)*tree)->boot_trees, (*tree)->rooted);
    //     for (i = 0; i < trees.size(); i++) {
    //         NodeVector taxa;
    //         // change the taxa name from ID to real name
    //         trees[i]->getOrderedTaxa(taxa);
    //         for (j = 0; j < taxa.size(); j++)
    //             taxa[j]->name = aln->getSeqName(taxa[j]->id);
    //         if (removed_seqs.size() > 0) {
    //             // reinsert removed seqs into each tree
    //             trees[i]->insertTaxa(removed_seqs, twin_seqs);
    //         }
    //         // now print to file
    //         for (j = 0; j < trees.tree_weights[i]; j++)
    //             if (params.print_ufboot_trees == 1)
    //                 trees[i]->printTree(out, WT_NEWLINE);
    //             else
    //                 trees[i]->printTree(out, WT_NEWLINE + WT_BR_LEN);
    //     }
    // }
    // cout << "UFBoot trees printed to " << filename << endl;
    // out.close();
}

void PhyloSuperTreeUnlinked::optimizeAlignments(Params &params) {
    for (auto tree = begin(); tree != end(); tree++) {
        IQTree* iqtree = (IQTree*)*tree;
        optimizeAlignment(iqtree, params);
    }
}

void PhyloSuperTreeUnlinked::computeInitialPLLTrees(Params &params) {
    for (auto tree = begin(); tree != end(); tree++) {
        IQTree* iqtree = (IQTree*)*tree;
        iqtree->pllInst->randomNumberSeed = params.ran_seed;
        if(params.maximum_parsimony){
            _pllComputeRandomizedStepwiseAdditionParsimonyTree(iqtree->pllInst, iqtree->pllPartitions, params.sprDist, iqtree);
        }
        else
            pllComputeRandomizedStepwiseAdditionParsimonyTree(iqtree->pllInst, iqtree->pllPartitions, params.sprDist);

        resetBranches(iqtree->pllInst);
        pllTreeToNewick(iqtree->pllInst->tree_string, iqtree->pllInst, iqtree->pllPartitions, iqtree->pllInst->start->back,
                PLL_TRUE, PLL_TRUE, PLL_FALSE, PLL_FALSE, PLL_FALSE, PLL_SUMMARIZE_LH, PLL_FALSE, PLL_FALSE);
        iqtree->readTreeString(string(iqtree->pllInst->tree_string));
        iqtree->initializeAllPartialPars();
        iqtree->clearAllPartialLH();
    }
}

void PhyloSuperTreeUnlinked::initCandidateUnlinkedTreeSet(Params &params, int numInitTrees) {
    for (auto tree = begin(); tree != end(); tree++) {
        IQTree* iqtree = (IQTree*)*tree;
        int numDup = initCandidateTreeSet(params, *iqtree, numInitTrees);
        assert(iqtree->candidateTrees.size() != 0);
        // cout << "Finish initializing candidate tree set. ";
        // cout << "Number of distinct locally optimal trees: " << iqtree->candidateTrees.size() << endl;
    }
    bestScore = getBestScore();
}

void PhyloSuperTreeUnlinked::reportUnlinkedPhyloAnalysis(Params &params, string &original_model, vector<ModelInfo> &model_info, StrVector &removed_seqs, StrVector &twin_seqs) {
    for (auto tree = begin(); tree != end(); tree++) {
        IQTree* iqtree = (IQTree*)*tree;
        reportPhyloAnalysis(params, original_model, *(iqtree->aln), *iqtree, model_info, removed_seqs, twin_seqs);
    }
}