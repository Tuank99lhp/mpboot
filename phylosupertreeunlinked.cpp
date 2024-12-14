//
//  phylosupertreeunlinked.cpp
//  tree
//
//  Created by Minh Bui on 2/5/18.
//

#include "phylosupertreeunlinked.h"
#include "iqtree.h"
#include "timeutil.h"

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
    for (int i = 0; i < size(); i++)
        at(i)->setAlignment(saln->partitions[i]);
}

void PhyloSuperTreeUnlinked::setParams(Params& params) {
    PhyloSuperTree::setParams(params);
    for (auto it = begin(); it != end(); it++)
        ((IQTree*)(*it))->setParams(params);
}

void PhyloSuperTreeUnlinked::mapTrees() {
    // do nothing here as partition trees are unlinked
}

void PhyloSuperTreeUnlinked::computeParsimonyTree(const char *out_prefix, Alignment *alignment) {
    assert(0);
}

bool PhyloSuperTreeUnlinked::isBifurcating(Node *node, Node *dad) {
    for (auto it = begin(); it != end(); it++)
        if (!(*it)->isBifurcating())
            return false;
    return true;
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
        (*it)->readTree(str, rooted);
        (*it)->assignLeafNames();
        // (*it)->resetCurScore();
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

double PhyloSuperTreeUnlinked::treeLength(Node *node, Node *dad) {
    double len = 0.0;
    for (iterator tree = begin(); tree != end(); tree++)
        len += (*tree)->treeLength();
    return len;
}

double PhyloSuperTreeUnlinked::treeLengthInternal(double epsilon, Node *node, Node *dad) {
    double len = 0.0;
    for (iterator tree = begin(); tree != end(); tree++)
        len += (*tree)->treeLengthInternal(epsilon);
    return len;
}

string PhyloSuperTreeUnlinked::doNNISearch(int &nniCount, int &nniSteps) {
    for (iterator tree = begin(); tree != end(); tree++)
        ((IQTree*)(*tree))->doNNISearch(nniCount, nniSteps);
}

double PhyloSuperTreeUnlinked::doTreeSearch() {
    double score = 0.0;
    for (iterator tree = begin(); tree != end(); tree++)
        score += ((IQTree*)(*tree))->doTreeSearch();
    return score;
}

void PhyloSuperTreeUnlinked::summarizeBootstrap(Params &params) {
    assert(0);
    // for (auto tree = begin(); tree != end(); tree++)
    //     ((IQTree*)*tree)->summarizeBootstrap(params);
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