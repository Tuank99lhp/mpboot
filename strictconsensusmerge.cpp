#include "strictconsensusmerge.h"

StrictConsensusMerge::StrictConsensusMerge(const StrVector &sourceTrees, const map<string, int> &seqNameToIndex) {
    for (int i = 0; i < sourceTrees.size(); ++i) {
        GeneTree *gt = new GeneTree(sourceTrees[i]);

        if (gt->leafNum < 4) {
            cout << "SCM: Tree " << i << " uninformative, less than 4 leaves" << endl;
            delete gt;
            continue; 
        }
        
        gt->setOriginalNodeIndex(seqNameToIndex);
        trees.push_back(gt);
    }

    cout << "SCM: Init SCM with " << trees.size() << " trees" << endl;
}

StrictConsensusMerge::~StrictConsensusMerge() {
    
}

int StrictConsensusMerge::getOverlap(GeneTree *tree1, GeneTree *tree2) {
    Split ans = *tree1->getLeavesSplit();
    ans *= *tree2->getLeavesSplit();
    return ans.countTaxa();
}

void StrictConsensusMerge::rerootOnLowestCommonIndexPath(GeneTree *tree, Split commonMask) {
    Split l = commonMask.lowestBitOnly();
    assert(commonMask.countTaxa() > 2);

    assert(tree->splitToNode.find(l) != tree->splitToNode.end());
    GeneNode *curNode = tree->splitToNode[l];
    assert(curNode->isLeaf());

    GeneNode *p = curNode->parent;
    while (p != NULL) {
        if ((*p->split * commonMask) != l) {
            break;
        }
        curNode = p;
        p = p->parent;
    }

    assert(p != NULL);
    assert((*curNode->split * commonMask) == l);

    Split withoutLowest = commonMask - l;
    tree->makeFirstChildOfRoot(curNode);
    
    while (true) {
        bool tmp = false;
        for (auto rootChild: tree->seedNode->getChildren()) {
            if (rootChild == curNode) {
                assert(tmp == false);
                tmp = true;
                continue;
            }

            assert(tmp);

            Split cm = *rootChild->split * withoutLowest;
            if (cm.countTaxa() > 0) {                
                if (cm == withoutLowest) {
                    
                    GeneNode *oldRoot = tree->seedNode;
                    assert(rootChild->parent == oldRoot);

                    tree->reroot(rootChild);
                    assert(oldRoot->parent == rootChild);
                    
                    tree->makeFirstChildOfRoot(oldRoot);
                    curNode = oldRoot;

                    break;
                } else {
                    return;
                }
            }
        }
        assert(tmp);
    }
}

void StrictConsensusMerge::collapsePathNotFound(
    map<Split, vector<pair<Split, GeneNode*>>> &treeRelevantSplits, 
    map<Split, vector<pair<Split, GeneNode*>>> &otherTreeRelevantSplits,
    GeneTree *tree
) {
    vector<Split> delMasked;
    for (const auto& [masked, path]: treeRelevantSplits) {
        if (otherTreeRelevantSplits.find(masked) == otherTreeRelevantSplits.end()) {
            for (auto it: path) {
                GeneNode *node = it.second;
                assert(node->parent != NULL);
                tree->collapseEdge(node->parent, node);
            }
            delMasked.push_back(masked);
        }
    }
    for (auto it: delMasked) {
        assert(treeRelevantSplits.erase(it) == 1);
    }
}

// Merge tree2 to tree1, tree2 should be deleted after the merge
void StrictConsensusMerge::pairwiseMerger(GeneTree *tree1, GeneTree *tree2) {
    assert(tree1->seedNode->degree() > 2 && tree2->seedNode->degree() > 2);
    assert(tree1->splitToNode.empty() == false && tree2->splitToNode.empty() == false);

    Split leavesIntersection = *tree1->getLeavesSplit();
    leavesIntersection *= *tree2->getLeavesSplit();

    int numCommonLeaves = leavesIntersection.countTaxa();
    if (numCommonLeaves < 2) {
        cout << "SCM: Trees must have at least 2 common leaves" << endl;
        assert(0);
    }

    if (numCommonLeaves == 2) {
        cout << "SCM with 2 leaves in common results in a polytomy" << endl;
        tree1->collapseSubtree(tree1->seedNode, false);
        tree2->collapseSubtree(tree2->seedNode, false);

        for (auto child: tree2->seedNode->getChildren()) {
            if (leavesIntersection.containTaxon(child->getOriginalIndex())) {
                continue;
            }

            tree2->seedNode->removeChild(child);
            tree1->seedNode->addChild(child);
        }
        tree1->resetAfterMergeSCM();
        return;
    }

    rerootOnLowestCommonIndexPath(tree1, leavesIntersection);
    rerootOnLowestCommonIndexPath(tree2, leavesIntersection);

    map<Split, vector<pair<Split, GeneNode*>>> tree1RelevantSplits;
    map<Split, vector<pair<Split, GeneNode*>>> tree2RelevantSplits;

    for (const auto& [s, e]: tree1->splitToNode) {
        assert(s == *e->split);
        Split masked = leavesIntersection * (*e->split);
        if (masked.countTaxa() > 0 && masked != leavesIntersection) {
            if (tree1RelevantSplits.find(masked) == tree1RelevantSplits.end()) {
                tree1RelevantSplits[masked] = vector<pair<Split, GeneNode*>>();
            }
            tree1RelevantSplits[masked].push_back(make_pair(s, e));
        }
    }

    for (const auto& [s, e]: tree2->splitToNode) {
        assert(s == *e->split);
        Split masked = leavesIntersection * (*e->split);
        if (masked.countTaxa() > 0 && masked != leavesIntersection) {
            if (tree2RelevantSplits.find(masked) == tree2RelevantSplits.end()) {
                tree2RelevantSplits[masked] = vector<pair<Split, GeneNode*>>();
            }
            tree2RelevantSplits[masked].push_back(make_pair(s, e));
        }
    }

    // NOTE: After this point, the trees's splits will be unchanged and can be unmatch the trees's topology

    for (auto& [_, path]: tree1RelevantSplits) {
        sort(path.begin(), path.end(), [](const pair<Split, GeneNode*> &a, const pair<Split, GeneNode*> &b) {
            return a.first.countTaxa() > b.first.countTaxa();
        });
        for (int i = 1; i < path.size(); ++i) {
            assert(path[i].second->parent == path[i - 1].second);
        }
    }

    for (auto& [_, path]: tree2RelevantSplits) {
        sort(path.begin(), path.end(), [](const pair<Split, GeneNode*> &a, const pair<Split, GeneNode*> &b) {
            return a.first.countTaxa() > b.first.countTaxa();
        });
        for (int i = 1; i < path.size(); ++i) {
            assert(path[i].second->parent == path[i - 1].second);
        }
    }

    collapsePathNotFound(tree1RelevantSplits, tree2RelevantSplits, tree1);
    collapsePathNotFound(tree2RelevantSplits, tree1RelevantSplits, tree2);

    for (auto child: tree2->seedNode->getChildren()) {
        if (child->split->overlap(leavesIntersection)) {
            continue;
        }
        tree2->seedNode->removeChild(child);
        tree1->seedNode->addChild(child);
    }

    for (const auto& [tree2Masked, tree2Path]: tree2RelevantSplits) {
        assert(tree1RelevantSplits.find(tree2Masked) != tree1RelevantSplits.end());
        auto tree1Path = tree1RelevantSplits[tree2Masked];

        GeneNode *tree1Head = tree1Path.back().second;
        GeneNode *tree2Head = tree2Path.back().second;

        for (auto child: tree2Head->getChildren()) {
            if (child->split->overlap(leavesIntersection)) {
                continue;
            }
            tree2Head->removeChild(child);
            tree1Head->addChild(child);
        }

        if (tree2Path.size() > 1) {
            if (tree1Path.size() > 1) {
                for (int i = 1; i < tree1Path.size() - 1; ++i) {
                    GeneNode *node = tree1Path[i].second;
                    assert(node->parent == tree1Path[0].second);
                    tree1->collapseEdge(node->parent, node);
                }
                
                GeneNode *midNode = tree1Path[0].second;
                
                for (int i = 1; i < tree2Path.size(); ++i) {
                    GeneNode *p = tree2Path[i - 1].second;
                    GeneNode *avoid = tree2Path[i].second;
                    assert(avoid->parent == p);

                    for (auto child: p->getChildren()) {
                        if (child == avoid) {
                            continue;
                        }
                        assert(child->split->overlap(leavesIntersection) == false);
                        p->removeChild(child);
                        midNode->addChild(child);
                    }
                }
            } else {
                GeneNode *tree1Head = tree1Path[0].second;
                GeneNode *tree1Tail = tree1Head->parent;

                GeneNode *tree2DeepestNode = tree2Path[0].second;
                GeneNode *tree2TipmostNode = tree2Path.back().second->parent;
                GeneNode *prevHead = tree2Path.back().second;

                tree2TipmostNode->removeChild(prevHead);
                
                tree2DeepestNode->parent->removeChild(tree2DeepestNode);
                tree1Tail->addChild(tree2DeepestNode);

                tree1Tail->removeChild(tree1Head);
                tree2TipmostNode->addChild(tree1Head);
            }
        }
    }

    tree1->resetAfterMergeSCM();
}

GeneTree *StrictConsensusMerge::getSCMTree() {
    if (trees.size() == 0) {
        cout << "SCM: No trees to merge" << endl;
        return NULL;
    }
    if (trees.size() == 1) {
        cout << "SCM: Only one tree, no need to merge" << endl;
        return trees[0];
    }
    cout << "SCM: Running SCM..." << endl;

    int nTrees = trees.size();
    vector<vector<pair<int, int>>> treePairsScore(nTrees, vector<pair<int, int>>());
    
    for (int i = 0; i < nTrees; ++i) {
        for (int j = i + 1; j < nTrees; ++j) {
            treePairsScore[i].push_back(make_pair(getOverlap(trees[i], trees[j]), j));
        }
        sort(treePairsScore[i].begin(), treePairsScore[i].end());
    }

    int numMergers = trees.size() - 1;
    for (int iter = 0; iter < numMergers; ++iter) {
        int tree1Index, tree2Index;
        tree1Index = tree2Index = -1;

        int maxIntersection = 3;
        for (int i = 0; i < treePairsScore.size(); ++i) {
            while (!treePairsScore[i].empty()) {
                int j = treePairsScore[i].back().second;
                if (trees[j] != NULL) {
                    break;
                }
                treePairsScore[i].pop_back();
            }
            if (!treePairsScore[i].empty()) {
                int intersection = treePairsScore[i].back().first;
                if (intersection > maxIntersection) {
                    maxIntersection = intersection;
                    tree1Index = i;
                    tree2Index = treePairsScore[i].back().second;
                }
            }
        }

        if (maxIntersection < 4) {
            cout << "SCM: Insufficient overlap for merger" << endl;
            assert(0);
        }

        assert(tree1Index != -1 && tree2Index != -1);

        pairwiseMerger(trees[tree1Index], trees[tree2Index]);

        // Set root as seedNode for deleting unnecessary node
        trees[tree2Index]->root = trees[tree2Index]->seedNode;
        delete trees[tree2Index];
        
        trees.push_back(trees[tree1Index]);
        trees[tree1Index] = trees[tree2Index] = NULL;
        
        treePairsScore.push_back(vector<pair<int, int>>());
        treePairsScore[tree1Index].clear();
        treePairsScore[tree2Index].clear();
        
        for (int j = 0; j < trees.size() - 1; ++j) {
            if (trees[j] == NULL) {
                continue;
            }
            treePairsScore.back().push_back(make_pair(getOverlap(trees.back(), trees[j]), j));
        }
        sort(treePairsScore.back().begin(), treePairsScore.back().end());
    }

    cout << "SCM: Merged trees successfully" << endl;
    // for (int i = 0; i < trees.size() - 1; ++i) {
    //     assert(trees[i] == NULL);
    // }
    // assert(trees.back() != NULL);
    return trees.back();
}