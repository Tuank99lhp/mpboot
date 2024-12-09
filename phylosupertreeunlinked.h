//
//  phylosupertreeunlinked.h
//  tree
//
//  Created by Minh Bui on 2/5/18.
//

#ifndef phylosupertreeunlinked_h
#define phylosupertreeunlinked_h

#include "phylosupertree.h"

/**
    Super-tree with separate partition tree topologies
 */
class PhyloSuperTreeUnlinked : public PhyloSuperTree {
public:
    /**
     constructors
     */
    PhyloSuperTreeUnlinked(Params &params);

    virtual bool isSuperTreeUnlinked() {
        return true;
    }

};

#endif /* phylosupertreeunlinked_h */
