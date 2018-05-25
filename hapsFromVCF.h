//
//  hapsFromVCF.h
//  fineRADstructure
//
//  Created by Milan Malinsky on 05/03/2018.
//  Copyright © 2018 Milan Malinsky. All rights reserved.
//

#ifndef hapsFromVCF_h
#define hapsFromVCF_h

#include <stdio.h>

int VCFhapsMain(int argc, char** argv);
void parseVCFoptions(int argc, char** argv);

inline bool isDNAonly(char base) {
    bool DNAonly;
    switch (base)
    {
        case 'A': DNAonly = true; break;
        case 'C': DNAonly = true; break;
        case 'G': DNAonly = true; break;
        case 'T': DNAonly = true; break;
        default: DNAonly = false;
    }
    return DNAonly;
}


#endif /* hapsFromVCF_h */
