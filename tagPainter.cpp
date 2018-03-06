//
//  main.cpp
//  RADTAGpainter
//
//  Created by Milan Malinsky on 29/01/2016.
//  Copyright (c) 2016 Milan Malinsky. All rights reserved.
//

#include "paintSql.h"
#include "hapsFromVCF.h"
#include "utils.h"

 

static const char *VERSION_MESSAGE =
"RADPainter Version " V "\n"
"Written by Milan Malinsky.\n"
"\n";


static const char *USAGE_MESSAGE =
"Program: " BIN "\n"
"Version: " V "\n"
"Contact: " THIS_AUTHOR " [" BUGREPORT "]\n"
"Usage: " BIN " <command> [options]\n\n"
"Commands:\n"
"           paint           Get co-ancestry matrix for fineStructure\n"
"           hapsFromVCF     (experimental) Get the haplotype format input from a VCF file\n"
"                           !NOTE! for now this only works for some VCF formats\n"
"\nReport bugs to " BUGREPORT "\n\n";


int main(int argc, char * argv[]) {
    // insert code here...
    if(argc <= 1)
    {
        std::cout << USAGE_MESSAGE;
        return 0;
    }
    else
    {
        std::string command(argv[1]);
        if(command == "help" || command == "--help" || command == "-h")
        {
            std::cout << USAGE_MESSAGE;
            return 0;
        }
        else if(command == "version" || command == "--version")
        {
            std::cout << VERSION_MESSAGE;
            return 0;
        }
        
        if(command == "paint")
            paintSqlMain(argc - 1, argv + 1);
        else if(command == "hapsFromVCF")
            VCFhapsMain(argc - 1, argv + 1);
        else
        {
            std::cerr << "Unrecognized command: " << command << "\n";
            return 1;
        }
        return 0;
    }
}
